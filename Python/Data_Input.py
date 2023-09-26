# Yuncong Ma, 9/18/2023
# Data Input module of pNet


#########################################
# Packages
import nibabel as nib
import numpy as np
import scipy.io
import scipy.io as sio
import os
import re
import json
import torch
import h5py
import time


def load_matlab_array(file_matlab: str, variable_name: str):
    """
    Load a single matlab variable with variable name into np array
    This support both matrix and cell vector

    :param file_matlab: string
    :param variable_name: string
    :return: data as nparray\

    By Yuncong Ma, 9/6/2023
    """
    try:
        matlab_data = sio.loadmat(file_matlab, variable_names=variable_name)
        data = np.array(matlab_data[variable_name])
    except NotImplementedError:
        matlab_data = h5py.File(file_matlab, 'r')
        data = np.array(matlab_data[variable_name]).T
    finally:
        ValueError('Cannot read the MATLAB file: '+str(file_matlab))

    return data


def load_matlab_single_array(file_matlab: str):
    """
    Load a matlab file with only one variable stored
    This support both matrix and cell vector

    :param file_matlab: string
    :return: data as nparray

    By Yuncong Ma, 9/6/2023
    """
    version = 0
    try:
        matlab_data = sio.loadmat(file_matlab)
    except NotImplementedError:
        matlab_data = h5py.File(file_matlab, 'r')
        version = 7.3
    finally:
        ValueError('Cannot read the MATLAB file: '+str(file_matlab))

    variable_names = matlab_data.keys()
    actual_variable_names = [name for name in variable_names if not name.startswith('__')]
    # In case there are more than one variable in the matlab file
    if len(actual_variable_names) > 1:
        print('The MATLAB file ' + file_matlab + ' contains more than one variable')
        print('This file contains ' + ', '.join(actual_variable_names))
        data = []
        return
    # Extract the content in the variable
    data = np.array(matlab_data[actual_variable_names[0]])
    if version == 7.3:
        data = data.T
    return data


def load_matlab_single_variable(file_matlab: str):
    """
    Load a matlab file with only one variable stored
    This support both matrix and cell vector

    :param file_matlab: string
    :return: data as its original format

    By Yuncong Ma, 9/24/2023
    """
    version = 0
    try:
        matlab_data = sio.loadmat(file_matlab)
    except NotImplementedError:
        version = 7.3
        matlab_data = h5py.File(file_matlab, 'r')
    finally:
        ValueError('Cannot read the MATLAB file: '+str(file_matlab))

    variable_names = matlab_data.keys()
    actual_variable_names = [name for name in variable_names if not name.startswith('__')]
    # In case there are more than one variables in the matlab file
    if len(actual_variable_names) > 1:
        print('The MATLAB file ' + file_matlab + ' contains more than one variable')
        print('This file contains ' + ', '.join(actual_variable_names))
        data = []
        return data
    # Extract the content in the variable
    data = matlab_data[actual_variable_names[0]]
    return data


def set_data_precision(data_precision: str):
    """
    Set the data format and eps
    Support single, float32, double, float64

    :param data_precision: 'float32' or 'float64' in python, 'single' or 'double' in MATLAB
    :return: np_float, np_eps

    By Yuncong Ma, 9/6/2023
    """
    if data_precision.lower() == 'single' or data_precision.lower() == 'float32':
        np_float = np.float32
        np_eps = np.finfo(np_float).eps
    elif data_precision.lower() == 'double' or data_precision.lower() == 'float64':
        np_float = np.float64
        np_eps = np.finfo(np_float).eps
    else:
        raise ValueError('Unknown data type settings: ' + data_precision)
    return np_float, np_eps


def set_data_precision_torch(data_precision: str):
    """
    Set the data format and eps
    Support single, float32, double, float64

    :param data_precision: 'float32' or 'float64' in python, 'single' or 'double' in MATLAB
    :return: torch_float, torch_eps

    By Yuncong Ma, 9/6/2023
    """
    if data_precision.lower() == 'single' or data_precision.lower() == 'torch.float32':
        torch_float = torch.float32
        torch_eps = torch.finfo(torch_float).eps
    elif data_precision.lower() == 'double' or data_precision.lower() == 'torch.float64':
        torch_float = torch.float64
        torch_eps = torch.finfo(torch_float).eps
    else:
        raise ValueError('Unknown data type settings: ' + data_precision)
    torch_eps = torch.tensor(torch_eps)
    return torch_float, torch_eps


def write_json_setting(setting, file_setting: str):
    """
    Write setting parameter in json format

    :param setting: a json based variable
    :param file_setting: Directory of a json setting file
    :return: none

    By Yuncong Ma, 9/26/2023
    """
    file_extension = os.path.splitext(file_setting)[1]

    if file_extension != '.json':
        raise ValueError('It is not a JSON file: '+file_setting)
    # save serialized json file
    with open(file_setting, 'w') as file:
        json.dump(setting, file, indent=4)


def load_json_setting(file_setting: str):
    """
    Load setting variable in a json file

    :param file_setting: Directory of a json setting file
    :return: Setting

    By Yuncong Ma, 9/6/2023
    """
    file_extension = os.path.splitext(file_setting)[1]

    if file_extension != '.json':
        raise ValueError('It is not a JSON file: '+file_setting)

    with open(file_setting, 'r') as file:
        json_string = file.read()
    Setting = json.loads(json_string)
    return Setting


def normalize_data(data, algorithm='vp', normalization='vmax'):
    """
    Normalize data by algorithm and normalization settings

    :param data: data in 2D matrix [dim_time, dim_space]
    :param algorithm: 'z' 'gp' 'vp'
    :param normalization: 'n2' 'n1' 'rn1' 'g' 'vmax'
    :return: data
    Consistent to MATLAB function normalize_data(X, algorithm, normalization)
    'vp' is to shift each vector to all non-negative
    'vmax' is to normalize each vector by its max value

    By Yuncong Ma, 9/6/2023
    """

    if len(data.shape) != 2:
        raise ValueError("data must be a 2D matrix")

    X = np.array(data)
    np_float, np_eps = set_data_precision(str(X.dtype))

    if algorithm.lower() == 'z':
        # standard score for each variable
        mVec = np.mean(X, axis=1)
        sVec = np.maximum(np.std(X, axis=1), np_eps)
        pX = (X - mVec[:, np.newaxis]) / sVec[:, np.newaxis]
    elif algorithm.lower() == 'gp':
        # remove negative value globally
        minVal = np.min(X)
        shiftVal = np.abs(np.minimum(minVal, 0))
        pX = X + shiftVal
    elif algorithm.lower() == 'vp':
        # remove negative value voxel-wisely
        minVal = np.min(X, axis=0)
        shiftVal = np.abs(np.minimum(minVal, 0))
        pX = X + np.tile(shiftVal,(X.shape[0],1))
    else:
        # do nothing
        print('  unknown preprocess parameters, no preprocess applied')
        pX = X

    if normalization.lower() == 'n2':
        # l2 normalization for each observation
        l2norm = np.sqrt(np.sum(pX ** 2, axis=1)) + np_eps
        pX = pX / l2norm[:, np.newaxis]
    elif normalization.lower() == 'n1':
        # l1 normalization for each observation
        l1norm = np.sum(pX, axis=1) + np_eps
        pX = pX / l1norm[:, np.newaxis]
    elif normalization.lower() == 'rn1':
        # l1 normalization for each variable
        l1norm = np.sum(pX, axis=0) + np_eps
        pX = pX / l1norm
    elif normalization.lower() == 'g':
        # global scale
        sVal = np.sort(pX, axis=None)
        perT = 0.001
        minVal = sVal[int(len(sVal) * perT)]
        maxVal = sVal[int(len(sVal) * (1 - perT))]
        pX[pX < minVal] = minVal
        pX[pX > maxVal] = maxVal
        pX = (pX - minVal) / max((maxVal - minVal), np_eps)
    elif normalization.lower() == 'vmax':
        cmin = np.tile(np.min(pX, axis=0), (X.shape[0], 1))
        cmax = np.tile(np.max(pX, axis=0), (X.shape[0], 1))
        pX = (pX - cmin) / np.maximum(cmax - cmin, np_eps)
    else:
        # do nothing
        pX = X
        print('  unknown normalization parameters, no normalization applied')

    if np.isnan(pX).any():
        raise ValueError('  nan exists, check the preprocessed data')

    return pX


def load_fmri_scan(file_scan_list: str, dataType: str, dataFormat: str, Reshape=False,
                   Brain_Mask=None, Normalization=None, Concatenation=True, logFile=None):
    """
    Load one or multiple fMRI scans, and concatenate them into a single 2D matrix along the time dimension
    Optional normalization can be added for each scan before concatenation

    :param file_scan_list: Directory of a single txt file storing fMRI file directories, or a directory of a single scan file
    :param dataType: 'Surface', 'Volume'
    :param dataFormat: 'HCP Surface (*.cifti, *.mat)', 'MGH Surface (*.mgh)', 'MGZ Surface (*.mgz)', 'Volume (*.nii, *.nii.gz, *.mat)'
    :param Reshape: False or True, whether to reshape 4D volume-based fMRI data to 2D
    :param Brain_Mask: None or a brain mask [X Y Z]
    :param Normalization: False, 'vp-vmax'
    :param Concatenation: True, False
    :param logFile: a log file to save the output
    :return: Data: a 2D or 4D NumPy array [dim_time dim_space]

    By Yuncong Ma, 9/25/2023
    """

    # Check setting
    check_data_type_format(dataType, dataFormat, logFile=logFile)

    # Suppress warning messages when loading CIFTI 2 formatted files
    nib.imageglobals.logger.setLevel(40)

    # setup log file
    if logFile is not None:
        logFile = open(logFile, 'w+')
        print(f'\nStart loading fMRI data at '+time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))+'\n', file=logFile, flush=True)

    if os.path.isfile(file_scan_list) and file_scan_list.endswith('.txt'):
        scan_list = [line.replace('\n', '') for line in open(file_scan_list, "r")]
    else:
        scan_list = [file_scan_list]

    Data = None
    for i in range(len(scan_list)):
        if len(scan_list[i]) == 0:
            break

        if logFile is not None:
            print(f' loading scan ' + scan_list[i], file=logFile)

        file_list = str.split(scan_list[i], ';')
        for file in file_list:
            if not os.path.isfile(file):
                raise ValueError('The file does not exist: ' + file)

        # Loading a single fMRI scan
        # 2D or 4D matrix with dimension definition [dim_space dim_time] or [X, Y, Z, T]

        # 'HCP Surface (*.cifti, *.mat)'
        if dataFormat == 'HCP Surface (*.cifti, *.mat)':
            if scan_list[i].endswith('.dtseries.nii'):
                cifti = nib.load(scan_list[i])  # [dim_time dim_space]
                cifti_data = cifti.get_fdata(dtype=np.float32)
                # Extract desired parts of the data
                if dataType == 'Surface':
                    scan_data = cifti_data[:, range(59412)]
                else:
                    raise ValueError('Unsupported data type ' + dataType + ' for data format HCP Surface')

            elif scan_list[i].endswith('.mat'):
                scan_data = load_matlab_single_array(scan_list[i])  # [dim_space dim_time]
                if scan_data.shape[0] < 59412:
                    raise ValueError('The MATLAB file contains a 2D matrix with the spatial dimension smaller than 59412 in file ' + scan_list[i])
                scan_data = scan_data[range(59412), :].T

            else:
                raise ValueError('Unsupported data format ' + scan_list[i])

        elif dataFormat == 'MGH Surface (*.mgh)':
            # need to split each line to two directories for left and right hemispheres
            if len(str.split(scan_list[i], ';')) != 2:
                raise ValueError("For MGH surface data format, directories of two hemisphere data need to be combined into one line with ';' as separator")

            # get files for two hemispheres
            file_L = str.split(scan_list[i], ';')[0]
            file_R = str.split(scan_list[i], ';')[1]

            # check extension
            if not file_L.endswith('.mgh') or not file_R.endswith('.mgh'):
                raise ValueError('For MGH surface format, the file extension should be .mgh')

            scan_data = np.array(nib.load(file_L).get_fdata(dtype=np.float32))  # [dim_space, _, _, dim_time]
            scan_data = np.squeeze(scan_data)
            scan_data = np.append(scan_data, np.squeeze(np.array(nib.load(file_R).get_fdata(dtype=np.float32))), axis=0)
            scan_data = scan_data.T

        elif dataFormat == 'MGZ Surface (*.mgz)':
            mgz = nib.load(scan_list[i])
            scan_data = mgz.get_fdata(dtype=np.float32)  # [dim_space dim_time]
            scan_data = np.squeeze(np.array(scan_data)).T

        elif dataFormat == 'Volume (*.nii, *.nii.gz, *.mat)':
            if len(scan_list) > 1 and Reshape is False and Concatenation:
                raise ValueError('4D fMRI data must be reshaped to 2D first before concatenation')

            if scan_list[i].endswith('.nii') or scan_list[i].endswith('.nii.gz'):
                nii = nib.load(scan_list[i])
                scan_data = nii.get_fdata(dtype=np.float32)
            elif scan_list[i].endswith('.mat'):
                scan_data = load_matlab_single_array(scan_list[i])  # [X Y Z dim_time]
            else:
                raise ValueError('Unsupported data format ' + scan_list[i])

            if Reshape:
                if Brain_Mask is None:
                    raise ValueError('Brain_Mask must be provided when Reshape is enabled for 4D fMRI data')
                scan_data = reshape_fmri_data(scan_data, dataType, Brain_Mask)

        else:
            raise ValueError('Unsupported data format ' + dataFormat)

        # scan_data should be in [dim_time dim_space]
        # Convert to NumPy array
        if not isinstance(scan_data, np.ndarray):
            scan_data = np.array(scan_data)
        if logFile is not None:
            print(f' loaded data size is ' + str(scan_data.shape), file=logFile)

        # Combine scans along the time dimension
        # The Data will be permuted to [dim_time dim_space] for both 2D and 4D matrices
        if i == 0:
            if dataType == 'Surface':
                if Normalization is not None and Normalization is not False:
                    if Normalization == 'vp-vmax':
                        scan_data = normalize_data(scan_data, 'vp', 'vmax')
                    else:
                        raise ValueError('Unsupported data normalization: ' + Normalization)
            elif dataType == 'Volume':
                if Normalization is not None and Normalization is not False:
                    if Normalization == 'vp-vmax':
                        scan_data = normalize_data(scan_data, 'vp', 'vmax')
                    else:
                        raise ValueError('Unsupported data normalization: ' + Normalization)
            Data = scan_data

        else:
            if dataType == 'Surface':
                if Data is None or len(Data.shape) != 2 or scan_data.shape[1] != Data.shape[1]:
                    raise ValueError('Scans have different spatial dimensions when loading scan: ' + scan_list[i])
                if Normalization is not None and Normalization is not False:
                    if Normalization == 'vp-vmax':
                        scan_data = normalize_data(scan_data, 'vp', 'vmax')
                    else:
                        raise ValueError('Unsupported data normalization: ' + Normalization)

            elif dataType == 'Volume':
                if Data is None or len(Data.shape) != 2 or scan_data.shape[1] != Data.shape[1]:
                    raise ValueError('Scans have different spatial dimensions when loading scan: ' + scan_list[i])
                if Reshape is False:
                    raise ValueError('4D volume-based fMRI scans need to be reshaped before concatenation')
                if Normalization is not None and Normalization is not False:
                    if Normalization == 'vp-vmax':
                        scan_data = normalize_data(scan_data, 'vp', 'vmax')
                    else:
                        raise ValueError('Unsupported data normalization: ' + Normalization)
            else:
                raise ValueError('Unknown dataType: ' + dataType)

            if Concatenation:
                Data = np.append(Data, scan_data, axis=0)
            else:
                raise ValueError('Only supports to concatenate data for output')

    if logFile is not None:
        print('\nConcatenated data is a 2D matrix with size ' + str(Data.shape), file=logFile)
    return Data


def compute_brain_surface(file_surfL: str, file_surfR: str, file_maskL: str, file_maskR: str, file_surfL_inflated=None, file_surfR_inflated=None,
                          maskValue=0, dataType='Surface', dataFormat='HCP Surface (*.cifti, *.mat)', logFile=None):
    """
    Prepare a brain surface variable to store surface shape (vertices and faces), and brain masks for useful vertices

    :param file_surfL: file that stores the surface shape information of the left hemisphere, including vertices and faces
    :param file_surfR: file that stores the surface shape information of the right hemisphere, including vertices and faces
    :param file_maskL: file that stores the mask information of the left hemisphere, a 1D 0-1 vector
    :param file_maskR: file that stores the mask information of the right hemisphere, a 1D 0-1 vector
    :param file_surfL_inflated: file that stores the inflated surface shape information of the left hemisphere, including vertices and faces
    :param file_surfR_inflated: file that stores the inflated surface shape information of the right hemisphere, including vertices and faces
    :param maskValue: 0 or 1, 0 means 0s in mask files are useful vertices, otherwise vice versa. maskValue=0 for medial wall in HCP data, and maskValue=1 for brain masks
    :param dataType: 'Surface'
    :param dataFormat: 'HCP Surface (*.cifti, *.mat)' or 'FreeSurfer Surface (*.)'
    :param logFile:
    :return: Brain_Surface: a structure with keys Data_Type, Data_Format, Shape (including L and R), Shape_Inflated (if used), Mask (including L and R)

    Yuncong Ma, 9/25/2023
    """

    if dataType == 'Surface' and dataFormat == 'HCP Surface (*.cifti, *.mat)':

        shapeL = nib.load(file_surfL)
        shapeR = nib.load(file_surfR)
        maskL = nib.load(file_maskL)
        maskR = nib.load(file_maskR)

        # Initialize Brain_Surface
        if file_surfL_inflated is not None and file_surfR_inflated is not None:
            shapeL_inflated = nib.load(file_surfL_inflated)
            shapeR_inflated = nib.load(file_surfR_inflated)
            Brain_Surface = {'Data_Type': dataType, 'Data_Format': dataFormat,
                             'Shape': {'L': {'vertices': [], 'faces': []}, 'R': {'vertices': [], 'faces': []}},
                             'Shape_Inflated': {'L': {'vertices': [], 'faces': []}, 'R': {'vertices': [], 'faces': []}},
                             'Brain_Mask': {'L': [], 'R': []}}
        else:
            Brain_Surface = {'Data_Type': 'Surface', 'Data_Format': dataFormat,
                             'Shape': {'L': {'vertices': [], 'faces': []}, 'R': {'vertices': [], 'faces': []}},
                             'Brain_Mask': {'L': [], 'R': []}}
        # Surface shape
        # Index starts from 1
        Brain_Surface['Shape']['L']['vertices'], Brain_Surface['Shape']['L']['faces'] = shapeL.darrays[0].data, shapeL.darrays[1].data + int(1)
        Brain_Surface['Shape']['R']['vertices'], Brain_Surface['Shape']['R']['faces'] = shapeR.darrays[0].data, shapeR.darrays[1].data + int(1)
        # Surface shape inflated
        if file_surfL_inflated is not None and file_surfR_inflated is not None:
            Brain_Surface['Shape_Inflated']['L']['vertices'], Brain_Surface['Shape_Inflated']['L']['faces'] = shapeL_inflated.darrays[0].data, shapeL_inflated.darrays[1].data + int(1)
            Brain_Surface['Shape_Inflated']['R']['vertices'], Brain_Surface['Shape_Inflated']['R']['faces'] = shapeR_inflated.darrays[0].data, shapeR_inflated.darrays[1].data + int(1)
        # Index for brain mask
        Brain_Surface['Brain_Mask']['L'] = maskL.darrays[0].data
        Brain_Surface['Brain_Mask']['R'] = maskR.darrays[0].data

        # Change 0 to 1 if the mask is to label unused vertices
        if maskValue == 0:
            Brain_Surface['Brain_Mask']['L'] = (Brain_Surface['Brain_Mask']['L'] == 0)
            Brain_Surface['Brain_Mask']['R'] = (Brain_Surface['Brain_Mask']['R'] == 0)

    else:
        raise ValueError('Unknown combination of Data_Type and Data_Surface: ' + dataType + ' : ' + dataFormat)

    if logFile is not None:
        print('\nBrain_Surface is created', file=logFile, flush=True)

    return Brain_Surface


def compute_brain_template(dataType: str, dataFormat: str, file_surfL=None, file_surfR=None, file_maskL=None, file_maskR=None,
                           file_mask_vol=None, file_overlayImage=None, maskValue=0, file_surfL_inflated=None, file_surfR_inflated=None,
                           logFile=None):
    """
    Prepare a brain surface variable to store surface shape (vertices and faces), and brain masks for useful vertices

    :param dataType: 'Surface', 'Volume'
    :param dataFormat: 'HCP Surface (*.cifti, *.mat)', 'FreeSurfer Surface (*.)', or 'Volume (*.nii, *.nii.gz, *.mat)'
    :param file_surfL: file that stores the surface shape information of the left hemisphere, including vertices and faces
    :param file_surfR: file that stores the surface shape information of the right hemisphere, including vertices and faces
    :param file_maskL: file that stores the mask information of the left hemisphere, a 1D 0-1 vector
    :param file_maskR: file that stores the mask information of the right hemisphere, a 1D 0-1 vector
    :param file_surfL_inflated: file that stores the inflated surface shape information of the left hemisphere, including vertices and faces
    :param file_surfR_inflated: file that stores the inflated surface shape information of the right hemisphere, including vertices and faces
    :param file_mask_vol: file of a mask file for volume-based data type
    :param file_overlayImage: file of a background image for visualizing volume-based results
    :param maskValue: 0 or 1, 0 means 0s in mask files are useful vertices, otherwise vice versa. maskValue=0 for medial wall in HCP data, and maskValue=1 for brain masks

    :param logFile:
    :return: Brain_Template: a structure with keys Data_Type, Data_Format, Shape (including L and R), Shape_Inflated (if used), Mask (including L and R) for surface type
                            a structure with keys Data_Type, Data_Format, Mask, Overlay_Image

    Yuncong Ma, 9/24/2023
    """

    if dataType == 'Volume':
        if file_mask_vol is None or file_overlayImage is None:
            raise ValueError('When data type is volume, both file_mask_vol and file_overlayImage are required')
        Brain_Mask = load_fmri_scan(file_mask_vol, dataType=dataType, dataFormat=dataFormat, Reshape=False, Normalization=None)
        Overlay_Image = load_fmri_scan(file_overlayImage, dataType=dataType, dataFormat=dataFormat, Reshape=False, Normalization=None)
        Brain_Mask = (Brain_Mask == maskValue)
        Brain_Template = {'Data_Type': dataType,
                          'Data_Format': dataFormat,
                          'Brain_Mask': Brain_Mask,
                          'Overlay_Image': Overlay_Image}

    elif dataType == 'Surface':
        if file_surfL is None or file_surfR is None or file_maskL is None or file_maskR is None:
            raise ValueError('When data type is surface, file_surfL, file_surfR, file_maskL and file_maskR are required')
        Brain_Template = \
            compute_brain_surface(file_surfL, file_surfR, file_maskL, file_maskR,
                                  file_surfL_inflated=file_surfL_inflated, file_surfR_inflated=file_surfR_inflated,
                                  maskValue=maskValue,
                                  dataType=dataType, dataFormat=dataFormat,
                                  logFile=logFile)

    else:
        raise ValueError('Unknown data type: ' + dataType)

    return Brain_Template


def setup_brain_template(dir_pnet_dataInput: str, Brain_Template):
    """
    Save the Brain_Template.mat

    :param dir_pnet_dataInput: the directory of the Data Input folder
    :param Brain_Template: a structure created by function compute_brain_template

    Yuncong Ma, 9/24/2023
    """

    # Use both matlab and json files for convenience
    scipy.io.savemat(os.path.join(dir_pnet_dataInput, 'Brain_Template.mat'), {'Brain_Template': Brain_Template})
    if Brain_Template['Data_Type'] == 'Volume':
        Brain_Template['Brain_Mask'] = Brain_Template['Brain_Mask'].tolist()
        Brain_Template['Overlay_Image'] = Brain_Template['Overlay_Image'].tolist()
        write_json_setting(Brain_Template, os.path.join(dir_pnet_dataInput, 'Brain_Template.json'))

    elif Brain_Template['Data_Type'] == 'Surface':
        Brain_Template['Brain_Mask']['L'] = Brain_Template['Brain_Mask']['L'].tolist()
        Brain_Template['Brain_Mask']['R'] = Brain_Template['Brain_Mask']['R'].tolist()
        Brain_Template['Shape']['L']['vertices'] = Brain_Template['Shape']['L']['vertices'].tolist()
        Brain_Template['Shape']['R']['vertices'] = Brain_Template['Shape']['R']['vertices'].tolist()
        Brain_Template['Shape']['L']['faces'] = Brain_Template['Shape']['L']['faces'].tolist()
        Brain_Template['Shape']['R']['faces'] = Brain_Template['Shape']['R']['faces'].tolist()
        if 'Shape_Inflated' in Brain_Template.keys():
            Brain_Template['Shape_Inflated']['L']['vertices'] = Brain_Template['Shape_Inflated']['L']['vertices'].tolist()
            Brain_Template['Shape_Inflated']['R']['vertices'] = Brain_Template['Shape_Inflated']['R']['vertices'].tolist()
            Brain_Template['Shape_Inflated']['L']['faces'] = Brain_Template['Shape_Inflated']['L']['faces'].tolist()
            Brain_Template['Shape_Inflated']['R']['faces'] = Brain_Template['Shape_Inflated']['R']['faces'].tolist()
        write_json_setting(Brain_Template, os.path.join(dir_pnet_dataInput, 'Brain_Template.json'))

    else:
        raise ValueError('Unsupported data type: ' + Brain_Template['Data_Type'])


def load_brain_template(dir_Brain_Template: str, logFile=None):
    """
    Load a brain template file

    :param dir_Brain_Template: directory of the brain_template file, in json format. Python cannot read the MATLAB version
    :param logFile: directory of a log file
    :return: Brain_Template: nested dictionary storing information and matrices of brain template. Matrices are converted to np.ndarray

    Yuncong Ma, 9/25/2023
    """

    Brain_Template = load_json_setting(dir_Brain_Template)

    # Check Brain_Template
    if 'Data_Type' not in Brain_Template.keys():
        raise ValueError('Cannot find Data_type in the Brain_Template file')

    # Convert list to np.ndarray
    if Brain_Template['Data_Type'] == 'Volume':
        if 'Brain_Mask' not in Brain_Template.keys() or 'Overlay_Image' not in Brain_Template.keys():
            raise ValueError('Cannot find Brain_Mask and Overlay_Image in the Brain_Template file')
        # convert to np.ndarray
        Brain_Template['Brain_Mask'] = np.array(Brain_Template['Brain_Mask'])
        Brain_Template['Overlay_Image'] = np.array(Brain_Template['Overlay_Image'])

    elif Brain_Template['Data_Type'] == 'Surface':
        # check sub-keys
        if 'Brain_Mask' not in Brain_Template.keys() or 'Shape' not in Brain_Template.keys():
            raise ValueError('Cannot find Brain_Mask and Shape in the Brain_Template file')
        if 'L' not in Brain_Template['Shape'].keys() or 'R' not in Brain_Template['Shape'].keys():
            raise ValueError("Cannot find L or R in Brain_Template['Shape']")
        if 'vertices' not in Brain_Template['Shape']['L'].keys() or 'vertices' not in Brain_Template['Shape']['L'].keys():
            raise ValueError("Cannot find vertices in Brain_Template['Shape']['L']")
        if 'faces' not in Brain_Template['Shape']['L'].keys() or 'faces' not in Brain_Template['Shape']['L'].keys():
            raise ValueError("Cannot find faces in Brain_Template['Shape']['L']")
        if 'vertices' not in Brain_Template['Shape']['R'].keys() or 'vertices' not in Brain_Template['Shape']['R'].keys():
            raise ValueError("Cannot find vertices in Brain_Template['Shape']['R']")
        if 'faces' not in Brain_Template['Shape']['R'].keys() or 'faces' not in Brain_Template['Shape']['R'].keys():
            raise ValueError("Cannot find faces in Brain_Template['Shape']['R']")
        # convert to np.ndarray
        Brain_Template['Brain_Mask']['L'] = np.array(Brain_Template['Brain_Mask']['L'])
        Brain_Template['Brain_Mask']['R'] = np.array(Brain_Template['Brain_Mask']['R'])
        Brain_Template['Shape']['L']['vertices'] = np.array(Brain_Template['Shape']['L']['vertices'])
        Brain_Template['Shape']['R']['vertices'] = np.array(Brain_Template['Shape']['R']['vertices'])
        Brain_Template['Shape']['L']['faces'] = np.array(Brain_Template['Shape']['L']['faces'])
        Brain_Template['Shape']['R']['faces'] = np.array(Brain_Template['Shape']['R']['faces'])
        if 'Shape_Inflated' in Brain_Template.keys():
            # check sub-keys
            if 'vertices' not in Brain_Template['Shape_Inflated']['L'].keys() or 'vertices' not in Brain_Template['Shape_Inflated']['L'].keys():
                raise ValueError("Cannot find vertices in Brain_Template['Shape_Inflated']['L']")
            if 'faces' not in Brain_Template['Shape_Inflated']['L'].keys() or 'faces' not in Brain_Template['Shape_Inflated']['L'].keys():
                raise ValueError("Cannot find faces in Brain_Template['Shape_Inflated']['L']")
            if 'vertices' not in Brain_Template['Shape_Inflated']['R'].keys() or 'vertices' not in Brain_Template['Shape_Inflated']['R'].keys():
                raise ValueError("Cannot find vertices in Brain_Template['Shape_Inflated']['R']")
            if 'faces' not in Brain_Template['Shape_Inflated']['R'].keys() or 'faces' not in Brain_Template['Shape_Inflated']['R'].keys():
                raise ValueError("Cannot find faces in Brain_Template['Shape_Inflated']['R']")
            # convert to np.ndarray
            Brain_Template['Shape_Inflated']['L']['vertices'] = np.array(Brain_Template['Shape_Inflated']['L']['vertices'])
            Brain_Template['Shape_Inflated']['R']['vertices'] = np.array(Brain_Template['Shape_Inflated']['R']['vertices'])
            Brain_Template['Shape_Inflated']['L']['faces'] = np.array(Brain_Template['Shape_Inflated']['L']['faces'])
            Brain_Template['Shape_Inflated']['R']['faces'] = np.array(Brain_Template['Shape_Inflated']['R']['faces'])

    else:
        raise ValueError('Unsupported data type: ' + Brain_Template['Data_Type'])

    return Brain_Template


def reshape_fmri_data(scan_data: np.ndarray, dataType: str, Brain_Mask: np.ndarray, logFile=None):
    """
    reshape_fmri_data(scan_data: np.ndarray, dataType: str, Brain_Mask: np.ndarray, logFile=None)
    Reshape 4D volume fMRI data [X Y Z dim_time] into 2D matrix [dim_time dim_space]
    Reshape 2D fMRI data back to 4D volume type

    :param scan_data: 4D or 2D matrix [X Y Z dim_time] [dim_time dim_space]
    :param dataType: 'Surface' or 'Volume'
    :param Brain_Mask: 3D matrix
    :param logFile:
    :return: reshaped_data: 2D matrix if input is 4D, vice versa

    Yuncong Ma, 9/13/2023
    """

    if dataType == 'Volume':
        if len(scan_data.shape) == 4:  # 4D fMRI data, reshape to 2D [dim_time, dim_space]
            if scan_data.shape[0:3] != Brain_Mask.shape:
                raise ValueError('The shapes of Brain_Mask and scan_data are not the same when scan_data is a 4D matrix')
            reshaped_data = np.reshape(scan_data, (np.prod(scan_data.shape[0:3]), scan_data.shape[3]), order='F').T   # Match colum based index used in MATLAB
            reshaped_data = reshaped_data[:, Brain_Mask.flatten('F') > 0]   # Match colum based index used in MATLAB

        elif len(scan_data.shape) == 2:  # 2D fMRI data, reshape back to 4D [X Y Z T]
            if scan_data.shape[1] != np.sum(Brain_Mask > 0):
                raise ValueError('The nodes in Brain_Mask and scan_data are not the same when scan_data is a 2D matrix')
            ps = (Brain_Mask > 0).flatten('F')   # Match colum based index used in MATLAB
            reshaped_data = np.zeros((scan_data.shape[0], len(ps)), like=scan_data)
            reshaped_data[:, ps] = scan_data
            reshaped_data = np.reshape(reshaped_data.T, (Brain_Mask.shape[0], Brain_Mask.shape[1], Brain_Mask.shape[2], scan_data.shape[0]), order='F')   # Match colum based index used in MATLAB

        else:
            raise ValueError('The scan_data needs to be a 2D or 4D matrix')

    else:
        reshaped_data = scan_data

    return reshaped_data


def reshape_FN(FN: np.ndarray, dataType: str, Brain_Mask: np.ndarray, logFile=None):
    """
    reshape_fmri_data(scan_data: np.ndarray, dataType: str, Brain_Mask: np.ndarray, logFile=None)
    If dataType is 'Volume'
    Reshape 4D FNs [X Y Z dim_time] into 2D matrix [dim_time dim_space], extracting voxels in Brain_Mask
    Reshape 2D FNs back to 4D for storage and visualization

    :param FN: 4D or 2D matrix [X Y Z K] [dim_space K]
    :param dataType: 'Surface' or 'Volume'
    :param Brain_Mask: 3D matrix
    :param logFile:
    :return: reshaped_FN: 2D matrix if input is 4D, vice versa

    Yuncong Ma, 9/13/2023
    """

    if dataType == 'Volume':
        if len(FN.shape) == 4:  # 4D FN [X Y Z K], reshape to 2D [dim_space, K]
            if FN.shape[0:3] != Brain_Mask.shape:
                raise ValueError('The shapes of Brain_Mask and FN are not the same when scan_data is a 4D matrix')
            reshaped_FN = np.reshape(FN, (np.prod(FN.shape[0:3]), FN.shape[3]), order='F')   # Match colum based index used in MATLAB
            reshaped_FN = reshaped_FN[Brain_Mask.flatten('F') > 0, :]   # Match colum based index used in MATLAB

        elif len(FN.shape) == 2:  # 2D FN [dim_space, K], reshape back to 4D [X Y Z K]
            if FN.shape[0] != np.sum(Brain_Mask > 0):
                raise ValueError('The nodes in Brain_Mask and scan_data are not the same when scan_data is a 2D matrix')
            ps = (Brain_Mask > 0).flatten('F')   # Match colum based index used in MATLAB
            reshaped_FN = np.zeros((len(ps), FN.shape[1]), like=FN)
            reshaped_FN[ps, :] = FN
            reshaped_FN = np.reshape(reshaped_FN, (Brain_Mask.shape[0], Brain_Mask.shape[1], Brain_Mask.shape[2], FN.shape[1]), order='F')   # Match colum based index used in MATLAB

        else:
            raise ValueError('The scan_data needs to be a 2D or 4D matrix')

    else:
        reshaped_FN = FN

    return reshaped_FN


def setup_result_folder(dir_pnet_result: str):
    """
    setup_result_folder(dir_pnet_result: str)
    Setup sub-folders in the pNet results folder to store setting and results
    Including, Data_Input, FN_Computation, Group_FN, Personalized_FN, Quality_Control, Statistics

    :param dir_pnet_result:
    :return: dir_pnet_dataInput, dir_pnet_FNC, dir_pnet_gFN, dir_pnet_pFN, dir_pnet_QC, dir_pnet_STAT

    Yuncong Ma, 9/13/2023
    """

    if not os.path.exists(dir_pnet_result):
        os.makedirs(dir_pnet_result)

    # Sub-folders
    dir_pnet_dataInput = os.path.join(dir_pnet_result, 'Data_Input')
    if not os.path.exists(dir_pnet_dataInput):
        os.makedirs(dir_pnet_dataInput)
    dir_pnet_FNC = os.path.join(dir_pnet_result, 'FN_Computation')
    if not os.path.exists(dir_pnet_FNC):
        os.makedirs(dir_pnet_FNC)
    dir_pnet_gFN = os.path.join(dir_pnet_result, 'Group_FN')
    if not os.path.exists(dir_pnet_gFN):
        os.makedirs(dir_pnet_gFN)
    dir_pnet_pFN = os.path.join(dir_pnet_result, 'Personalized_FN')
    if not os.path.exists(dir_pnet_pFN):
        os.makedirs(dir_pnet_pFN)
    dir_pnet_QC = os.path.join(dir_pnet_result, 'Quality_Control')
    if not os.path.exists(dir_pnet_QC):
        os.makedirs(dir_pnet_QC)
    dir_pnet_STAT = os.path.join(dir_pnet_result, 'Statistics')
    if not os.path.exists(dir_pnet_STAT):
        os.makedirs(dir_pnet_STAT)

    return dir_pnet_dataInput, dir_pnet_FNC, dir_pnet_gFN, dir_pnet_pFN, dir_pnet_QC, dir_pnet_STAT


def check_data_type_format(dataType: str,  dataFormat: str, logFile=None):
    """
    check_data_type_format(dataType: str,  dataFormat: str, logFile=None)
    Check setting for dataType and dataFormat

    :param dataType: 'Volume', 'Surface'
    :param dataFormat: 'HCP Surface (*.cifti, *.mat)', 'MGH Surface (*.mgh)', 'MGZ Surface (*.mgz)', 'Volume (*.nii, *.nii.gz, *.mat)'
    :param logFile:

    Yuncong Ma, 9/25/2023
    """

    if dataType not in ('Volume', 'Surface'):
        if logFile is None:
            raise ValueError("Data type should be 'Volume' or 'Surface'")
        else:
            logFile = open(logFile, 'a')
            print("Data type should be 'Volume' or 'Surface'", file=logFile, flush=True)

    if dataFormat not in ('HCP Surface (*.cifti, *.mat)', 'MGH Surface (*.mgh)', 'MGZ Surface (*.mgz)', 'Volume (*.nii, *.nii.gz, *.mat)'):
        if logFile is None:
            raise ValueError("Data format should be 'HCP Surface (*.cifti, *.mat)', 'MGH Surface (*.mgh)', 'MGZ Surface (*.mgz)', or 'Volume (*.nii, *.nii.gz, *.mat)'")
        else:
            logFile = open(logFile, 'a')
            print("Data format should be 'HCP Surface (*.cifti, *.mat)', 'MGH Surface (*.mgh)', 'MGZ Surface (*.mgz)', or 'Volume (*.nii, *.nii.gz, *.mat)'", file=logFile, flush=True)


def setup_scan_info(dir_pnet_dataInput: str, dataType: str, dataFormat: str, file_scan: str, file_subject_ID=None, file_subject_folder=None, file_group_ID=None, scan_info='Automatic', Combine_Scan=False, logFile=None):
    """
    setup_scan_info(dir_pnet_dataInput: str, file_scan: str, file_subject_ID=None, file_subject_folder=None, file_group=None, scan_info='Automatic',Combine_Scan=False)
    Set up a few txt files for labeling scans
    file_scan contains the directories of all fMRI scans
    file_subject_ID contains the subject ID information for each corresponding fMRI scans in file_scan
    file_subject_folder contains the sub-folder names for each corresponding fMRI scans in file_scan, in order to store results for each fMRI scan or combined scans
    file_Group_ID contains the group ID for each corresponding fMRI scans in file_scan, in order to do desired sampling based on the group information

    
    :param dir_pnet_dataInput: directory of the Data_Input folder or a directory to save Scan_List.txt, Subject_ID.txt, Subject_Folder.txt and Group_ID.txt
    :param dataType: 'Surface', 'Volume'
    :param dataFormat: 'HCP Surface (*.cifti, *.mat)', 'MGH Surface (*.mgh)', 'MGZ Surface (*.mgz)', or 'Volume (*.nii, *.nii.gz, *.mat)'
    :param file_scan: a txt file that stores directories of all fMRI scans
    :param file_subject_ID: a txt file that store subject ID information corresponding to fMRI scan in file_scan
    :param file_subject_folder: a txt file that store subject folder names corresponding to fMRI scans in file_scan
    :param file_group_ID: a txt file that store group information corresponding to fMRI scan in file_scan
    :param scan_info: 'Automatic' or 'Manual', 'Manual' requires manual input of file_subject_ID, file_subject_folder and file_group
    :param Combine_Scan: False or True, whether to combine multiple scans for the same subject
    :param logFile: directory of a txt formatted log file

    Yuncong Ma, 9/25/2023
    """

    # Check setting
    check_data_type_format(dataType, dataFormat, logFile=logFile)

    # output setting file
    setting = {'Data_Type': dataType, 'Data_Format': dataFormat,
               'file_scan': file_scan, 'file_subject_ID': file_subject_ID, 'file_subject_folder': file_subject_folder,
               'file_group_ID': file_group_ID, 'scan_info': scan_info, 'Combine_Scan': Combine_Scan}
    write_json_setting(setting, os.path.join(dir_pnet_dataInput, 'Setting.json'))

    if scan_info == 'Manual':
        # file_subject_ID is required for Manual setting
        if file_subject_ID is None:
            raise ValueError('When scan_info is set to "Manual", file_subject_ID is required')
        if file_subject_folder is None:
            raise ValueError('When scan_info is set to "Manual", file_subject_folder is required')

        # copy files to dir_pnet_dataInput
        if os.path.join(dir_pnet_dataInput, 'Scan_List.txt') != file_scan:
            dest_file_scan = os.path.join(dir_pnet_dataInput, 'Scan_List.txt')
            dest_file_scan = open(dest_file_scan, 'w')
            [print(line.replace('\n', ''), file=dest_file_scan) for line in open(file_scan, 'r')]
            dest_file_scan.close()
        if os.path.join(dir_pnet_dataInput, 'Subject_ID.txt') != file_subject_ID:
            dest_file_subject_ID = os.path.join(dir_pnet_dataInput, 'Subject_ID.txt')
            dest_file_subject_ID = open(dest_file_subject_ID, 'w')
            [print(line.replace('\n', ''), file=dest_file_subject_ID) for line in open(file_subject_ID, 'r')]
            dest_file_subject_ID.close()
        if os.path.join(dir_pnet_dataInput, 'Subject_Folder.txt') != file_subject_folder:
            dest_file_subject_folder = os.path.join(dir_pnet_dataInput, 'Subject_Folder.txt')
            dest_file_subject_folder = open(dest_file_subject_folder, 'w')
            [print(line.replace('\n', ''), file=dest_file_subject_folder) for line in open(file_subject_folder, 'r')]
            dest_file_subject_folder.close()
        if file_group_ID is not None and os.path.join(dir_pnet_dataInput, 'Group_ID.txt') != file_group_ID:
            dest_file_group_ID = os.path.join(dir_pnet_dataInput, 'Group_ID.txt')
            dest_file_group_ID = open(dest_file_group_ID, 'w')
            [print(line.replace('\n', ''), file=dest_file_group_ID) for line in open(file_group_ID, 'r')]
            dest_file_group_ID.close()

    elif scan_info == 'Automatic':
        # read the scan list file
        list_scan = [line.replace('\n', '') for line in open(file_scan, 'r')]
        N_Scan = len(list_scan)
        # Automatically extract the common directory in list_scan as the root_directory to generate subject_ID and subject_folder
        # Subject ID is extracted as the first level sub-folder names
        # Subject folder will be either subject_ID/ numbers starting from 1 to N, or subject_ID when Combine_Scan is enabled
        if file_subject_ID is None and file_subject_folder is None:
            common_prefix = os.path.commonprefix(list_scan)
            root_directory = os.path.dirname(common_prefix)
            list_subject_ID = [os.path.normpath(list_scan[i][len(root_directory)+1 : -1]).split(os.path.sep)[0] for i in range(N_Scan)]
            list_subject_ID_unique = np.unique(np.array(list_subject_ID))
            N_Subject = list_subject_ID_unique.shape[0]
            list_subject_folder = list_subject_ID.copy()
            if not Combine_Scan:
                for i in range(N_Subject):
                    ps = [j for j, x in enumerate(list_subject_ID) if x == list_subject_ID_unique[i]]
                    for j in range(len(ps)):
                        list_subject_folder[ps[j]] = os.path.join(list_subject_folder[ps[j]], str(j+1))

            # Output
            # copy files to dir_pnet_dataInput
            if os.path.join(dir_pnet_dataInput, 'Scan_List.txt') != file_scan:
                dest_file_scan = os.path.join(dir_pnet_dataInput, 'Scan_List.txt')
                dest_file_scan = open(dest_file_scan, 'w')
                [print(line.replace('\n', ''), file=dest_file_scan) for line in open(file_scan, 'r')]
                dest_file_scan.close()
            file_subject_ID = os.path.join(dir_pnet_dataInput, 'Subject_ID.txt')
            file_subject_ID = open(file_subject_ID, 'w')
            for i in range(len(list_subject_ID)):
                print(list_subject_ID[i], file=file_subject_ID)
            file_subject_ID.close()
            file_subject_folder = os.path.join(dir_pnet_dataInput, 'Subject_Folder.txt')
            file_subject_folder = open(file_subject_folder, 'w')
            for i in range(len(list_subject_folder)):
                print(list_subject_folder[i], file=file_subject_folder)
            file_subject_folder.close()










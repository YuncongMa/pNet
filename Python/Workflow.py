# Yuncong Ma, 12/1/2023
# pNet
# Provide examples of running the whole workflow of pNet

#########################################
# Packages

import argparse

# Module
# This script builds the five modules of pNet
# Functions for modules of pNet
from Brain_Template import Brain_Template
from Data_Input import *
from FN_Computation import *
from FN_Computation_torch import *
from Computation_Environment import *
from Quality_Control import *
from Quality_Control_torch import *
from Visualization import *
from Server import *


def workflow(dir_pnet_result: str,
             file_scan: str,
             dataType='Surface', dataFormat='HCP Surface (*.cifti, *.mat)',
             file_subject_ID=None, file_subject_folder=None, file_group_ID=None,
             file_Brain_Template=None,
             templateFormat='HCP',
             file_surfL=None, file_surfR=None, file_maskL=None, file_maskR=None,
             file_mask_vol=None, file_overlayImage=None,
             maskValue=0,
             file_surfL_inflated=None, file_surfR_inflated=None,
             K=17, Combine_Scan=False,
             file_gFN=None,
             samplingMethod='Subject', sampleSize=10, nBS=50,
             maxIter=1000, minIter=30, meanFitRatio=0.1, error=1e-8, normW=1,
             Alpha=2, Beta=30, alphaS=0, alphaL=0, vxI=0, ard=0, eta=0, nRepeat=5,
             Parallel=False, Computation_Mode='CPU_Torch', N_Thread=1,
             dataPrecision='double',
             outputFormat='Both',
             synchronized_view=True or bool,
             synchronized_colorbar=True or bool):
    """
    Run the workflow of pNet, including Data Input, FN Computation, Quality Control and Visualization
    This function is for running pNet in a single job

    :param dir_pnet_result: directory of the pNet result folder
    :param dataType: 'Surface', 'Volume', 'Surface-Volume'
    :param dataFormat: 'HCP Surface (*.cifti, *.mat)', 'MGH Surface (*.mgh)', 'MGZ Surface (*.mgz)', 'Volume (*.nii, *.nii.gz, *.mat)', 'HCP Surface-Volume (*.cifti)', 'HCP Volume (*.cifti)'

    :param file_scan: a txt file that stores directories of all fMRI scans
    :param file_subject_ID: a txt file that store subject ID information corresponding to fMRI scan in file_scan
    :param file_subject_folder: a txt file that store subject folder names corresponding to fMRI scans in file_scan
    :param file_group_ID: a txt file that store group information corresponding to fMRI scan in file_scan

    :param file_Brain_Template: file directory of a brain template file in json format
    :param templateFormat: 'HCP', 'FreeSurfer', '3D Matrix'
    :param file_surfL: file that stores the surface shape information of the left hemisphere, including vertices and faces
    :param file_surfR: file that stores the surface shape information of the right hemisphere, including vertices and faces
    :param file_maskL: file that stores the mask information of the left hemisphere, a 1D 0-1 vector
    :param file_maskR: file that stores the mask information of the right hemisphere, a 1D 0-1 vector
    :param file_surfL_inflated: file that stores the inflated surface shape information of the left hemisphere, including vertices and faces
    :param file_surfR_inflated: file that stores the inflated surface shape information of the right hemisphere, including vertices and faces
    :param file_mask_vol: file of a mask file for volume-based data type
    :param file_overlayImage: file of a background image for visualizing volume-based results
    :param maskValue: 0 or 1, 0 means 0s in mask files are useful vertices, otherwise vice versa. maskValue=0 for medial wall in HCP data, and maskValue=1 for brain masks

    :param K: number of FNs
    :param Combine_Scan: False or True, whether to combine multiple scans for the same subject

    :param file_gFN: None or a directory of a precomputed gFN in .mat format
    :param samplingMethod: 'Subject' or 'Group_Subject'. Uniform sampling based subject ID, or group and then subject ID
    :param sampleSize: number of subjects selected for each bootstrapping run
    :param nBS: number of runs for bootstrap

    :param maxIter: maximum iteration number for multiplicative update
    :param minIter: minimum iteration in case fast convergence
    :param meanFitRatio: a 0-1 scaler, exponential moving average coefficient, used for the initialization of U when using group initialized V
    :param error: difference of cost function for convergence
    :param normW: 1 or 2, normalization method for W used in Laplacian regularization
    :param Alpha: hyper parameter for spatial sparsity
    :param Beta: hyper parameter for Laplacian sparsity
    :param alphaS: internally determined, the coefficient for spatial sparsity based Alpha, data size, K, and gNb
    :param alphaL: internally determined, the coefficient for Laplacian sparsity based Beta, data size, K, and gNb
    :param vxI: flag for using the temporal correlation between nodes (vertex, voxel)
    :param ard: 0 or 1, flat for combining similar clusters
    :param eta: a hyper parameter for the ard regularization term
    :param nRepeat: Any positive integer, the number of repetition to avoid poor initialization

    :param Parallel: False or True, whether to enable parallel computation
    :param Computation_Mode: 'CPU_Numpy', 'CPU_Torch'
    :param N_Thread: positive integers, used for parallel computation

    :param dataPrecision: 'double' or 'single'

    :param outputFormat: 'MAT', 'Both', 'MAT' is to save results in FN.mat and TC.mat for functional networks and time courses respectively. 'Both' is for both matlab format and fMRI input file format

    :param synchronized_view: True or False, whether to synchronize view centers for volume data between gFNs and pFNs
    :param synchronized_colorbar: True or False, whether to synchronize color bar between gFNs and pFNs

    Yuncong Ma, 11/29/2023
    """

    # Check setting
    check_data_type_format(dataType, dataFormat)

    # setup all sub-folders in the pNet result folder
    dir_pnet_dataInput, dir_pnet_FNC, dir_pnet_gFN, dir_pnet_pFN, dir_pnet_QC, dir_pnet_STAT = setup_result_folder(dir_pnet_result)

    # ============== Data Input ============== #
    # setup dataInput
    setup_scan_info(
        dir_pnet_dataInput=dir_pnet_dataInput,
        dataType=dataType, dataFormat=dataFormat,
        file_scan=file_scan, file_subject_ID=file_subject_ID,
        file_subject_folder=file_subject_folder, file_group_ID=file_group_ID,
        Combine_Scan=Combine_Scan
    )
    # setup brain template
    # Volume and surface data types require different inputs to compute the brain template
    if file_Brain_Template is None:
        if dataType == 'Volume':
            setup_brain_template(
                dir_pnet_dataInput,
                dataType=dataType, dataFormat=dataFormat,
                templateFormat=templateFormat,
                file_mask_vol=file_mask_vol, file_overlayImage=file_overlayImage,
                maskValue=maskValue
            )
        elif dataType == 'Surface':
            setup_brain_template(
                dir_pnet_dataInput,
                dataType=dataType, dataFormat=dataFormat,
                templateFormat=templateFormat,
                file_surfL=file_surfL, file_surfR=file_surfR,
                file_maskL=file_maskL, file_maskR=file_maskR,
                maskValue=maskValue,
                file_surfL_inflated=file_surfL_inflated, file_surfR_inflated=file_surfR_inflated
            )
        elif dataType == 'Surface-Volume':
            setup_brain_template(
                dir_pnet_dataInput,
                dataType=dataType, dataFormat=dataFormat,
                templateFormat=templateFormat,
                file_surfL=file_surfL, file_surfR=file_surfR,
                file_maskL=file_maskL, file_maskR=file_maskR,
                file_mask_vol=file_mask_vol, file_overlayImage=file_overlayImage,
                maskValue=maskValue,
                file_surfL_inflated=file_surfL_inflated, file_surfR_inflated=file_surfR_inflated
            )

    else:
        setup_brain_template(dir_pnet_dataInput, file_Brain_Template)
    # ============================================= #

    # ============== FN Computation ============== #
    # setup parameters for FN computation
    setup_NMF_setting(
        dir_pnet_result,
        K=K,
        Combine_Scan=Combine_Scan,
        file_gFN=file_gFN,
        samplingMethod=samplingMethod, sampleSize=sampleSize, nBS=nBS,
        maxIter=maxIter, minIter=minIter, meanFitRatio=meanFitRatio, error=error, normW=normW,
        Alpha=Alpha, Beta=Beta, alphaS=alphaS, alphaL=alphaL,
        vxI=vxI, ard=ard, eta=eta,
        nRepeat=nRepeat,
        Parallel=Parallel, Computation_Mode=Computation_Mode, N_Thread=N_Thread,
        dataPrecision=dataPrecision,
        outputFormat=outputFormat
    )
    # perform FN computation
    if Computation_Mode == 'CPU_Numpy':
        run_FN_Computation(dir_pnet_result)
    elif Computation_Mode == 'CPU_Torch':
        run_FN_Computation_torch(dir_pnet_result)
    # ============================================= #

    # ============== Quality Control ============== #
    # perform quality control
    if Computation_Mode == 'CPU_Numpy':
        run_quality_control(dir_pnet_result)
    elif Computation_Mode == 'CPU_Torch':
        run_quality_control_torch(dir_pnet_result)
    # ============================================= #

    # =============== Visualization =============== #
    setup_Visualization(dir_pnet_result, synchronized_view=synchronized_view, synchronized_colorbar=synchronized_colorbar)
    run_Visualization(dir_pnet_result)
    # ============================================= #


def workflow_simple(dir_pnet_result: str,
                    dataType: str, dataFormat: str,
                    file_scan: str,
                    file_Brain_Template: str,
                    K=17,
                    Combine_Scan=False,
                    file_gFN=None):
    """
    Run the workflow of pNet, including Data Input, FN Computation, Quality Control and Visualization
    This is a minimal version of run_workflow for fast deployment using a single job

    :param dir_pnet_result: directory of the pNet result folder
    :param dataType: 'Surface', 'Volume'
    :param dataFormat: 'HCP Surface (*.cifti, *.mat)', 'MGH Surface (*.mgh)', 'MGZ Surface (*.mgz)', 'Volume (*.nii, *.nii.gz, *.mat)', 'HCP Surface-Volume (*.cifti)', 'HCP Volume (*.cifti)'
    :param file_scan: a txt file that stores directories of all fMRI scans
    :param file_Brain_Template: file directory or content of a brain template file in json format
    :param K: number of FNs
    :param Combine_Scan: False or True, whether to combine multiple scans for the same subject
    :param file_gFN: directory of a precomputed gFN in .mat format

    Yuncong Ma, 11/22/2023
    """

    # setup all sub-folders in the pNet result folder
    dir_pnet_dataInput, dir_pnet_FNC, dir_pnet_gFN, dir_pnet_pFN, dir_pnet_QC, dir_pnet_STAT = setup_result_folder(dir_pnet_result)

    # ============== Data Input ============== #
    # setup dataInput
    setup_scan_info(
        dir_pnet_dataInput=dir_pnet_dataInput,
        dataType=dataType, dataFormat=dataFormat,
        file_scan=file_scan,
        Combine_Scan=Combine_Scan
    )
    # setup brain template
    setup_brain_template(dir_pnet_dataInput, file_Brain_Template)
    # ============================================= #

    # ============== FN Computation ============== #
    # setup parameters for FN computation
    setup_NMF_setting(
        dir_pnet_result,
        K=K,
        Combine_Scan=Combine_Scan,
        file_gFN=file_gFN
    )
    # perform FN computation
    run_FN_Computation_torch(dir_pnet_result)
    # ============================================= #

    # ============== Quality Control ============== #
    # perform quality control
    run_quality_control_torch(dir_pnet_result)
    # ============================================= #

    # =============== Visualization =============== #
    setup_Visualization(dir_pnet_result)
    run_Visualization(dir_pnet_result)
    # ============================================= #


def guide_YN(prompt: str, skip=False, default_value='Y'):
    """
    terminal guidance for choosing yes or no
    Users can type Y, y, Yes, yes, N, n, No, or no

    :param prompt: a string for prompt
    :param skip: False or True to skip the setting
    :param default_value: a default value when skip is enabled
    :return: input_YN, 'Y' or 'N'

    Yuncong Ma, 10/5/2023
    """

    input_YN = None
    while input_YN is None:
        print(prompt + "\n[Y/N]")
        time.sleep(0.1)
        if skip is True:
            print(f"Enter to use default value {default_value}")
            time.sleep(0.1)

        input_YN = input("User Input > ")

        if input_YN is None or len(input_YN) == 0:
            if skip is True:
                input_YN = default_value
            else:
                print('Unknown choice, try again')

        elif input_YN in ('Y', 'y', 'Yes', 'yes'):
            input_YN = 'Y'
        elif input_YN in ('N', 'n', 'No', 'no'):
            input_YN = 'N'
        else:
            print('Unknown choice, try again\nUser Input > ')
            input_YN = None
    return input_YN


def guide_dir(prompt: str):
    """
    terminal guidance for getting a directory

    :param prompt: a string for prompt
    :return: input_dir

    Yuncong Ma, 10/5/2023
    """

    input_dir = None
    while input_dir is None:
        print(prompt)
        time.sleep(0.1)
        input_dir = input("User Input > ")
        if input_dir is None:
            print('Wrong setup, try again\nUser Input >')
    return input_dir


def guide_file(prompt: str, existed=True, extension=None):
    """
    terminal guidance for setting up an existing or a new file

    :param prompt: a string for prompt
    :param existed: True or False, check the existence of the file if True.
    :param extension: None or a str, or a tuple of strings, specifying the file extension
    :return: input_file

    Yuncong Ma, 10/5/2023
    """

    input_file = None
    while input_file is None:
        print("# "+prompt)
        time.sleep(0.1)
        input_file = input("User Input > ")
        if input_file is None:
            print('Wrong setup, try again\nUser Input >')
            continue

        if existed is True and not os.path.isfile(input_file):
            print('Cannot find this file, please try again')
            input_file = None
            continue

        if extension is not None:
            if isinstance(extension, str):
                if not input_file.endswith(extension):
                    print('Please provide a file directory with the required extension')
                    input_file = None
            elif isinstance(extension, tuple):
                flagMatch = False
                for i in enumerate(extension):
                    if input_file.endswith(i):
                        flagMatch = True
                        break
                if flagMatch is False:
                    print('Please provide a file directory with the required extension')
                    input_file = None
            else:
                raise ValueError('The extension needs to be None, a str, or a tuple of strings')

    return input_file


def guide_number(prompt: str, data_type='Int', data_range=None, skip=False, default_value=0):
    """
    terminal guidance for setting up a value

    :param prompt: a string for prompt
    :param data_type: 'Int' or 'Float'
    :param data_range: None or (1, 2)
    :param skip: False or True to skip the setting
    :param default_value: a default value when skip is enabled
    :return: input_value

    Yuncong Ma, 10/5/2023
    """

    input_value = None
    while input_value is None:
        print("# "+prompt)
        time.sleep(0.1)

        if skip is True:
            print(f"Enter to use default value {default_value}")
            time.sleep(0.1)

        input_value = input("User Input > ")

        if input_value is None:
            if skip is True:
                input_value = default_value
                print(f'Set to the default value {default_value}')
                return input_value
            print('Wrong setup, try again\nUser Input >')
        else:
            input_value = float(input_value)
            if data_type == 'Int':
                input_value = int(input_value)
            elif data_type == 'Float':
                input_value = float(input_value)
            if data_range is not None:
                if data_range[0] <= input_value <= data_range[1]:
                    return input_value
                else:
                    print(f'The value should be within {data_range[0]} and {data_range[1]}, try again\nUser Input >')
                    input_value = None
            else:
                return input_value

    return input_value


def guide_choice(prompt: str, list_choice: tuple, skip=False, default_value=None):
    """
    Terminal guidance for getting a directory
    Use sequence ratio to find the most similar option in the list_choice to the user input

    :param prompt: a string for prompt
    :param list_choice: a list of choices
    :param skip: False or True to skip the setting
    :param default_value: None or value
    :return: choice

    Yuncong Ma, 10/5/2023
    """

    choice = None
    while choice is None:
        print("# "+prompt)
        time.sleep(0.1)
        for i in range(len(list_choice)):
            print(str(i+1) + ". " + list_choice[i])
        time.sleep(0.1)
        if skip is False:
            print("Choose by number")
            time.sleep(0.1)
        else:
            print(f"Choose by number or enter to use default {default_value}")
            time.sleep(0.1)

        choice = input("User Input > ")

        if choice in list_choice:
            return choice

        else:
            if choice is not None:
                choice = int(float(choice))
                if 1 <= choice <= len(list_choice):
                    choice = list_choice[choice-1]
                else:
                    choice = None
                    print('Wrong choice, try again\nUser Input >')
            else:
                if skip is True:
                    if default_value is None:
                        choice = list_choice[0]
                    else:
                        choice = default_value
                        print(f'Set to the default value {default_value}')
                    return choice
                print('Wrong choice, try again\nUser Input >')
                choice = None
    return choice


def workflow_guide():
    """
    This is a step-by-step guidance for configuring a workflow of pNet in command line
    It will generate a Python script to run the desired workflow with comments
    Yuncong Ma, 11/7/2023
    """

    print('This is a step-by-step guidance for setting up a workflow of pNet')

    # Setup result folder
    dir_pnet_result = guide_dir('Set up a directory for storing pNet results:')

    # ============== Data Input ============== #
    print('# ============== Data Input ============== # ')
    # setup dataInput
    print('setup dataInput')
    dataType = guide_choice("Choose a data type:", ('Surface', 'Volume', 'Surface-Volume'))
    if dataType == 'Surface':
        dataFormat = guide_choice("Choose a data format:", ('HCP Surface (*.cifti, *.mat)', 'MGH Surface (*.mgh)', 'MGZ Surface (*.mgz)'))
    elif dataType == 'Volume':
        dataFormat = 'Volume (*.nii, *.nii.gz, *.mat)'
        print("data format is automatically set to 'Volume (*.nii, *.nii.gz, *.mat)'")
    else:  # Surface-Volume
        dataFormat = 'HCP Surface-Volume (*.cifti)'
        print("data format is automatically set to 'HCP Surface-Volume (*.cifti)'")
    file_scan = guide_file("Provide a txt formatted file containing all fMRI scans:", existed=True, extension='.txt')
    Choice = guide_YN("Do you have a txt formatted file containing subject ID information for each corresponding scan?", skip=True, default_value='N')
    if Choice == 'Y':
        file_subject_ID = guide_file("Provide a txt formatted file containing subject ID information for each corresponding scan:", existed=True, extension='.txt')
        Choice = guide_YN("Do you have a txt formatted file containing subject folder information for each corresponding scan?")
        if Choice == 'Y':
            file_subject_folder = guide_file("Provide a txt formatted file containing subject folder information for each corresponding scan:", existed=True, extension='.txt')
        else:
            file_subject_folder = None
        Choice = guide_YN("Do you have a txt formatted file containing group ID for each corresponding scan?")
        if Choice == 'Y':
            file_group_ID = guide_file("Provide a txt formatted file containing group ID for each corresponding scan:", existed=True, extension='.txt')
        else:
            file_group_ID = None
    else:
        file_subject_ID = None
        file_subject_folder = None
    Choice = guide_YN("Do you want to concatenate multiple scans for the same subject?", skip=True, default_value='N')
    if Choice == 'Y':
        Combine_Scan = True
    else:
        Combine_Scan = False
    # setup brain template
    file_mask_vol = None
    file_overlayImage = None
    maskValue = 0
    file_surfL = None
    file_surfR = None
    file_maskL = None
    file_maskR = None
    file_surfL_inflated = None
    file_surfR_inflated = None
    Choice = guide_YN("Would you like to select a built-in brain template file?", skip=True, default_value='Y')
    if Choice == 'Y':
        file_Brain_Template = guide_choice("Select a built-in brain template:", ('HCP Surface', 'MNI Volume', 'FreeSurfer_fsaverage5', 'HCP Surface-Volume', 'HCP Subcortical Volume'))
        if file_Brain_Template == 'HCP Surface':
            file_Brain_Template = Brain_Template.file_HCP_surf
        elif file_Brain_Template == 'MNI Volume':
            file_Brain_Template = Brain_Template.file_MNI_vol
        elif file_Brain_Template == 'FreeSurfer_fsaverage5':
            file_Brain_Template = Brain_Template.file_FS_surf
        elif file_Brain_Template == 'HCP Surface-Volume':
            file_Brain_Template = Brain_Template.file_HCP_surf_vol
        elif file_Brain_Template == 'HCP Subcortical Volume':
            file_Brain_Template = Brain_Template.file_HCP_vol
    else:
        Choice = guide_YN("Would you like to select a customized brain template file?")
        if Choice == 'Y':
            file_Brain_Template = guide_file("Set up the directory of the brain template file in json format:", existed=True, extension='.json')
        else:
            file_Brain_Template = None
            # Volume and surface data types require different inputs to compute the brain template
            if dataType == 'Volume':
                templateFormat = guide_choice("Choose a template format:", ('3D Matrix', 'HCP'))
                file_mask_vol = guide_file("Set up the directory of a brain mask:", existed=True, extension=('.mat', '.nii', '.nii.gz'))
                file_overlayImage = guide_file("Set up the directory of a high resolution T1/T2 image as the overlay background:", existed=True, extension=('.mat', '.nii', '.nii.gz'))
                maskValue = guide_number("What is the value used for labeling useful voxels in the brain mask file?", 'Int')
            elif dataType == 'Surface':
                templateFormat = guide_choice("Choose a template format:", ('HCP', 'FreeSurfer'))
                file_surfL = guide_file("Set up the directory of the left hemisphere brain shape file (ex. Conte69.L.inflated.32k_fs_LR.surf.gii):", existed=True, extension='.surf.gii')
                file_surfR = guide_file("Set up the directory of the left hemisphere brain shape file (ex. Conte69.R.inflated.32k_fs_LR.surf.gii):", existed=True, extension='.surf.gii')
                Choice = guide_YN("Would you like to load an inflated brain surface shape?", skip=True, default_value='N')
                if Choice == 'Y':
                    file_surfL = guide_file("Set up the directory of the left hemisphere brain shape file (ex. Conte69.L.very_inflated.32k_fs_LR.surf.gii):", existed=True, extension='.surf.gii')
                    file_surfR = guide_file("Set up the directory of the left hemisphere brain shape file (ex. Conte69.R.very_inflated.32k_fs_LR.surf.gii):", existed=True, extension='.surf.gii')
                else:
                    file_surfL_inflated = None
                    file_surfR_inflated = None
                file_maskL = guide_file("Set up the directory of the left hemisphere brain mask file (ex. medial_wall.L.32k_fs_LR.func.gii):", existed=True, extension='.surf.gii')
                file_maskR = guide_file("Set up the directory of the left hemisphere brain mask file (ex. medial_wall.R.32k_fs_LR.func.gii):", existed=True, extension='.surf.gii')
                maskValue = guide_number("Set up the value used for labeling useful voxels in the brain mask file?", 'Int')
            elif dataType == 'Surface-Volume':
                templateFormat = guide_choice("Choose a template format:", ('HCP', 'FreeSurfer'))
                file_surfL = guide_file("Set up the directory of the left hemisphere brain shape file (ex. Conte69.L.inflated.32k_fs_LR.surf.gii):", existed=True, extension='.surf.gii')
                file_surfR = guide_file("Set up the directory of the left hemisphere brain shape file (ex. Conte69.R.inflated.32k_fs_LR.surf.gii):", existed=True, extension='.surf.gii')
                Choice = guide_YN("Would you like to load an inflated brain surface shape?", skip=True, default_value='N')
                if Choice == 'Y':
                    file_surfL = guide_file("Set up the directory of the left hemisphere brain shape file (ex. Conte69.L.very_inflated.32k_fs_LR.surf.gii):", existed=True, extension='.surf.gii')
                    file_surfR = guide_file("Set up the directory of the left hemisphere brain shape file (ex. Conte69.R.very_inflated.32k_fs_LR.surf.gii):", existed=True, extension='.surf.gii')
                else:
                    file_surfL_inflated = None
                    file_surfR_inflated = None
                file_maskL = guide_file("Set up the directory of the left hemisphere brain mask file (ex. medial_wall.L.32k_fs_LR.func.gii):", existed=True, extension='.surf.gii')
                file_maskR = guide_file("Set up the directory of the left hemisphere brain mask file (ex. medial_wall.R.32k_fs_LR.func.gii):", existed=True, extension='.surf.gii')
                file_mask_vol = guide_file("Set up the directory of a brain mask:", existed=True, extension=('.mat', '.nii', '.nii.gz'))
                file_overlayImage = guide_file("Set up the directory of a high resolution T1/T2 image as the overlay background:", existed=True, extension=('.mat', '.nii', '.nii.gz'))
                maskValue = guide_number("Set up the value used for labeling useful voxels in the brain mask file?", 'Int')
    # ============================================= #

    # ============== FN Computation ============== #
    print('# ============== FN Computation ============== # ')
    # setup parameters for FN computation
    K = guide_number("How many functional networks (default 17)?", 'Int', (2, 1000), skip=True, default_value=17)
    Choice = guide_YN("Do you want to load a precomputed group-level functional networks?", skip=True, default_value='N')
    if Choice == 'Y':
        Compute_gFN = True
        file_gFN = guide_file('Set up the file directory of the precomputed group-level functional networks in matlab format?', existed=True, extension='.mat')
    else:
        Compute_gFN = False
        file_gFN = None
    Choice_simple = guide_YN("Would you like to customize the model parameters for computing personalized functional networks?", skip=True, default_value='N')
    if Choice_simple == 'Y':
        if Compute_gFN is True and file_group_ID is not None:
            samplingMethod = guide_choice("Choose an option for the bootstrap sampling method:", ('Subject', 'Group_Subject'), skip=True)
        else:
            samplingMethod = 'Subject'
        sampleSize = guide_number("Set up the sample size for each bootstrap (ex. 10)?", 'Int', (1, 10000), skip=True, default_value='Automatic')
        nBS = guide_number("Set up the number of bootstrap runs (default 50)?", 'Int', (1, 10000), skip=True, default_value=50)
        maxIter = guide_number("Set up the max iteration number (default 1000)?", 'Int', (1, 1000000), skip=True, default_value=1000)
        minIter = guide_number("Set up the minimum iteration number (default 30)?", 'Int', (1, 1000000), skip=True, default_value=30)
        meanFitRatio = guide_number("Set up the meanFitRatio (default 0.1)?", 'Float', (0.0001, 0.9999), skip=True, default_value=0.1)
        error = guide_number("Set up the error for iteration convergence (default 0.00000001)?", 'Float', (0.0000000000001, 0.1), skip=True, default_value=0.00000001)
        Alpha = guide_number("Set up the Alpha (default 2)?", 'Float', (0.0001, 100000), skip=True, default_value=2)
        Beta = guide_number("Set up the Beta (default 30)?", 'Float', (0.0001, 100000), skip=True, default_value=30)
        nRepeat = guide_number("Set up the number of repeat (default 5)?", 'int', (1, 1000), skip=True, default_value=5)
        Computation_Mode = guide_choice("Choose an option for the computation mode (default CPU_Torch):", ('CPU_Numpy', 'CPU_Torch'), skip=True, default_value='CPU_Torch')
        dataPrecision = guide_choice("Choose an option for data precision (default double):", ('single', 'double'), skip=True, default_value='double')
    outputFormat = guide_YN("Output pFN results matching to the fMRI data format?", skip=True)
    # ============================================= #

    # Generate a python script for the workflow
    print('# ============================ #')
    file_script = guide_file("Setup a Python file directory (*.py) of this customized workflow:", existed=False, extension='.py')

    file_script = open(file_script, 'w')
    print('# Customized Python script for pNet workflow, built at ' + time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())), file=file_script)
    print('# Generated by running python Worflow_guide()', file=file_script)
    print('# To run this python code, use the terminal command line below', file=file_script)
    print(f"# python {file_script.name}", file=file_script)
    print('\n# Load packages', file=file_script)
    print('import pNet\n', file=file_script)

    print('# setup and run a customized workflow', file=file_script)
    if Choice_simple == 'N' and file_Brain_Template is not None:
        print(f"pNet.workflow_simple(", file=file_script)
        print(f"    dir_pnet_result='{dir_pnet_result}',", file=file_script)
        print(f"    dataType='{dataType}',", file=file_script)
        print(f"    file_scan='{file_scan}',", file=file_script)
        print(f"    dataFormat='{dataFormat}',", file=file_script)
        print(f"    file_Brain_Template='{file_Brain_Template}',", file=file_script)
        print(f"    K={K},", file=file_script)
        print(f"    Combine_Scan={Combine_Scan}", file=file_script)  # True or False
        print(")\n", file=file_script)

    else:
        print(f"pNet.workflow(", file=file_script)
        print(f"    dir_pnet_result='{dir_pnet_result}',", file=file_script)
        print(f"    dataType='{dataType}',", file=file_script)
        print(f"    dataFormat='{dataFormat}',", file=file_script)
        print(f"    file_scan='{file_scan}',", file=file_script)
        print(f"    file_subject_ID='{file_subject_ID}',", file=file_script)
        print(f"    file_subject_folder='{file_subject_folder}',", file=file_script)
        print(f"    file_group_ID='{file_group_ID}',", file=file_script)
        if file_Brain_Template is not None:
            print(f"    file_Brain_Template='{file_Brain_Template}',", file=file_script)
            if dataType == 'Surface':
                print(f"    templateFormat='{templateFormat}',", file=file_script)
                print(f"    file_surfL='{file_surfL}',", file=file_script)
                print(f"    file_surfR='{file_surfR}',", file=file_script)
                print(f"    file_maskL='{file_maskL}',", file=file_script)
                print(f"    file_maskR='{file_maskR}',", file=file_script)
                if file_surfL_inflated is not None:
                    print(f"    file_surfL_inflated='{file_surfL_inflated}',", file=file_script)
                    print(f"    file_surfR_inflated='{file_surfR_inflated}',", file=file_script)
            elif dataType == 'Volume':
                print(f"    templateFormat='{templateFormat}',", file=file_script)
                print(f"    file_mask_vol='{file_mask_vol}',", file=file_script)
                print(f"    file_overlayImage='{file_overlayImage}',", file=file_script)
            elif dataType == 'Surface-Volume':
                print(f"    templateFormat='{templateFormat}',", file=file_script)
                print(f"    file_surfL='{file_surfL}',", file=file_script)
                print(f"    file_surfR='{file_surfR}',", file=file_script)
                print(f"    file_maskL='{file_maskL}',", file=file_script)
                print(f"    file_maskR='{file_maskR}',", file=file_script)
                if file_surfL_inflated is not None:
                    print(f"    file_surfL_inflated='{file_surfL_inflated}',", file=file_script)
                    print(f"    file_surfR_inflated='{file_surfR_inflated}',", file=file_script)
                print(f"    file_mask_vol='{file_mask_vol}',", file=file_script)
                print(f"    file_overlayImage='{file_overlayImage}',", file=file_script)
            print(f"    maskValue={maskValue},", file=file_script)
        print(f"    K={K},", file=file_script)
        print(f"    Combine_Scan={Combine_Scan},", file=file_script)  # True or False
        if file_gFN is not None:
            print(f"    file_gFN='{file_gFN}',", file=file_script)
        if Choice_simple == 'Y':
            print(f"    samplingMethod='{samplingMethod}',", file=file_script)
            print(f"    sampleSize={sampleSize},", file=file_script)
            print(f"    nBS={nBS},", file=file_script)
            print(f"    maxIter={maxIter},", file=file_script)
            print(f"    minIter={minIter},", file=file_script)
            print(f"    meanFitRatio={meanFitRatio},", file=file_script)
            print(f"    error={error},", file=file_script)
            print(f"    Alpha={Alpha},", file=file_script)
            print(f"    Beta={Beta},", file=file_script)
            print(f"    nRepeat={nRepeat},", file=file_script)
            print(f"    Computation_Mode='{Computation_Mode}',", file=file_script)
            print(f"    dataPrecision='{dataPrecision}'", file=file_script)
        print(f"    outputFormat='{outputFormat}'", file=file_script)
        print(")\n", file=file_script)

    file_script.close()
    print('Customized workflow script is generated successfully, please open to check the details.')


def workflow_server(dir_pnet_result: str,
                    # data input
                    file_scan: str,
                    dataType='Surface', dataFormat='HCP Surface (*.cifti, *.mat)',
                    file_subject_ID=None, file_subject_folder=None, file_group_ID=None,
                    file_Brain_Template=None,
                    templateFormat='HCP',
                    file_surfL=None, file_surfR=None, file_maskL=None, file_maskR=None,
                    file_mask_vol=None, file_overlayImage=None,
                    maskValue=0,
                    file_surfL_inflated=None, file_surfR_inflated=None,
                    # FN computation
                    K=17, Combine_Scan=False,
                    file_gFN=None,
                    samplingMethod='Subject', sampleSize=10, nBS=50,
                    maxIter=1000, minIter=200, meanFitRatio=0.1, error=1e-8, normW=1,
                    Alpha=2, Beta=30, alphaS=0, alphaL=0, vxI=0, ard=0, eta=0, nRepeat=5,
                    outputFormat='Both',
                    Computation_Mode='CPU_Torch',
                    dataPrecision='double',
                    # visualization
                    synchronized_view=True,
                    synchronized_colorbar=False,
                    # server
                    dir_pnet=None,
                    dir_env=None,
                    dir_python=None,
                    submit_command='qsub -terse -j y',
                    thread_command='-pe threaded ',
                    memory_command='-l h_vmem=',
                    log_command='-o ',
                    computation_resource=dict(memory_bootstrap='50G', thread_bootstrap=4,
                                              memory_fusion='10G', thread_fusion=4,
                                              memory_pFN='10G', thread_pFN=1,
                                              memory_qc='10G', thread_qc=1,
                                              memory_visualization='10G', thread_visualization=1)
                    ):
    """
    Run the workflow of pNet, including Data Input, FN Computation, Quality Control and Visualization
    This function is for running pNet using multiple jobs to facilitate computation in a server environment
    This script can be re-run to restart the desired workflow from where it stops

    :param dir_pnet_result: directory of the pNet result folder
    :param dataType: 'Surface', 'Volume', 'Surface-Volume'
    :param dataFormat: 'HCP Surface (*.cifti, *.mat)', 'MGH Surface (*.mgh)', 'MGZ Surface (*.mgz)', 'Volume (*.nii, *.nii.gz, *.mat)', 'HCP Surface-Volume (*.cifti)', 'HCP Volume (*.cifti)'

    :param file_scan: a txt file that stores directories of all fMRI scans
    :param file_subject_ID: a txt file that store subject ID information corresponding to fMRI scan in file_scan
    :param file_subject_folder: a txt file that store subject folder names corresponding to fMRI scans in file_scan
    :param file_group_ID: a txt file that store group information corresponding to fMRI scan in file_scan

    :param file_Brain_Template: file directory of a brain template file in json format
    :param templateFormat: 'HCP', 'FreeSurfer', '3D Matrix'
    :param file_surfL: file that stores the surface shape information of the left hemisphere, including vertices and faces
    :param file_surfR: file that stores the surface shape information of the right hemisphere, including vertices and faces
    :param file_maskL: file that stores the mask information of the left hemisphere, a 1D 0-1 vector
    :param file_maskR: file that stores the mask information of the right hemisphere, a 1D 0-1 vector
    :param file_surfL_inflated: file that stores the inflated surface shape information of the left hemisphere, including vertices and faces
    :param file_surfR_inflated: file that stores the inflated surface shape information of the right hemisphere, including vertices and faces
    :param file_mask_vol: file of a mask file for volume-based data type
    :param file_overlayImage: file of a background image for visualizing volume-based results
    :param maskValue: 0 or 1, 0 means 0s in mask files are useful vertices, otherwise vice versa. maskValue=0 for medial wall in HCP data, and maskValue=1 for brain masks

    :param K: number of FNs
    :param Combine_Scan: False or True, whether to combine multiple scans for the same subject

    :param file_gFN: None or a directory of a precomputed gFN in .mat format
    :param samplingMethod: 'Subject' or 'Group_Subject'. Uniform sampling based subject ID, or group and then subject ID
    :param sampleSize: number of subjects selected for each bootstrapping run
    :param nBS: number of runs for bootstrap

    :param maxIter: maximum iteration number for multiplicative update
    :param minIter: minimum iteration in case fast convergence
    :param meanFitRatio: a 0-1 scaler, exponential moving average coefficient, used for the initialization of U when using group initialized V
    :param error: difference of cost function for convergence
    :param normW: 1 or 2, normalization method for W used in Laplacian regularization
    :param Alpha: hyper parameter for spatial sparsity
    :param Beta: hyper parameter for Laplacian sparsity
    :param alphaS: internally determined, the coefficient for spatial sparsity based Alpha, data size, K, and gNb
    :param alphaL: internally determined, the coefficient for Laplacian sparsity based Beta, data size, K, and gNb
    :param vxI: flag for using the temporal correlation between nodes (vertex, voxel)
    :param ard: 0 or 1, flat for combining similar clusters
    :param eta: a hyper parameter for the ard regularization term
    :param nRepeat: Any positive integer, the number of repetition to avoid poor initialization
    :param outputFormat: 'MAT', 'Both', 'MAT' is to save results in FN.mat and TC.mat for functional networks and time courses respectively. 'Both' is for both matlab format and fMRI input file format

    :param Computation_Mode: 'CPU_Numpy', 'CPU_Torch'
    :param dataPrecision: 'double' or 'single'

    :param synchronized_view: True or False, whether to synchronize view centers for volume data between gFNs and pFNs
    :param synchronized_colorbar: True or False, whether to synchronize color bar between gFNs and pFNs

    :param dir_pnet: directory of the pNet toolbox
    :param dir_env: directory of the desired virtual environment
    :param dir_python: absolute directory to the python folder, ex. /Users/YuncongMa/.conda/envs/pnet/bin/python
    :param submit_command: command to submit a server job
    :param thread_command: command to setup number of threads for each job
    :param memory_command: command to setup memory allowance for each job
    :param log_command: command to specify the logfile
    :param computation_resource: a dict to specify the number of threads and memory allowance for jobs in each predefined step

    Yuncong Ma, 12/11/2023
    """

    print('Start to run pNet workflow in server mode', flush=True)

    # Check setting
    check_data_type_format(dataType, dataFormat)

    if dir_pnet is None:
        raise ValueError('Require a valid setting for dir_pnet')
    if dir_python is None:
        raise ValueError('Require a valid setting for dir_python')

    # setup all sub-folders in the pNet result folder
    dir_pnet_dataInput, dir_pnet_FNC, dir_pnet_gFN, dir_pnet_pFN, dir_pnet_QC, dir_pnet_STAT = setup_result_folder(dir_pnet_result)

    # Generate a python script for the workflow
    # create script folder
    dir_script = os.path.join(dir_pnet_dataInput, 'Script')
    os.makedirs(dir_script, exist_ok=True)
    file_script = open(os.path.join(dir_pnet_dataInput, 'Script', 'server_job_workflow.py'), 'w')
    print('# Customized Python script for pNet workflow in server mode\n# Use corresponding bash script to submit the job\n# Created on ' + time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())), file=file_script)
    print('# This Python script can be re-run to restart the workflow from where it stops', file=file_script)
    print('\n# Load packages', file=file_script)
    print('# setup and run a customized workflow\n', file=file_script)
    print('import sys\nimport os\n', file=file_script, flush=True)
    print(f"dir_pnet = '{dir_pnet}'", file=file_script, flush=True)
    print(f"sys.path.append(os.path.join(dir_pnet, 'Python'))", file=file_script, flush=True)
    print('import pNet\n', file=file_script, flush=True)

    print("\n# ============== Parameter ============== #", file=file_script)
    print(f"dir_pnet_result = '{dir_pnet_result}'", file=file_script)
    print('\n# data input', file=file_script)
    print(f"dataType = '{dataType}'", file=file_script)
    print(f"dataFormat = '{dataFormat}'", file=file_script)
    print(f"file_scan = '{file_scan}'", file=file_script)
    if file_subject_ID is None:
        print(f"file_subject_ID = None", file=file_script)
    else:
        print(f"file_subject_ID = '{file_subject_ID}'", file=file_script)
    if file_subject_folder is None:
        print(f"file_subject_folder = None", file=file_script)
    else:
        print(f"file_subject_folder = '{file_subject_folder}'", file=file_script)
    if file_group_ID is None:
        print(f"file_group_ID = None", file=file_script)
    else:
        print(f"file_group_ID = '{file_group_ID}'", file=file_script)
    if file_Brain_Template is not None:
        print(f"file_Brain_Template = '{file_Brain_Template}'", file=file_script)
    else:
        if dataType == 'Surface':
            print(f"templateFormat = '{templateFormat}'", file=file_script)
            print(f"file_surfL = '{file_surfL}'", file=file_script)
            print(f"file_surfR = '{file_surfR}'", file=file_script)
            print(f"file_maskL = '{file_maskL}'", file=file_script)
            print(f"file_maskR = '{file_maskR}'", file=file_script)
            if file_surfL_inflated is not None:
                print(f"file_surfL_inflated = '{file_surfL_inflated}'", file=file_script)
                print(f"file_surfR_inflated = '{file_surfR_inflated}'", file=file_script)
        elif dataType == 'Volume':
            print(f"templateFormat = '{templateFormat}'", file=file_script)
            print(f"file_mask_vol='{file_mask_vol}'", file=file_script)
            print(f"file_overlayImage = '{file_overlayImage}'", file=file_script)
        elif dataType == 'Surface-Volume':
            print(f"templateFormat = '{templateFormat}'", file=file_script)
            print(f"file_surfL = '{file_surfL}'", file=file_script)
            print(f"file_surfR = '{file_surfR}'", file=file_script)
            print(f"file_maskL = '{file_maskL}'", file=file_script)
            print(f"file_maskR = '{file_maskR}'", file=file_script)
            if file_surfL_inflated is not None:
                print(f"file_surfL_inflated = '{file_surfL_inflated}'", file=file_script)
                print(f"file_surfR_inflated = '{file_surfR_inflated}'", file=file_script)
            print(f"file_mask_vol = '{file_mask_vol}'", file=file_script)
            print(f"file_overlayImage = '{file_overlayImage}'", file=file_script)
        print(f"maskValue = {maskValue}", file=file_script)
    print('\n# FN computation', file=file_script)
    print(f"K = {K}", file=file_script)
    print(f"Combine_Scan = {Combine_Scan}", file=file_script)  # True or False
    if file_gFN is None:
        print(f"file_gFN = None", file=file_script)
    else:
        print(f"file_gFN = '{file_gFN}'", file=file_script)
    print(f"samplingMethod = '{samplingMethod}'", file=file_script)
    print(f"sampleSize = {sampleSize}", file=file_script)
    print(f"nBS = {nBS}", file=file_script)
    print(f"maxIter = {maxIter}", file=file_script)
    print(f"minIter = {minIter}", file=file_script)
    print(f"meanFitRatio = {meanFitRatio}", file=file_script)
    print(f"error = {error}", file=file_script)
    print(f"normW = {normW}", file=file_script)
    print(f"Alpha = {Alpha}", file=file_script)
    print(f"Beta = {Beta}", file=file_script)
    print(f"alphaS = {alphaS}", file=file_script)
    print(f"alphaL = {alphaL}", file=file_script)
    print(f"vxI = {vxI}", file=file_script)
    print(f"eta = {eta}", file=file_script)
    print(f"ard = {ard}", file=file_script)
    print(f"nRepeat = {nRepeat}", file=file_script)
    print(f"Computation_Mode = '{Computation_Mode}'", file=file_script)
    print(f"dataPrecision = '{dataPrecision}'", file=file_script)
    print(f"outputFormat = '{outputFormat}'", file=file_script)
    print('# visualization', file=file_script)
    print(f"synchronized_view = {synchronized_view}", file=file_script)
    print(f"synchronized_colorbar = {synchronized_colorbar}", file=file_script)
    print('\n# server', file=file_script)
    print(f"dir_env = '{dir_env}'", file=file_script)
    print(f"dir_python = '{dir_python}'", file=file_script)
    print(f"submit_command = '{submit_command}'", file=file_script)
    print(f"thread_command = '{thread_command}'", file=file_script)
    print(f"memory_command = '{memory_command}'", file=file_script)
    print(f"log_command = '{log_command}'", file=file_script)

    print(f"computation_resource = {computation_resource}", file=file_script)

    # main job
    print('\n# Main job\n# The following part is imported from /pNet/Python/Workflow_Server_Template.py', file=file_script)

    file_pnet_workflow_server_template = os.path.join(dir_pnet, 'Python', 'Workflow_Server_Template.py')
    [print(line.replace('\n', ''), file=file_script) for line in open(file_pnet_workflow_server_template, 'r')]

    # =============== Server
    setup_server(
        dir_pnet=dir_pnet,
        dir_env=dir_env,
        dir_pnet_result=dir_pnet_result,
        dir_python=dir_python,
        submit_command=submit_command,
        thread_command=thread_command,
        memory_command=memory_command,
        log_command=log_command,
        computation_resource=computation_resource
    )

    # submit bash job
    submit_bash_job(dir_pnet_result=dir_pnet_result,
                    python_command=None,
                    bashFile=os.path.join(dir_script, 'server_job_workflow.sh'),
                    pythonFile=os.path.join(dir_script, 'server_job_workflow.py'),
                    logFile=os.path.join(dir_script, 'server_job_workflow.log'),
                    memory='10G',
                    n_thread=1,
                    create_python_file=False)
    print('Workflow job is submitted', flush=True)

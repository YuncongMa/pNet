# Yuncong Ma, 9/29/2023
# pNet
# Provide examples of running the whole workflow of pNet

#########################################
# Packages

# Example
from Example import Example

# Module
# This script builds the five modules of pNet
# Functions for modules of pNet
from Data_Input import *
from FN_Computation import *
from FN_Computation_torch import *
from Computation_Environment import *
from Quality_Control import *
from Quality_Control_torch import *
import Levenshtein


def workflow(dir_pnet_result: str,
             file_scan: str,
             dataType='Surface', dataFormat='HCP Surface (*.cifti, *.mat)',
             file_subject_ID=None, file_subject_folder=None, file_group=None,
             file_Brain_Template=None,
             file_surfL=None, file_surfR=None, file_maskL=None, file_maskR=None,
             file_mask_vol=None, file_overlayImage=None,
             maskValue=0,
             file_surfL_inflated=None, file_surfR_inflated=None,
             K=17, Combine_Scan=False,
             Compute_gFN=True, file_gFN=None,
             samplingMethod='Subject', sampleSize=10, nBS=50,
             maxIter=1000, minIter=30, meanFitRatio=0.1, error=1e-6, normW=1,
             Alpha=2, Beta=30, alphaS=0, alphaL=0, vxI=0, ard=0, eta=0, nRepeat=5,
             Parallel=False, Computation_Mode='CPU_Torch', N_Thread=1,
             dataPrecision='double'):
    """
    Run the workflow of pFN, including Data Input, FN Computation, and Quality Control

    :param dir_pnet_result: directory of the pNet result folder
    :param dataType: 'Surface', 'Volume'
    :param dataFormat: 'HCP Surface (*.cifti, *.mat)', 'MGH Surface (*.mgh)', 'MGZ Surface (*.mgz)', 'Volume (*.nii, *.nii.gz, *.mat)'

    :param file_scan: a txt file that stores directories of all fMRI scans
    :param file_subject_ID: a txt file that store subject ID information corresponding to fMRI scan in file_scan
    :param file_subject_folder: a txt file that store subject folder names corresponding to fMRI scans in file_scan
    :param file_group: a txt file that store group information corresponding to fMRI scan in file_scan

    :param file_Brain_Template: file directory of a brain template file in json format
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

    :param Compute_gFN: True or False, whether to compute gFNs from the provided data or load a precomputed gFN set
    :param file_gFN: directory of a precomputed gFN in .mat format
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

    Yuncong Ma, 9/28/2023
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
        file_subject_folder=file_subject_folder, file_group_ID=file_group,
        Combine_Scan=Combine_Scan
    )
    # setup brain template
    # Volume and surface data types require different inputs to compute the brain template
    if file_Brain_Template is None:
        if dataType == 'Volume':
            setup_brain_template(
                dir_pnet_dataInput,
                dataType=dataType, dataFormat=dataFormat,
                file_mask_vol=file_mask_vol, file_overlayImage=file_overlayImage,
                maskValue=maskValue,
                logFile=None
            )
        elif dataType == 'Surface':
            setup_brain_template(
                dir_pnet_dataInput,
                dataType=dataType, dataFormat=dataFormat,
                file_surfL=file_surfL, file_surfR=file_surfR,
                file_maskL=file_maskL, file_maskR=file_maskR,
                maskValue=maskValue,
                file_surfL_inflated=file_surfL_inflated, file_surfR_inflated=file_surfR_inflated,
                logFile=None
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
        Compute_gFN=Compute_gFN, file_gFN=file_gFN,
        samplingMethod=samplingMethod, sampleSize=sampleSize, nBS=nBS,
        maxIter=maxIter, minIter=minIter, meanFitRatio=meanFitRatio, error=error, normW=normW,
        Alpha=Alpha, Beta=Beta, alphaS=alphaS, alphaL=alphaL,
        vxI=vxI, ard=ard, eta=eta,
        nRepeat=nRepeat,
        Parallel=Parallel, Computation_Mode=Computation_Mode, N_Thread=N_Thread,
        dataPrecision=dataPrecision
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


def workflow_simple(dir_pnet_result: str,
                    dataType: str, dataFormat: str,
                    file_scan: str,
                    file_Brain_Template: str,
                    K=17,
                    Combine_Scan=False,
                    Compute_gFN=True,
                    file_gFN=None):
    """
    Run the workflow of pFN, including Data Input, FN Computation, and Quality Control
    This is a minimal version of run_workflow for fast deployment

    :param dir_pnet_result: directory of the pNet result folder
    :param dataType: 'Surface', 'Volume'
    :param dataFormat: 'HCP Surface (*.cifti, *.mat)', 'MGH Surface (*.mgh)', 'MGZ Surface (*.mgz)', 'Volume (*.nii, *.nii.gz, *.mat)'
    :param file_scan: a txt file that stores directories of all fMRI scans
    :param file_Brain_Template: file directory or content of a brain template file in json format
    :param K: number of FNs
    :param Combine_Scan: False or True, whether to combine multiple scans for the same subject
    :param Compute_gFN: True or False, whether to compute gFNs from the provided data or load a precomputed gFN set
    :param file_gFN: directory of a precomputed gFN in .mat format

    Yuncong Ma, 9/28/2023
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
        Compute_gFN=Compute_gFN,
        file_gFN=file_gFN
    )
    # perform FN computation
    run_FN_Computation_torch(dir_pnet_result)
    # ============================================= #

    # ============== Quality Control ============== #
    # perform quality control
    run_quality_control_torch(dir_pnet_result)
    # ============================================= #


def guide_YN(prompt: str):
    """
    terminal guidance for choosing yes or no

    :param prompt: a string for prompt
    :return: input_YN

    Yuncong Ma, 9/29/2023
    """

    input_YN = None
    while input_YN is None:
        input_YN = input(prompt)
        if input_YN is None:
            print('Unknown choice, try again')
        elif input_YN in ('Y', 'y', 'Yes', 'yes'):
            input_YN = 'Y'
        elif input_YN in ('N', 'n', 'No', 'no'):
            input_YN = 'N'
        else:
            print('Unknown choice, try again')
            input_YN = None
    return input_YN


def guide_dir(prompt: str):
    """
    terminal guidance for getting a directory

    :param prompt: a string for prompt
    :return: input_dir

    Yuncong Ma, 9/29/2023
    """

    input_dir = None
    while input_dir is None:
        input_dir = input(prompt)
        if input_dir is None:
            print('Wrong setup, try again')
    return input_dir


def guide_file(prompt: str):
    """
    terminal guidance for setting up a file

    :param prompt: a string for prompt
    :return: input_file

    Yuncong Ma, 9/29/2023
    """

    input_file = None
    while input_file is None:
        input_file = input(prompt)
        if input_file is None:
            print('Wrong setup, try again')
        elif not os.path.isfile(input_file):
            print('Cannot find this file, please try again')
    return input_file


def guide_choice(prompt: str, list_choice: tuple):
    """
    Terminal guidance for getting a directory
    Use sequence ratio to find the most similar option in the list_choice to the user input

    :param prompt: a string for prompt
    :param list_choice: a list of choices
    :return: choice

    Yuncong Ma, 9/29/2023
    """

    choice = None
    while choice is None:
        choice = input(prompt + 'Type in keywords to select an option\n')
        if choice in list_choice:
            return choice

        else:
            if choice is not None:
                score = [Levenshtein.seqratio(x, choice) for i, x in enumerate(list_choice)]
                ps = np.argmax(score)
                if isinstance(ps, np.int64):
                    choice = list_choice[ps]
                    print('Selected the most relevant option: ' + choice)
                else:
                    print('Find more than one potential option matched, try to give more information')
            else:
                print('Wrong setup, try again')
                choice = None
    return choice


def workflow_guide():
    """
    This is a step-by-step guidance for setting up a workflow of pNet in terminal
    It will generate a Python script to run the desired workflow
    Yuncong Ma, 9/29/2023
    """

    print('This is a step-by-step guidance for setting up a workflow of pNet')

    # Setup result folder
    dir_pnet_result = guide_dir('Please type in a directory for storing pNet results:\n')

    # ============== Data Input ============== #
    print('# ============== Data Input ============== # ')
    # setup dataInput
    print('setup dataInput')
    dataType = guide_choice("Choose a data type: 'Surface', 'Volume'\n", ('Surface', 'Volume'))
    dataFormat = guide_choice("Choose a data format: \n"
                              "'HCP Surface (*.cifti, *.mat)', 'MGH Surface (*.mgh)', 'MGZ Surface (*.mgz)', 'Volume (*.nii, *.nii.gz, *.mat)'\n",
                              ('HCP Surface (*.cifti, *.mat)', 'MGH Surface (*.mgh)', 'MGZ Surface (*.mgz)', 'Volume (*.nii, *.nii.gz, *.mat)'))
    file_scan = guide_file("Provide a txt formatted file containing all fMRI scans:\n")
    Choice = guide_YN("Do you have a txt formatted file containing subject ID information for each corresponding scan: \n[Y/N]")
    if Choice == 'Y':
        file_subject_ID = guide_file("Provide a txt formatted file containing subject ID information for each corresponding scan:\n")
    Choice = guide_YN("Do you have a txt formatted file containing subject folder information for each corresponding scan")
    if Choice == 'Y':
        file_subject_folder = guide_file("Provide a txt formatted file containing subject folder information for each corresponding scan:\n[Y/N]\n")
    Choice = guide_YN("Do you want to concatenate multiple scans for the same subject:\n[Y/N]\n")
    if Choice == 'Y':
        Combine_Scan = True
    else:
        Combine_Scan = False
    # setup brain template
    Choice = guide_YN("Would you like to select a built-in brain template file: \n[Y/N]")
    if Choice == 'Y':
        file_Brain_Template = guide_file("Provide a txt formatted file containing all fMRI scans:\n")
    else:
        Choice = guide_YN("Would you like to select a customized brain template file: \n[Y/N]")
        if Choice == 'Y':
            file_Brain_Template = guide_file("What is the directory of the brain template file:\n")
        else:
            # Volume and surface data types require different inputs to compute the brain template
            if dataType == 'Volume':
                file_mask_vol = guide_file("What is the directory of a brain mask:\n")
                file_overlayImage = guide_file("What is the directory of a high resolution T1/T2 image as the overlay background:\n")
                maskValue = guide_number()


    # Volume and surface data types require different inputs to compute the brain template
    if file_Brain_Template is None:
        if dataType == 'Volume':
            setup_brain_template(
                dir_pnet_dataInput,
                dataType=dataType, dataFormat=dataFormat,
                file_mask_vol=file_mask_vol, file_overlayImage=file_overlayImage,
                maskValue=maskValue,
                logFile=None
            )
        elif dataType == 'Surface':
            setup_brain_template(
                dir_pnet_dataInput,
                dataType=dataType, dataFormat=dataFormat,
                file_surfL=file_surfL, file_surfR=file_surfR,
                file_maskL=file_maskL, file_maskR=file_maskR,
                maskValue=maskValue,
                file_surfL_inflated=file_surfL_inflated, file_surfR_inflated=file_surfR_inflated,
                logFile=None
            )

    else:
        setup_brain_template(dir_pnet_dataInput, file_Brain_Template)
    # ============================================= #

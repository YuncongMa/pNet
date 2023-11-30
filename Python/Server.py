# Yuncong Ma, 11/30/2023
# pNet
# Provide functions to submit jobs to server environment


#########################################
# Packages
import sys
from datetime import datetime
import os
from sys import platform

from Data_Input import check_data_type_format, setup_scan_info, setup_brain_template
from Data_Input import setup_result_folder, write_json_setting, load_json_setting
from FN_Computation import *
from FN_Computation_torch import *
from Quality_Control_torch import *
from Visualization import *


#########################################

def setup_server(dir_pnet: str,
                 dir_pnet_result: str,
                 dir_python: str,
                 submit_command='qsub -terse -j y',
                 thread_command='-pe threaded ',
                 memory_command='-l h_vmem=',
                 log_command='-o ',
                 computation_resource=None or dict):

    """
    Setup server environment and commands to submit jobs

    :param dir_pnet: directory of the pNet toolbox
    :param dir_pnet_result: directory of pNet result folder
    :param dir_python: absolute directory to the python folder, ex. /Users/YuncongMa/.conda/envs/pnet/bin/python
    :param submit_command: command to submit a server job
    :param thread_command: command to setup number of threads for each job
    :param memory_command: command to setup memory allowance for each job
    :param log_command: command to specify the logfile
    :param computation_resource: None or a dict which specifies the number of threads and memory for different processes
    :return: None

    Yuncong Ma, 11/28/2023
    """

    dir_pnet_dataInput, _, _, _, _, _ = setup_result_folder(dir_pnet_result)
    setting = {'dir_pnet': dir_pnet, 'dir_python': dir_python, 'submit_command': submit_command, 'thread_command': thread_command, 'memory_command': memory_command, 'log_command': log_command}

    write_json_setting(setting, os.path.join(dir_pnet_dataInput, 'Server_Setting.json'))


def submit_bash_job(dir_pnet_result: str,
                    python_command: str,
                    memory=50, n_thread=4,
                    logFile=None,
                    bashFile=None,
                    pythonFile=None):
    """
    submit a bash job to the desired server environment

    :param dir_pnet_result: directory of pNet result folder
    :param python_command: the Python function to run, with dir_pnet_result as a preset variable
    :param memory: a real number in GB
    :param n_thread: number of threads to use
    :param logFile: full directory of a log file
    :param bashFile: full directory of the bash file to generate
    :param pythonFile: full directory of the python file to generate
    :return: None

    Yuncong Ma, 11/30/2023
    """

    # load server setting
    dir_pnet_dataInput, _, _, _, _, _ = setup_result_folder(dir_pnet_result)
    setting = load_json_setting(os.path.join(dir_pnet_dataInput, 'Server_Setting.json'))
    dir_pnet = setting['dir_pnet']
    dir_python = setting['dir_python']
    submit_command = setting['submit_command']
    thread_command = setting['thread_command']
    memory_command = setting['memory_command']
    log_command = setting['log_command']

    # current date and time
    now = datetime.now()
    date_time = now.strftime("%m/%d/%Y, %H:%M:%S")

    # create a new bash file
    if os.path.isfile(bashFile):
        os.remove(bashFile)

    bashFile = open(bashFile, 'w')

    # header
    print('#!/bin/sh\n', file=bashFile, flush=True)
    print('# This bash script is to run a pNet job in the desired server environment', file=bashFile, flush=True)
    print(f'# created on {date_time}\n', file=bashFile, flush=True)
    print(f'# Use command to submit this job: {submit_command} {thread_command}{n_thread} {memory_command}{memory} {log_command}{logFile} {bashFile.name}\n', file=bashFile, flush=True)
    print(r'echo -e "Start time : `date +%F-%H:%M:%S`\n" ', file=bashFile, flush=True)
    print(f'\n{dir_python} {pythonFile}\n', file=bashFile, flush=True)
    print(r'echo -e "Finished time : `date +%F-%H:%M:%S`\n" ', file=bashFile, flush=True)
    bashFile.close()
    bashFile = bashFile.name

    # create a Python job file
    if os.path.isfile(pythonFile):
        os.remove(pythonFile)

    pythonFile = open(pythonFile, 'w')
    print('# This python file is to run a pNet job', file=pythonFile, flush=True)
    print(f'# created on {date_time}\n', file=pythonFile, flush=True)
    print('import sys\nimport os\n', file=pythonFile, flush=True)
    print(f"dir_pnet = '{dir_pnet}'", file=pythonFile, flush=True)
    print(f"sys.path.append(os.path.join(dir_pnet, 'Python'))", file=pythonFile, flush=True)
    print('import pNet\n', file=pythonFile, flush=True)

    print(f"dir_pnet_result = '{dir_pnet_result}'\n", file=pythonFile, flush=True)
    print(f"{python_command}\n", file=pythonFile, flush=True)

    pythonFile.close()

    # execute a shell command to submit a server job, only for linux based systems
    if platform == "linux":
        os.system(f'{submit_command} {thread_command}{n_thread} {memory_command}{memory} {log_command}{logFile} {bashFile}')


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
                    synchronized_view=True or bool,
                    synchronized_colorbar=True or bool,
                    # server
                    dir_pnet=None,
                    dir_python=None,
                    submit_command='qsub -terse -j y',
                    thread_command='-pe threaded ',
                    memory_command='-l h_vmem=',
                    log_command='-o ',
                    computation_resource=dict(memory_boostrap='50G', memory_fusion='10G', memory_pFN='10G', memory_qc='10G', memory_visualization='10G',
                                              thread_bootstrap=4, thread_fusion=4, thread_pFN=1, thread_qc=1, thread_visualization=1)
                    ):
    """
    Run the workflow of pNet, including Data Input, FN Computation, Quality Control and Visualization
    This function is for running pNet using multiple jobs to facilitate computation in a server environment

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
    :param dir_python: absolute directory to the python folder, ex. /Users/YuncongMa/.conda/envs/pnet/bin/python
    :param submit_command: command to submit a server job
    :param thread_command: command to setup number of threads for each job
    :param memory_command: command to setup memory allowance for each job
    :param log_command: command to specify the logfile
    :param computation_resource: a dict to specify the number of threads and memory allowance for jobs in each predefined step

    Yuncong Ma, 11/30/2023
    """

    print('Start to run pNet workflow in sever mode', flush=True)

    # Check setting
    check_data_type_format(dataType, dataFormat)

    if dir_pnet is None:
        raise ValueError('Require a valid setting for dir_pnet')
    if dir_python is None:
        raise ValueError('Require a valid setting for dir_python')

    # setup all sub-folders in the pNet result folder
    dir_pnet_dataInput, dir_pnet_FNC, dir_pnet_gFN, dir_pnet_pFN, dir_pnet_QC, dir_pnet_STAT = setup_result_folder(dir_pnet_result)

    # ============== Setup ============== #
    # ============== Data Input
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

    # ============== FN Computation
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
        Computation_Mode=Computation_Mode,
        dataPrecision=dataPrecision,
        outputFormat=outputFormat
    )

    # =============== Visualization
    setup_Visualization(dir_pnet_result, synchronized_view=synchronized_view, synchronized_colorbar=synchronized_colorbar)

    # =============== Server
    setup_server(dir_pnet,
                 dir_pnet_result=dir_pnet_result,
                 dir_python=dir_python,
                 submit_command=submit_command,
                 thread_command=thread_command,
                 memory_command=memory_command,
                 log_command=log_command,
                 computation_resource=computation_resource)
    # ================================= #
    print('All setups are finished', flush=True)

    # ============== Run ============== #
    # create script folder
    dir_script = os.path.join(dir_pnet_dataInput, 'Script')
    os.makedirs(dir_script, exist_ok=True)
    # submit bash job
    submit_bash_job(dir_pnet_result=dir_pnet_result,
                    python_command='pNet.workflow_server_main(dir_pnet_result)',
                    bashFile=os.path.join(dir_script, 'Workflow.sh'),
                    pythonFile=os.path.join(dir_script, 'Workflow.py'),
                    logFile=os.path.join(dir_script, 'Workflow_Log.o'),
                    memory='10G',
                    n_thread=1)
    print('Workflow job is submitted', flush=True)
    # ================================= #


def workflow_server_main(dir_pnet_result: str):
    """
    run the main server job for the workflow

    :param dir_pnet_result: directory of the pNet result folder
    :return:

    Yuncong Ma, 11/30/2023
    """

    # ============== FN Computation
    run_FN_Computation_torch_server(dir_pnet_result)

    # ============== Quality Control
    run_quality_control_torch_server(dir_pnet_result)

    # =============== Visualization
    run_Visualization_server(dir_pnet_result)

    return

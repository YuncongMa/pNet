# Yuncong Ma, 2/2/2024
# Test server mode with generated bash and python scripts

# ======= Terminal command ======= #
# Activate a conda env with required packages for running pNet
# $ source activate /cbica/home/mayun/.conda/envs/pnet
# Run the customized python script for a pNet workflow for cluster computation
# $ python /cbica/home/mayun/Projects/NiChart/Script/Workflow_OASIS_server.py

# basic python packages
import os
import sys

# ======= Server mode ======= #
# setup the directory of the pNet toolbox folder
dir_pnet = '/cbica/home/mayun/Projects/NiChart/pNet'
sys.path.append(os.path.join(dir_pnet, 'Python'))
import pNet

# Setup the directory of the Conda Python environment
dir_env = '/Users/yuncongma/anaconda3/envs/pnet'
dir_python = '/Users/yuncongma/anaconda3/envs/pnet/bin/python'

# Setup server commands
submit_command = 'qsub -terse -j y'
thread_command = '-pe threaded '
memory_command = '-l h_vmem='
log_command = '-o '

# ======= Parameters for pNet ======= #
# Setup the result folder
dir_pnet_result = '/Users/yuncongma/Documents/Document/fMRI/Myworks/Nichart/OASIS/Test_FN17_Server'

# Setup directory of the raw or preprocessed data
dir_raw_data = '/Users/yuncongma/Documents/Document/fMRI/Myworks/Nichart/OASIS'
# Setup a local folder which stores txt files for scan information
dir_oasis_local = '/Users/yuncongma/Documents/Document/fMRI/Myworks/Nichart/OASIS'
# A txt file for directory of scans
file_scan = os.path.join(dir_oasis_local, 'Scan_List.txt')
# A txt file for corresponding subject ID information of each scan (prefered to be set)
file_subject_ID = os.path.join(dir_oasis_local, 'Subject_ID.txt')
# A txt file for corresponding subject folder information of each scan (prefered to be set)
# This file helps to combine multiple scans for the same subject, or organize results
file_subject_folder = os.path.join(dir_oasis_local, 'Subject_Folder.txt')
# Whether to combine multiple scans for the same subject
Combine_Scan = False
# A file for pre-existing brain template generated by pNet
file_Brain_Template = os.path.join(dir_oasis_local, 'Brain_Template.json.zip')

# Setup the data type and format
dataType = 'Volume'
dataFormat = 'Volume (*.nii, *.nii.gz, *.mat)'

# Setup the number of functional networks
K = 17
# Setup number of scans loaded for each bootstrap run for estimating gFNs
sampleSize = 100
# Setup number of runs for bootstraps
nBS = 50

# setup computation resource requirement
computation_resource = \
    dict(memory_bootstrap='50G', thread_bootstrap=4,
         memory_fusion='10G', thread_fusion=4,
         memory_pFN='10G', thread_pFN=1,
         memory_qc='10G', thread_qc=1,
         memory_visualization='10G', thread_visualization=1)

# ======= Run pNet workflow ======= #
pNet.workflow_cluster(
    dir_pnet_result=dir_pnet_result,
    dataType=dataType,
    dataFormat=dataFormat,
    file_Brain_Template=file_Brain_Template,
    file_scan=file_scan,
    file_subject_ID=file_subject_ID,
    file_subject_folder=file_subject_folder,
    K=K,
    Combine_Scan=Combine_Scan,
    sampleSize=sampleSize,
    nBS=nBS,
    dir_env=dir_env,
    dir_python=dir_python,
    dir_pnet=dir_pnet,
    submit_command=submit_command,
    thread_command=thread_command,
    memory_command=memory_command,
    log_command=log_command
)
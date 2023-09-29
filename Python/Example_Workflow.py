# Yuncong Ma, 9/29/2023
# Examples of setting up a workflow of pNet
# Those examples provide minimal settings to run a workflow
# Each example has a brief description

# load pNet toolbox
import pNet

# Choose the example
Example = 3


if Example == 1:
    # This example is to perform pNet using surface-based fMRI in HCP format
    # 1. Specify the result folder directory in dir_pnet_result
    # 2. Provide a txt formatted scan list file, such as Scan_List.txt
    # 3. Use a prepared brain template file provided in pNet
    # 4. Choose the desired number of FNs

    # Setup
    dir_pnet_result = '<User>/Test_FN17_HCP_Workflow'  # Change <User> to a desired directory
    file_scan = '/Volumes/Scratch_0/pNet/Example/HCP_Surface/Data/Scan_List.txt'
    file_Brain_Template = pNet.Brain_Template.file_HCP_surf
    K = 17

    # Run pNet workflow
    pNet.workflow_simple(
        dir_pnet_result=dir_pnet_result,
        dataType='Surface',
        dataFormat='HCP Surface (*.cifti, *.mat)',
        file_scan=file_scan,
        file_Brain_Template=file_Brain_Template,
        K=K
    )

elif Example == 2:
    # This example is to perform pNet using volume-based fMRI in NIFTI format, with co-registration to MNI space
    # 1. Specify the result folder directory in dir_pnet_result
    # 2. Provide a txt formatted scan list file, such as Scan_List.txt
    # 3. Use a prepared brain template file provided in pNet
    # 4. Choose the desired number of FNs

    # Setup
    dir_pnet_result = '<User>/Test_FN17_UKBB_Workflow'   # Change <User> to a desired directory
    file_scan = '/Volumes/Scratch_0/pNet/Test/Test_FN17_UKBB/Data_Input/Scan_List.txt'
    file_Brain_Template = pNet.Brain_Template.file_MNI_vol
    K = 17

    # Run pNet workflow
    pNet.workflow_simple(
        dir_pnet_result=dir_pnet_result,
        dataType='Volume',
        dataFormat='Volume (*.nii, *.nii.gz, *.mat)',
        file_scan=file_scan,
        file_Brain_Template=file_Brain_Template,
        K=K
    )

elif Example == 3:
    # This example use step-by-step guidance to set up a workflow of pNet
    pNet.workflow_guide()


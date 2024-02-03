# Yuncong Ma, 10/4/2023
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
    # data type is surface, fixed for this example
    dataType = 'Surface'
    # data format is HCP surface, usually in CIFTI format, but can also be store as a 2D matrix in MAT file, fixed
    dataFormat = 'HCP Surface (*.cifti, *.mat)'
    # directory of the pNet result folder. Change <User> to a desired directory
    dir_pnet_result = '<User>/Test_FN17_HCP_Workflow'
    # a txt file storing directory of each fMRI scan file, required to provide
    file_scan = '/Volumes/Scratch_0/pNet/Example/HCP_Surface/Data/Scan_List.txt'
    # a built-in brain template file, made for the HCP surface data (59412 vertices for cortical gray matter), optional to change
    file_Brain_Template = pNet.Brain_Template.file_HCP_surf
    # number of FNs, can be changed to any positive integer number
    K = 17

    # Run pNet workflow
    pNet.workflow_simple(
        dir_pnet_result=dir_pnet_result,
        dataType=dataType,
        dataFormat=dataFormat,
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
    # data type is volume, fixed for this example
    dataType = 'Volume'
    # data format is NIFTI, which stores a 4D matrix, fixed for this example
    dataFormat = 'Volume (*.nii, *.nii.gz, *.mat)'
    # directory of the pNet result folder. Change <User> to a desired directory
    dir_pnet_result = '<User>/Test_FN17_UKBB_Workflow'
    # a txt file storing directory of each fMRI scan file, required to provide
    file_scan = '/Volumes/Scratch_0/pNet/Test/Test_FN17_UKBB/Data_Input/Scan_List.txt'
    # a built-in brain template file, MNI standard space (2mm isotropic), made for the HCP surface data, optional to change
    file_Brain_Template = pNet.Brain_Template.file_MNI_vol
    # number of FNs, can be changed to any positive integer number
    K = 17

    # Run pNet workflow
    pNet.workflow_simple(
        dir_pnet_result=dir_pnet_result,
        dataType=dataType,
        dataFormat=dataFormat,
        file_scan=file_scan,
        file_Brain_Template=file_Brain_Template,
        K=K
    )

elif Example == 3:
    # This example use step-by-step guidance to set up a workflow of pNet
    pNet.workflow_guide()


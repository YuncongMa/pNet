# Yuncong Ma, 9/24/2023

import pNet

Test = 2

if Test == 1:
    file_scan = '/Volumes/Scratch_0/pNet/Test/Test_FN17_UKBB/Data_Input/Scan_List.txt'
    file_mask_vol = '/Volumes/Scratch_0/pNet/Brain_Template/MNI_Volume/Brain_Mask.mat'
    file_overlayImage = '/Volumes/Scratch_0/pNet/Example/UKBB_Volume/Data/Overlay_Image.mat'

    pNet.run_workflow(dir_pnet_result='/Volumes/Scratch_0/pNet/Test/Test_FN17_UKBB_Workflow',
                      dataType='Volume', dataFormat='Volume (*.nii, *.nii.gz, *.mat)',
                      file_scan=file_scan, scan_info='Automatic',
                      file_mask_vol=file_mask_vol, file_overlayImage=file_overlayImage, maskValue=1,
                      K=17, Combine_Scan=False,
                      Compute_gFN=True, samplingMethod='Subject', sampleSize=5, nBS=10,
                      maxIter=200, minIter=30, meanFitRatio=0.1, error=1e-6, normW=1,
                      Alpha=2, Beta=30, alphaS=0, alphaL=0, vxI=0, ard=0, eta=0, nRepeat=5,
                      Parallel=False, Computation_Mode='CPU_Numpy', N_Thread=1, dataPrecision='double')

elif Test == 2:
    file_scan = '/Volumes/Scratch_0/pNet/Test/Test_FN17_UKBB/Data_Input/Scan_List.txt'
    file_mask_vol = '/Volumes/Scratch_0/pNet/Brain_Template/MNI_Volume/Brain_Mask.mat'
    file_overlayImage = '/Volumes/Scratch_0/pNet/Example/UKBB_Volume/Data/Overlay_Image.mat'

    pNet.run_workflow(dir_pnet_result='/Volumes/Scratch_0/pNet/Test/Test_FN17_UKBB_Workflow',
                      dataType='Volume', dataFormat='Volume (*.nii, *.nii.gz, *.mat)',
                      file_scan=file_scan, scan_info='Automatic',
                      file_mask_vol=file_mask_vol, file_overlayImage=file_overlayImage, maskValue=1,
                      K=17, Combine_Scan=False,
                      Compute_gFN=True, samplingMethod='Subject', sampleSize=5, nBS=10,
                      maxIter=200, minIter=30, meanFitRatio=0.1, error=1e-6, normW=1,
                      Alpha=2, Beta=30, alphaS=0, alphaL=0, vxI=0, ard=0, eta=0, nRepeat=5,
                      Parallel=False, Computation_Mode='CPU_Torch', N_Thread=1, dataPrecision='double')

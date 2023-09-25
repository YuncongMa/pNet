# Yuncong Ma, 9/24/2023

import os
import pNet

Test = 4

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

elif Test == 3:
    file_scan = '/Users/yuncongma/Documents/Document/fMRI/Myworks/pNet/Example/HCP_Surface/Data/Scan_List.txt'
    dir_template = '/Users/yuncongma/Documents/Document/fMRI/Myworks/Template/HCP'
    file_surfL = os.path.join(dir_template, '32k_ConteAtlas_v2', 'Conte69.L.inflated.32k_fs_LR.surf.gii')
    file_surfR = os.path.join(dir_template, '32k_ConteAtlas_v2', 'Conte69.R.inflated.32k_fs_LR.surf.gii')
    file_maskL = os.path.join(dir_template, 'Gordon/Gordon2016Surface_parcellation_distribute-20agwt4',
                              'medial_wall.L.32k_fs_LR.func.gii')
    file_maskR = os.path.join(dir_template, 'Gordon/Gordon2016Surface_parcellation_distribute-20agwt4',
                              'medial_wall.R.32k_fs_LR.func.gii')

    pNet.run_workflow(dir_pnet_result='/Users/yuncongma/Documents/Document/fMRI/Myworks/pNet/Test/Test_FN17_HCP_Workflow',
                      dataType='Surface', dataFormat='HCP Surface (*.cifti, *.mat)',
                      file_scan=file_scan, scan_info='Automatic',
                      file_surfL=file_surfL, file_surfR=file_surfR, file_maskL=file_maskL, file_maskR=file_maskR, maskValue=0,
                      K=17, Combine_Scan=False,
                      Compute_gFN=True, samplingMethod='Subject', sampleSize=5, nBS=10,
                      maxIter=200, minIter=30, meanFitRatio=0.1, error=1e-6, normW=1,
                      Alpha=2, Beta=30, alphaS=0, alphaL=0, vxI=0, ard=0, eta=0, nRepeat=5,
                      Parallel=False, Computation_Mode='CPU_Torch', N_Thread=1, dataPrecision='double')

elif Test == 4:
    file_scan = '/Users/yuncongma/Documents/Document/fMRI/Myworks/pNet/Test/Test_FN17_UKBB/Data_Input/Scan_List.txt'
    file_mask_vol = '/Users/yuncongma/Documents/Document/fMRI/Myworks/pNet/Brain_Template/MNI_Volume/Brain_Mask.mat'
    file_overlayImage = '/Users/yuncongma/Documents/Document/fMRI/Myworks/pNet/Example/UKBB_Volume/Data/Overlay_Image.mat'

    pNet.run_workflow(dir_pnet_result='/Users/yuncongma/Documents/Document/fMRI/Myworks/pNet/Test/Test_FN17_UKBB_Workflow',
                      dataType='Volume', dataFormat='Volume (*.nii, *.nii.gz, *.mat)',
                      file_scan=file_scan, scan_info='Automatic',
                      file_mask_vol=file_mask_vol, file_overlayImage=file_overlayImage, maskValue=1,
                      K=17, Combine_Scan=False,
                      Compute_gFN=True, samplingMethod='Subject', sampleSize=5, nBS=10,
                      maxIter=200, minIter=30, meanFitRatio=0.1, error=1e-6, normW=1,
                      Alpha=2, Beta=30, alphaS=0, alphaL=0, vxI=0, ard=0, eta=0, nRepeat=5,
                      Parallel=False, Computation_Mode='CPU_Torch', N_Thread=1, dataPrecision='double')

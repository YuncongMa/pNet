# Yuncong Ma, 10/7/2023
# A script to generate built-in brain template files
# Some data are missing due to data license, but can be downloaded online

import pNet
import os


# =============== HCP Surface =============== #
# Use shape file
dir_template = os.path.join(pNet.dir_brain_template, 'HCP')
file_surfL = os.path.join(dir_template, '32k_ConteAtlas_v2', 'Conte69.L.inflated.32k_fs_LR.surf.gii')
file_surfR = os.path.join(dir_template, '32k_ConteAtlas_v2', 'Conte69.R.inflated.32k_fs_LR.surf.gii')
file_maskL = os.path.join(dir_template, 'Gordon/Gordon2016Surface_parcellation_distribute-20agwt4', 'medial_wall.L.32k_fs_LR.func.gii')
file_maskR = os.path.join(dir_template, 'Gordon/Gordon2016Surface_parcellation_distribute-20agwt4', 'medial_wall.R.32k_fs_LR.func.gii')
file_surfL_inflated = os.path.join(dir_template, '32k_ConteAtlas_v2', 'Conte69.L.very_inflated.32k_fs_LR.surf.gii')
file_surfR_inflated = os.path.join(dir_template, '32k_ConteAtlas_v2', 'Conte69.L.very_inflated.32k_fs_LR.surf.gii')
#
pNet.setup_brain_template(os.path.join(pNet.dir_brain_template, 'HCP_Surface'),
                          dataType='Surface', dataFormat='HCP Surface (*.cifti, *.mat)',
                          file_surfL=file_surfL, file_surfR=file_surfR,
                          file_maskL=file_maskL, file_maskR=file_maskR, maskValue=0,
                          file_surfL_inflated=file_surfL_inflated, file_surfR_inflated=file_surfR_inflated,
                          logFile=os.path.join(pNet.dir_brain_template, 'HCP_Surface', 'Brain_Template.log')
                          )


# =============== HCP Surface-Volume =============== #
# Use a real HCP fMRI data to generate the volume parts
pNet.setup_cifti_volume(os.path.join(pNet.dir_example, 'HCP_Surface-Volume/Data/100206/rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii'),
                        os.path.join(pNet.dir_brain_template, 'HCP_Surface_Volume/CIFTI_Volume.nii.gz'))
# Use shape files to generate the surface part and combine the previously generated volume parts
pNet.setup_brain_template(os.path.join(pNet.dir_brain_template, 'HCP_Surface_Volume'),
                          dataType='Surface-Volume', dataFormat='HCP Surface-Volume (*.cifti)',
                          file_surfL=file_surfL, file_surfR=file_surfR,
                          file_maskL=file_maskL, file_maskR=file_maskR, maskValue=0,
                          file_surfL_inflated=file_surfL_inflated, file_surfR_inflated=file_surfR_inflated,
                          file_mask_vol=os.path.join(pNet.dir_brain_template, 'HCP_Surface_Volume/CIFTI_Volume.nii.gz'),
                          file_overlayImage=os.path.join(pNet.dir_brain_template, 'HCP_Surface_Volume/T1.nii.gz'),
                          logFile=os.path.join(pNet.dir_brain_template, 'HCP_Surface_Volume/Brain_Template.log')
                          )

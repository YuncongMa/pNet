# Yuncong Ma, 9/7/2023
# Test gNb computation


#########################################
# Packages
import nibabel as nib
import numpy as np
import scipy
import scipy.io as sio
import os
import re
import time

import pNet

dataType = 'Surface'

if dataType == 'Surface':
    # Directory of template files
    dir_template = '/Users/yuncongma/Documents/Document/fMRI/Myworks/HCP'

    file_surfL = os.path.join(dir_template, '32k_ConteAtlas_v2', 'Conte69.L.inflated.32k_fs_LR.surf.gii')
    file_surfR = os.path.join(dir_template, '32k_ConteAtlas_v2', 'Conte69.R.inflated.32k_fs_LR.surf.gii')
    file_maskL = os.path.join(dir_template, 'Gordon/Gordon2016Surface_parcellation_distribute-20agwt4', 'medial_wall.L.32k_fs_LR.func.gii')
    file_maskR = os.path.join(dir_template, 'Gordon/Gordon2016Surface_parcellation_distribute-20agwt4', 'medial_wall.R.32k_fs_LR.func.gii')

    Brain_Surface = pNet.compute_Brain_Surface(file_surfL, file_surfR, file_maskL, file_maskR, maskValue=0, dataType='Surface', dataFormat='HCP Surface (*.cifti, *.mat)')

    sio.savemat('Brain_Surface.mat', {'Brain_Surface': Brain_Surface})

    gNb = pNet.compute_gNb(Brain_Surface)

else:
    Brain_Mask = pNet.load_matlab_single_array('/Users/yuncongma/Documents/Document/fMRI/Myworks/pNet/Brain_Template/MNI_Volume/Brain_Mask.mat')

    Brain_Template = {'Data_Type': 'Volume', 'Data_Format': 'Volume (*.nii, *.nii.gz, *.mat)', 'Mask': Brain_Mask, 'Overlay_Image': ()}
    gNb = pNet.compute_gNb(Brain_Template)

# Output
sio.savemat('gNb.mat', {'gNb': gNb})





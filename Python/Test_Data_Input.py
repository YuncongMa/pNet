# Yuncong Ma, 9/12/2023
# Test the data input module

#########################################
# Packages
import nibabel as nib
import numpy as np
import scipy.io as sio
import os
import re
from scipy.io import savemat
from scipy.sparse import spdiags, lil_matrix
from scipy.spatial import distance
import time

import pNet



Brain_Mask = pNet.load_matlab_single_array('/Users/yuncongma/Documents/Document/fMRI/Myworks/pNet/Brain_Template/MNI_Volume/Brain_Mask.mat')
scan_data = pNet.load_fmri_scan('/Users/yuncongma/Documents/Document/fMRI/Myworks/pNet/Example/UKBB_Volume/Data/1054816/filtered_func_data_clean.nii.gz',
                                dataType='Volume', dataFormat='Volume (*.nii, *.nii.gz, *.mat)', Reshape=True, Brain_Mask=Brain_Mask)
print(scan_data.shape)

sio.savemat('scan.mat', {'scan_data': scan_data})






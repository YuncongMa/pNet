# Yuncong Ma, 9/25/2023

import numpy as np
import nibabel as nib
import scipy
# import pNet
import json

import pNet

gNb = pNet.load_matlab_single_array('/Users/yuncongma/Documents/Document/fMRI/Myworks/pNet/Test/Test_FN17_HCP_Workflow/FN_Computation/gNb.mat')

print(np.unique(gNb[:, 0]).shape[0])
print(np.max(gNb[:, 0]))

print(pNet.Brain_Template.dir_HCP_surf)

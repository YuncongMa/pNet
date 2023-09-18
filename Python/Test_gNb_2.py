# Yuncong Ma, 9/6/2023
# Test gNb computation


#########################################
# Packages
import nibabel as nib
import numpy as np
import scipy.io as sio
import os
import re
import time

import pNet


# Directory of template files
dir_template = '/Users/yuncongma/Documents/Document/fMRI/Myworks/HCP'

shapeL = nib.load(os.path.join(dir_template, '32k_ConteAtlas_v2', 'Conte69.L.inflated.32k_fs_LR.surf.gii'))
shapeR = nib.load(os.path.join(dir_template, '32k_ConteAtlas_v2', 'Conte69.R.inflated.32k_fs_LR.surf.gii'))
MW_L = nib.load(os.path.join(dir_template, 'Gordon/Gordon2016Surface_parcellation_distribute-20agwt4', 'medial_wall.L.32k_fs_LR.func.gii'))
MW_R = nib.load(os.path.join(dir_template, 'Gordon/Gordon2016Surface_parcellation_distribute-20agwt4', 'medial_wall.R.32k_fs_LR.func.gii'))

# Create Brain_Surface
# Index starts from 1
Brain_Surface = {'Version': 'HCP', 'Shape': {'L': {'vertices': [], 'faces': []}, 'R': {'vertices': [], 'faces': []}}, 'MW': {'L': [], 'R': []}}
# Surface shape
Brain_Surface['Shape']['L']['vertices'], Brain_Surface['Shape']['L']['faces'] = shapeL.darrays[0].data + int(1), shapeL.darrays[1].data + int(1)
Brain_Surface['Shape']['R']['vertices'], Brain_Surface['Shape']['R']['faces'] = shapeR.darrays[0].data + int(1), shapeR.darrays[1].data + int(1)
# Index for medial wall
Brain_Surface['MW']['L'] = MW_L.darrays[0].data
Brain_Surface['MW']['R'] = MW_R.darrays[0].data

sio.savemat('Brain_Surface.mat', {'Brain_Surface': Brain_Surface})

# Create gNb using matrix format
# Exclude the medial wall
gNb_L = list()
gNb_R = np.zeros([0], dtype=np.int64)
# Number of vertices
Nv_L = Brain_Surface['Shape']['L']['vertices'].shape[0]  # left hemisphere
Nv_R = Brain_Surface['Shape']['R']['vertices'].shape[0]
# Number of faces
Nf_L = Brain_Surface['Shape']['L']['faces'].shape[0]  # left hemisphere
Nf_R = Brain_Surface['Shape']['R']['faces'].shape[0]
# Use index starting from 1
vL = np.sort(np.where(Brain_Surface['MW']['L'] == 0)[0]) + int(1)  # left hemisphere
vR = np.sort(np.where(Brain_Surface['MW']['R'] == 0)[0]) + int(1)

# Map the index from all vertices to useful ones
mapL = np.zeros((len(vL), 2), dtype=np.int64)
mapL[:, 0] = vL
mapL[:, 1] = range(1, 1+len(vL))
#for i in range(len(vL)):
#    gNb_L[gNb_L == mapL[i, 0]] = mapL[i, 1]

# Get gNb of all useful vertices in the left hemisphere
for i in range(1, 1+Nv_L):
    if Brain_Surface['MW']['L'][i-1] == 1:
        continue

    ps = np.where(np.sum(Brain_Surface['Shape']['L']['faces'] == i, axis=0))
    temp = np.intersect1d(Brain_Surface['Shape']['L']['faces'][ps, :].flatten(), vL)
    temp = np.setdiff1d(temp, i)
    gNb_L.append(temp)
    if i==10:
        break
    # gNb_L = np.append(gNb_L, temp)
print(gNb_L)
sio.savemat('gNb_L.mat', {'gNb_L': gNb_L})
# Right hemisphere



# right hemisphere
mapR = np.zeros((len(vR), 2), dtype=np.int64)
mapR[:, 0] = vR
mapR[:, 1] = range(1, 1+len(vR))
for i in range(len(vR)):
    gNb_R[gNb_R == mapR[i, 0]] = mapR[i, 1]
#gNb_R = mapR[mapR[(gNb_R.flatten() - 1).astype(int), 0], 1]
#gNb_R = np.reshape(gNb_R, (np.round(len(gNb_R)/2), 2))

# Combine two hemispheres
# Shift index in right hemisphere by the number of useful vertices in left hemisphere
gNb = np.concatenate((gNb_L, gNb_R + len(vL)), axis=0)

# Output
sio.savemat('gNb.mat', {'gNb': gNb})





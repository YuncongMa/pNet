# Yuncong Ma, 9/14/2023
# Test the workflow of pFN computation using UKBB volume-based data


#########################################
# Packages
import nibabel as nib
import numpy as np
import scipy
import torch
import scipy.io as sio
import os
import re
import time

import pNet


# Prepare work directory in pNet result folder
dir_pnet = '/Volumes/Scratch_0/pNet'
dir_pnet_result = os.path.join(dir_pnet, 'Test/Test_FN17_UKBB')
dir_pnet_dataInput, dir_pnet_FNC, dir_pnet_gFN, dir_pnet_pFN, dir_pnet_QC, dir_pnet_STAT = pNet.setup_result_folder(dir_pnet_result)
dir_pnet_BS = os.path.join(dir_pnet_FNC, 'BootStrapping')
if not os.path.exists(dir_pnet_BS):
    os.makedirs(dir_pnet_BS)

# Prepare dataset: Scan_List, Subject_ID, Subject_Folder
dir_raw_data = os.path.join(dir_pnet, 'Example/UKBB_Volume/Data')
file_scan = os.path.join(dir_pnet_dataInput, 'Scan_List.txt')
file_subject_ID = os.path.join(dir_pnet_dataInput, 'Subject_ID.txt')
file_subject_folder = os.path.join(dir_pnet_dataInput, 'Subject_Folder.txt')
FID = open(file_scan, 'w')
N_Scan = 0
for root, dirs, files in os.walk(dir_raw_data):
    for file in files:
        if file.endswith('.nii.gz'):
            print(os.path.join(root, file), file=FID)
            N_Scan += 1
FID.close()
# Do not use readline to avoid the maximum txt file reading limit
list_scan = [line.replace('\n', '') for line in open(file_scan, 'r')]

FID = open(file_subject_ID, 'w')
for i in range(N_Scan):
    temp = re.findall(r'\d+', list_scan[i])
    subject_ID = temp[-1]
    print(subject_ID, file=FID)
FID.close()
list_subject_ID = [line.replace('\n', '') for line in open(file_subject_ID, 'r')]

FID = open(file_subject_folder, 'w')
for i in range(N_Scan):
    temp = re.findall(r'\d+', list_scan[i])
    subject_folder = os.path.join(temp[-1])
    print(subject_folder, file=FID)
FID.close()
list_subject_folder = [line.replace('\n', '') for line in open(file_subject_folder, 'r')]

# Create bootstrapped scan lists
BS = 3
N_BS = 10
pNet.bootstrap_scan(dir_pnet_BS, file_scan, file_subject_ID, file_subject_folder, BS=BS, N_BS=N_BS)

# Setup primary parameter
setting = {'K': 17, 'Data_Type': 'Volume', 'Data_Format': 'Volume (*.nii, *.nii.gz, *.mat)'}
Brain_Mask = pNet.load_matlab_single_array(os.path.join(dir_pnet, 'Brain_Template/MNI_Volume/Brain_Mask.mat'))
Brain_Template = {'Data_Type': 'Volume', 'Data_Format': 'Volume (*.nii, *.nii.gz, *.mat)', 'Brain_Mask': Brain_Mask}

# Additional intermediate parameters
gNb = pNet.compute_gNb(Brain_Template)
scipy.io.savemat(os.path.join(dir_pnet_FNC, 'gNb.mat'), {'gNb': gNb})

K = setting['K']
Data_Type = setting['Data_Type']
Data_Format = setting['Data_Format']
for i in range(1, N_BS+1):
    print('Iter ' + str(i) + '')
    logFile = os.path.join(dir_pnet_BS, str(i), 'Log.log')
    file_scan_list = os.path.join(dir_pnet_BS, str(i), 'Scan_List.txt')

    print('Load data')
    Data = pNet.load_fmri_scan(file_scan_list, Data_Type, Data_Format, Reshape=True, Brain_Mask=Brain_Mask, Normalization='vp-vmax', logFile=logFile)
    print(Data.shape)

    print('gFN_NMF')
    FN_BS = pNet.gFN_NMF_torch(Data, K, gNb, logFile=logFile, maxIter=500)
    FN_BS = pNet.reshape_FN(FN_BS.numpy(), dataType=Data_Type, Brain_Mask=Brain_Mask)
    sio.savemat(os.path.join(dir_pnet_BS, str(i), 'FN.mat'), {"FN": FN_BS})

FN_BS = np.empty(N_BS, dtype=torch.Tensor)
for i in range(1, N_BS+1):
    FN_BS[i-1] = torch.tensor(pNet.reshape_fmri_data(pNet.load_matlab_single_array(os.path.join(dir_pnet_BS, str(i), 'FN.mat')), dataType=Data_Type, Brain_Mask=Brain_Mask))
gFN_BS = torch.concatenate(FN_BS, dim=1)
print(gFN_BS.shape)

print('gFN_Fusion')
logFile = os.path.join(dir_pnet_gFN, 'Log.log')
gFN = pNet.gFN_fusion_NCut_torch(gFN_BS, K, logFile=logFile)
gFN = pNet.reshape_FN(gFN.numpy(), dataType=Data_Type, Brain_Mask=Brain_Mask)
print(gFN.shape)
sio.savemat(os.path.join(dir_pnet_gFN, 'FN.mat'), {"FN": gFN})

gFN = pNet.reshape_FN(gFN.numpy(), dataType=Data_Type, Brain_Mask=Brain_Mask)
print('\npFN_NMF')
for i in range(1, N_Scan+1):
    print(' Processing ' + str(i))
    dir_pnet_pFN_indv = os.path.join(dir_pnet_pFN, list_subject_folder[i-1])
    if not os.path.exists(dir_pnet_pFN_indv):
        os.makedirs(dir_pnet_pFN_indv)
    FID = open(os.path.join(dir_pnet_pFN_indv, 'Scan_List.txt'), 'w')
    print(list_scan[i-1], file=FID)
    FID.close()
    logFile = os.path.join(dir_pnet_pFN_indv, 'Log.log')
    Data = pNet.load_fmri_scan(list_scan[i-1], Data_Type, Data_Format, logFile=logFile)
    _, pFN = pNet.pFN_NMF_torch(Data, gFN, gNb, maxIter=1000, logFile=logFile)
    pFN = pNet.reshape_FN(pFN.numpy(), dataType=Data_Type, Brain_Mask=Brain_Mask)
    sio.savemat(os.path.join(dir_pnet_pFN_indv, 'FN.mat'), {"pFN": pFN})

# Yuncong Ma, 9/9/2023
# Test the workflow of pFN computation


#########################################
# Packages
import nibabel as nib
import numpy as np
import torch
import scipy.io as sio
import os
import re
import time

import pNet

# Prepare dataset: Scan_List, Subject_ID, Subject_Folder
dir_raw_data = '/Volumes/Scratch_0/Scratch/HCP_XCPD'
file_raw_data = os.path.join(dir_raw_data, 'Scan_List.txt')
file_subject_ID = os.path.join(dir_raw_data, 'Subject_ID.txt')
file_subject_folder = os.path.join(dir_raw_data, 'Subject_Folder.txt')
FID = open(file_raw_data, 'w')
N_Scan = 0
for root, dirs, files in os.walk(dir_raw_data):
    for file in files:
        if file.endswith('dtseries.nii'):
            print(os.path.join(root, file), file=FID)
            N_Scan += 1
FID.close()
# Do not use readline to avoid the maximum txt file reading limit
list_raw_data = [line.replace('\n', '') for line in open(file_raw_data, 'r')]

FID = open(file_subject_ID, 'w')
for i in range(N_Scan):
    temp = re.findall(r'\d+', list_raw_data[i])
    subject_ID = temp[-3]
    print(subject_ID, file=FID)
FID.close()
list_subject_ID = [line.replace('\n', '') for line in open(file_subject_ID, 'r')]

FID = open(file_subject_folder, 'w')
for i in range(N_Scan):
    temp = re.findall(r'\d+', list_raw_data[i])
    subject_folder = os.path.join(temp[-3], 'REST' + temp[-2])
    print(subject_folder, file=FID)
FID.close()
list_subject_folder = [line.replace('\n', '') for line in open(file_subject_folder, 'r')]


# Prepare work directory in pNet result folder
dir_pnet_result = '/Volumes/Scratch_0/Scratch/Test_FN17_HCP_XCPD'
if not os.path.exists(dir_pnet_result):
    os.makedirs(dir_pnet_result)
dir_pnet_BS = os.path.join(dir_pnet_result, 'FN_Computation', 'BootStrapping')
if not os.path.exists(dir_pnet_BS):
    os.makedirs(dir_pnet_BS)
dir_pnet_gFN = os.path.join(dir_pnet_result, 'Group_FN')
if not os.path.exists(dir_pnet_gFN):
    os.makedirs(dir_pnet_gFN)
dir_pnet_pFN = os.path.join(dir_pnet_result, 'Personalized_FN')
if not os.path.exists(dir_pnet_pFN):
    os.makedirs(dir_pnet_pFN)

# Create bootstrapped scan lists
N_BS = 10
BS = 20
for i in range(1, N_BS+1):
    if not os.path.exists(os.path.join(dir_pnet_BS, str(i))):
        os.makedirs(os.path.join(dir_pnet_BS, str(i)))
    FID = open(os.path.join(dir_pnet_BS, str(i), 'Scan_List.txt'), 'w')
    ps = np.sort(np.random.choice(N_Scan, BS, replace=False))
    for j in ps:
        print(list_raw_data[j], file=FID)
    FID.close()


# Setup primary parameter
setting = {'K': 17, 'Data_Type': 'Surface', 'Data_Format': 'HCP Surface (*.cifti, *.mat)'}
gNb = pNet.load_matlab_single_array(os.path.join(dir_pnet_result, 'FN_Computation', 'gNb.mat'))

K = setting['K']
Data_Type = setting['Data_Type']
Data_Format = setting['Data_Format']
for i in ():  # range(1, N_BS+1):
    print('Iter ' + str(i) + '')
    logFile = os.path.join(dir_pnet_BS, str(i), 'Log.log')
    file_scan_list = os.path.join(dir_pnet_BS, str(i), 'Scan_List.txt')

    print('Load data')
    Data = pNet.load_fmri_scan(file_scan_list, Data_Type, Data_Format, Normalization='vp-vmax', logFile=logFile)
    print(Data.shape)

    print('gFN_NMF')
    FN_BS = pNet.gFN_NMF_torch(Data, K, gNb, logFile=logFile, maxIter=1000)
    sio.savemat(os.path.join(dir_pnet_BS, str(i), 'FN.mat'), {"FN": FN_BS.numpy()})

FN_BS = [torch.tensor(pNet.load_matlab_single_array(os.path.join(dir_pnet_BS, str(i), 'FN.mat'))) for i in range(1, N_BS+1)]
gFN_BS = torch.concatenate(FN_BS, dim=1)
print(gFN_BS.shape)

print('gFN_Fusion')
logFile = os.path.join(dir_pnet_gFN, 'Log.log')
gFN = pNet.gFN_fusion_NCut_torch(gFN_BS, K, logFile=logFile)
print(gFN.shape)
sio.savemat(os.path.join(dir_pnet_gFN, 'FN.mat'), {"FN": gFN.numpy()})


print('\npFN_NMF')
for i in range(1, N_Scan+1):
    print(' Processing ' + str(i))
    dir_pnet_pFN_indv = os.path.join(dir_pnet_pFN, list_subject_folder[i-1])
    if not os.path.exists(dir_pnet_pFN_indv):
        os.makedirs(dir_pnet_pFN_indv)
    FID = open(os.path.join(dir_pnet_pFN_indv, 'Scan_List.txt'), 'w')
    print(list_raw_data[i-1], file=FID)
    FID.close()
    logFile = os.path.join(dir_pnet_pFN_indv, 'Log.log')
    Data = pNet.load_fmri_scan(list_raw_data[i-1], Data_Type, Data_Format, logFile=logFile)
    _, pFN = pNet.pFN_NMF_torch(Data, gFN, gNb, maxIter=1000, logFile=logFile)
    sio.savemat(os.path.join(dir_pnet_pFN_indv, 'FN.mat'), {"pFN": pFN.numpy()})

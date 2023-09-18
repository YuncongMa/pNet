# Yuncong Ma, 9/10/2023
# Test the quality control

import os
import re
import numpy as np
import scipy
import time
import sys
import pNet


# Prepare work directory in pNet result folder
dir_pnet_result = '/Volumes/Scratch_0/Scratch/Test_FN17_HCP_XCPD'
dir_pnet_dataInput = os.path.join(dir_pnet_result, 'Data_Input')
dir_pnet_gFN = os.path.join(dir_pnet_result, 'Group_FN')
dir_pnet_pFN = os.path.join(dir_pnet_result, 'Personalized_FN')
dir_pnet_QC = os.path.join(dir_pnet_result, 'Quality_Control')
if not os.path.exists(dir_pnet_QC):
    os.makedirs(dir_pnet_QC)

setting = {'K': 17, 'Data_Type': 'Surface', 'Data_Format': 'HCP Surface (*.cifti, *.mat)'}

K = setting['K']
Data_Type = setting['Data_Type']
Data_Format = setting['Data_Format']

# Setup parameters
combineFlag = 0

file_scan = os.path.join(dir_pnet_dataInput, 'Scan_List.txt')
file_subject_ID = os.path.join(dir_pnet_dataInput, 'Subject_ID.txt')
file_subject_folder = os.path.join(dir_pnet_dataInput, 'Subject_Folder.txt')

list_scan = np.array([line.replace('\n', '') for line in open(file_scan, 'r')])
list_subject_ID = np.array([line.replace('\n', '') for line in open(file_subject_ID, 'r')])
subject_ID_unique = np.unique(list_subject_ID)
N_Subject = subject_ID_unique.shape[0]
list_subject_folder = np.array([line.replace('\n', '') for line in open(file_subject_folder, 'r')])
subject_folder_unique = np.unique(list_subject_folder)

# Load gFNs
gFN = pNet.load_matlab_single_array(os.path.join(dir_pnet_gFN, 'FN.mat'))  # [dim_space, K]

#
dataPrecision = 'double'
np_float, np_eps = pNet.set_data_precision(dataPrecision)

file_Final_Report = open(os.path.join(dir_pnet_QC, 'Final_Report.txt'), 'w')
flag_QC = 0
print('\nStart QC at '+time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))+'\n', file=file_Final_Report, flush=True)
if combineFlag == 0:
    N_pFN = list_scan.shape[0]
    for i in range(N_pFN):
        print('  Processing ' + str(i+1))
        dir_pFN_indv = os.path.join(dir_pnet_pFN, list_subject_folder[i])
        pFN = pNet.load_matlab_single_array(os.path.join(dir_pFN_indv, 'FN.mat'))

        # Spatial correspondence
        Spatial_Correspondence = pNet.mat_corr(gFN, pFN, dataPrecision=dataPrecision)
        Delta_Spatial_Correspondence = np.diag(Spatial_Correspondence) - np.max(Spatial_Correspondence - np.diag(2 * np.ones(K)), axis=0)

        # Miss match between gFNs and pFNs
        # Index starts from 1
        if np.min(Delta_Spatial_Correspondence) >= 0:
            Miss_Match = np.empty((0, ))
        else:
            ps = np.where(Delta_Spatial_Correspondence < 0)[0]
            ps2 = np.argmax(Spatial_Correspondence, axis=0)
            Miss_Match = np.concatenate((ps[:, np.newaxis] + 1, ps2[ps, np.newaxis] + 1), axis=1)

        # Functional homogeneity
        dir_scan = list_scan[i]
        scan_data = pNet.load_fmri_scan(dir_scan, dataType=Data_Type, dataFormat=Data_Format, Normalization=None).astype(np_float)
        pFN_signal = scan_data @ pFN / np.sum(pFN, axis=0, keepdims=True)
        Corr_FH = pNet.mat_corr(pFN_signal, scan_data, dataPrecision=dataPrecision)
        Functional_Homogeneity = np.sum(Corr_FH.T * pFN,axis=0) / np.sum(pFN, axis=0)
        # Use gFN as control
        gFN_signal = scan_data @ gFN / np.sum(pFN, axis=0, keepdims=True)
        Corr_FH = pNet.mat_corr(gFN_signal, scan_data, dataPrecision=dataPrecision)
        Functional_Homogeneity_Control = np.sum(Corr_FH.T * gFN, axis=0) / np.sum(gFN, axis=0)

        # Finalize results
        Result = {'Spatial_Correspondence': Spatial_Correspondence,
                  'Delta_Spatial_Correspondence': Delta_Spatial_Correspondence,
                  'Miss_Match': Miss_Match,
                  'Functional_Homogeneity': Functional_Homogeneity,
                  'Functional_Homogeneity_Control': Functional_Homogeneity_Control}
        if Miss_Match.shape[0] > 0:
            flag_QC = 1
            print(' ' + str(Miss_Match.shape[0]) + ' miss matched FNs in sub folder: ' + list_subject_folder[i], file=file_Final_Report, flush=True)

        # Save results
        dir_pFN_indv_QC = os.path.join(dir_pnet_QC, list_subject_folder[i])
        if not os.path.exists(dir_pFN_indv_QC):
            os.makedirs(dir_pFN_indv_QC)
        scipy.io.savemat(os.path.join(dir_pFN_indv_QC, 'Result.mat'), {'Result': Result})

else:
    N_pFN = subject_ID_unique.shape[0]

# Finish the final report
if flag_QC == 0:
    print(f'\n All scans passed QC\n', file=file_Final_Report)
print('\nFinished QC at '+time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))+'\n', file=file_Final_Report, flush=True)
file_Final_Report.close()








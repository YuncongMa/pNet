# Yuncong Ma, 9/8/2023
# Test the bootstrapping function

import os
import re
import numpy as np
import pNet


# Prepare dataset: Scan_List, Subject_ID, Subject_Folder
dir_raw_data = '/Users/yuncongma/Documents/Document/fMRI/Myworks/pNet/Test/HCP_XCPD'
file_scan = os.path.join(dir_raw_data, 'Scan_List.txt')
file_subject_ID = os.path.join(dir_raw_data, 'Subject_ID.txt')
file_subject_folder = os.path.join(dir_raw_data, 'Subject_Folder.txt')
FID = open(file_scan, 'w')
N_Scan = 0
for root, dirs, files in os.walk(dir_raw_data):
    for file in files:
        if file.endswith('dtseries.nii'):
            print(os.path.join(root, file), file=FID)
            N_Scan += 1
FID.close()
# Do not use readline to avoid the maximum txt file reading limit
list_scan = [line.replace('\n', '') for line in open(file_scan, 'r')]

FID = open(file_subject_ID, 'w')
for i in range(N_Scan):
    temp = re.findall(r'\d+', list_scan[i])
    subject_ID = temp[-3]
    print(subject_ID, file=FID)
FID.close()


FID = open(file_subject_folder, 'w')
for i in range(N_Scan):
    temp = re.findall(r'\d+', list_scan[i])
    subject_folder = os.path.join(temp[-3], 'REST' + temp[-2])
    print(subject_folder, file=FID)
FID.close()

# Prepare work directory in pNet result folder
dir_pnet_result = '/Users/yuncongma/Documents/Document/fMRI/Myworks/pNet/Test/Test_FN17_HCP_XCPD'
if not os.path.exists(dir_pnet_result):
    os.mkdir(dir_pnet_result)
dir_pnet_BS = os.path.join(dir_pnet_result, 'FN_Computation', 'BootStrapping')
if not os.path.exists(dir_pnet_BS):
    os.makedirs(dir_pnet_BS)

dir_output = '/Users/yuncongma/Documents/Document/fMRI/Myworks/pNet/Test/Test_FN17_HCP_XCPD/FN_Computation/BootStrapping'
pNet.bootstrap_scan(dir_output, file_scan, file_subject_ID, file_subject_folder, file_group=None, BS=2, N_BS=10, combineFLag=1, samplingMethod='Subject')
raise ValueError('Stop')
# Lists for input
list_scan = np.array([line.replace('\n', '') for line in open(file_scan, 'r')])
list_subject_ID = np.array([line.replace('\n', '') for line in open(file_subject_ID, 'r')])
subject_ID_unique = np.unique(list_subject_ID)
N_Subject = subject_ID_unique.shape[0]
list_subject_folder = np.array([line.replace('\n', '') for line in open(file_subject_folder, 'r')])
list_group = ''
if list_group is not None and len(list_group) > 0:
    group_unique = np.unique(list_subject_ID)

# Create bootstrapped scan lists
N_BS = 10
BS = 2
samplingMethod = 'Subject'  # 'Subject', or 'Group_Subject'
combineFLag = 1  # Whether to combine multiple fMRI scans for the same subject

# check parameter
if BS >= N_Subject:
    raise ValueError('The number of randomly selected subjects should be fewer than the total number of subjects')
if samplingMethod == 'Group_Subject' and (list_group is None or len(list_group) == 0):
    raise ValueError('Group information is absent')

for i in range(1, N_BS+1):
    if not os.path.exists(os.path.join(dir_pnet_BS, str(i))):
        os.mkdir(os.path.join(dir_pnet_BS, str(i)))
    List_BS = np.empty(BS, dtype=list)

    # Randomly select subjects
    if samplingMethod == 'Subject':
        ps = np.sort(np.random.choice(N_Subject, BS, replace=False))
        for j in range(BS):
            if combineFLag == 1:
                # Get all scans from the selected subject
                temp = list_scan[np.where(np.compare_chararrays(list_subject_ID, subject_ID_unique[ps[j]], '==', False))[0]]
                List_BS[j] = str.join('\n', temp)
            else:
                # Choose one scan from the selected subject
                temp = list_scan[np.where(np.compare_chararrays(list_subject_ID, subject_ID_unique[ps[j]], '==', False))[0]]
                ps2 = np.random.choice(temp.shape[0], 1)  # list
                List_BS[j] = temp[ps2[0]]  # transform to string list

    if samplingMethod == 'Group_Subject':
        break

    # Write the Scan_List.txt file
    FID = open(os.path.join(dir_pnet_BS, str(i), 'Scan_List.txt'), 'w')
    for j in range(BS):
        print(List_BS[j], file=FID)
    FID.close()







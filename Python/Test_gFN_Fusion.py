# Yuncong Ma, 9/1/2023
# Test the fusion of gFNs

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
import h5py
import scipy

import pNet

folder_gFN = os.path.join(pNet.Example.HCP_surf.dir_pnet, 'Group_FN', 'BootStrapping')
file_output = os.path.join(pNet.Example.HCP_surf.dir_pnet, 'Group_FN', 'FN.mat')
file_setting = os.path.join(pNet.Example.HCP_surf.dir_pnet, 'FN_Computation', 'Setting.json')

# Load setting
setting = pNet.load_json_setting(file_setting)

# Setup parameters
K = setting['K']
BS_Repetition = setting['GroupFN']['BootStrap']['Repetition']
NCut_MaxTrial = 100
dataPrecision = 'double'
logFile = 'Log_gFN_fusion.Log'

# setup log file
logFile = open(logFile, 'w+')
print(f'\nStart NCut for gFN fusion using PyTorch at '+time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))+'\n', file=logFile)


# Setup data precision and eps
np_float, np_eps = pNet.set_data_precision(dataPrecision)

# Setup file directories of gFN results
list_gFN = [os.path.join(folder_gFN, str(i), 'FN.mat') for i in range(1, BS_Repetition + 1)]

# Load gFNs from bootstrapping
# Concatenate FN along the FN dimension, gFN_BS is [dim_space, K * n_BS]
gFN_BS = [(pNet.load_matlab_single_array(list_gFN[i])).astype(np_float) for i in range(0, BS_Repetition)]
gFN_BS = np.concatenate(gFN_BS, axis=1)

########################################################################################################################
# clustering by NCut

# Get similarity between samples
corrVal = pNet.mat_corr(gFN_BS)  # similarity between FNs, [K * n_BS, K * n_BS]
corrVal[np.isnan(corrVal)] = -1
nDis = 1 - corrVal  # Transform Pearson correlation to non-negative values similar to distance
triuInd = np.triu(np.ones(nDis.shape), 1)  # Index of upper triangle
nDisVec = nDis[triuInd == 1]  # Take the values in upper triangle
# Make all distance non-negative and normalize their distribution
nW = np.exp(-np.power(nDis, 2) / (np.power(np.median(nDisVec), 2)))  # Transform distance values using exp(-X/std^2) with std as the median value
nW[np.isnan(nW)] = 0
sumW = np.sum(nW, axis=1)  # total distance for each FN
sumW[sumW == 0] = 1  # In case two FNs are the same
# Construct Laplacian matrix
D = np.diag(sumW)
L = np.sqrt(np.linalg.inv(D)) @ nW @ np.linalg.inv(np.sqrt(D))  # A way to normalize nW based on the total distance of each FN
L = (L + L.T) / 2  # Ensure L is symmetric. Computation error may result in asymmetry

eVal, Ev = scipy.sparse.linalg.eigs(L.astype(np.float64), K, which='LR')  # Get first K eigenvectors, sign of vectors may be different to MATLAB results
Ev = np.real(Ev)
# Correct the sign of eigenvectors to make them same as derived from MATLAB
temp = np.sign(np.sum(Ev, axis=0))  # Use the total value of each eigenvector to reset its sign
temp[temp == 0.0] = 1.0
Ev = Ev * np.tile(temp, (Ev.shape[0], 1))  # Reset the sign of each eigenvector
normvect = np.sqrt(np.diag(Ev @ Ev.T))  # Get the norm of each eigenvector
normvect[normvect == 0.0] = 1  # Incase all 0 eigenvector
Ev = np.linalg.solve(np.diag(normvect), Ev)  # Use linear solution to normalize Ev satisfying normvect * Ev_new = Ev_old

# Multiple trials to get reproducible results
Best_C = []
Best_NCutValue = np.inf
for i in range(NCut_MaxTrial):

    EigenVectors = Ev
    n, k = EigenVectors.shape  # n is K * n_BS, k is K

    vm = np.sqrt(np.sum(EigenVectors**2, axis=1, keepdims=True))  # norm of each row
    EigenVectors = EigenVectors / np.tile(vm, (1, k))  # normalize eigenvectors to ensure each FN vector's norm = 1

    R = np.zeros((k, k))
    ps = np.random.randint(0, n, 1)  # Choose a random row in eigenvectors
    R[:, 0] = EigenVectors[ps, :]  # This randomly selected row in eigenvectors is used as an initial center

    c_index = np.zeros(n)
    c = np.zeros(n)  # Total distance to different rows in R [K * n_BS, 1]
    c_index[0] = ps  # Store the index of selected samples
    c[ps] = np.inf

    for j in range(2, k+1):  # Find another K-1 rows in eigenvectors which have the minimum similarity to previous selected rows, similar to initialization in k++
        c += np.abs(EigenVectors @ R[:, j-2])
        ps = np.argmin(c)
        c_index[j-1] = ps
        c[ps] = np.inf
        R[:, j-1] = EigenVectors[ps, :].T

    lastObjectiveValue = 0
    exitLoop = 0
    nbIterationsDiscretisation = 0
    nbIterationsDiscretisationMax = 20

    while exitLoop == 0:
        nbIterationsDiscretisation += 1

        EigenVectorsR = EigenVectors @ R
        n, k = EigenVectorsR.shape
        J = np.argmax(EigenVectorsR, axis=1)  # Assign each sample to K centers of R based on highest similarity
        EigenvectorsDiscrete = scipy.sparse.csr_matrix((np.ones(n), (np.arange(n), J)), shape=(n, k))  # Generate a 0-1 matrix with each row containing only one 1

        U, S, Vh = scipy.linalg.svd(EigenvectorsDiscrete.T @ EigenVectors, full_matrices=False)  # Economy-size decomposition

        S = np.diag(S)  # To match MATLAB svd results
        V = Vh.T  # To match MATLAB svd results
        NcutValue = 2 * (n - np.trace(S))

        # escape the loop when converged or meet max iteration
        if abs(NcutValue - lastObjectiveValue) < np_eps or nbIterationsDiscretisation > nbIterationsDiscretisationMax:
            exitLoop = 1
            print(f'\nReach stop criterion of NCut, NcutValue = '+str(NcutValue)+'\n', file=logFile)
        else:
            print(f' NcutValue = '+str(NcutValue), file=logFile)
            lastObjectiveValue = NcutValue
            R = V @ U.T  # Update R which stores the new centers

    C = np.argmax(EigenvectorsDiscrete.toarray(), axis=1)  # Assign each sample to K centers in R

    if len(np.unique(C)) < K:  # Check whether there are empty results
        print(f'Found empty results in iteration '+str(i+1)+'\n', file=logFile)
    else:  # Update the best result
        if NcutValue < Best_NCutValue:
            Best_NCutValue = NcutValue
            Best_C = C

if len(set(Best_C)) < K:  # In case even the last trial has empty results
    raise ValueError('Cannot generate non-empty gFNs\n')
    Flag = 1
    Message = "Cannot generate non-empty FN"

print(f'Best NCut value = '+str(Best_NCutValue)+'\n', file=logFile)
########################################################################################################################

# Get centroid
C = Best_C
gFN = np.zeros((gFN_BS.shape[0], K))
for ki in range(K):
    if np.sum(C == ki) > 1:
        candSet = gFN_BS[:, C == ki]  # Get the candidate set of FNs assigned to cluster ki
        corrW = np.abs(pNet.mat_corr(candSet))  # Get the similarity between candidate FNs
        corrW[np.isnan(corrW)] = 0
        mInd = np.argmax(np.sum(corrW, axis=0), axis=0)  # Find the FN with the highest total similarity to all other FNs
        gFN[:, ki] = candSet[:, mInd]
    elif np.sum(C == ki) == 1:
        mInd = C.index(ki)
        gFN[:, ki] = gFN_BS[:, mInd]

gFN = gFN / np.maximum(np.tile(np.max(gFN, axis=0), (gFN.shape[0], 1)), np_eps)  # Normalize each FN by its max value

FN = gFN

print(f'\nFinished at '+time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))+'\n', file=logFile)

savemat('gFN_initV_FN.mat', {'FN': FN})









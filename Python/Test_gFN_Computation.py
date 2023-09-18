# Yuncong Ma, 9/13/2023
# Test the computation of gFNs

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


os.system('export OMP_NUM_THREADS=1')
os.system('export OPENBLAS_NUM_THREADS=1')

#########################################
# Example data from HCP surf Test_FN17
# gNb uses index starting from 1, but it starts from 0 in Python
gNb = pNet.load_matlab_single_array(os.path.join(pNet.Example.HCP_surf.dir_pnet, 'FN_Computation', 'gNb.mat'))
Data = pNet.load_matlab_single_array(os.path.join(pNet.Example.HCP_surf.dir_data, '100206', '1', 'LR', 'Image.mat'))
Data = Data.T
Setting = pNet.load_json_setting(os.path.join(pNet.Example.HCP_surf.dir_pnet, 'FN_Computation', 'Setting.json'))

# Setup parameters
K = 17
maxIter = 1000
minIter = 30
error = 1e-6
normW = 1
eta = 0
Alpha = 2
Beta = 30
vxI = 0
ard = 0
dataPrecision = 'double'
logFile = 'Log_gFN_NMF.log'

# Internal parameters
alphaS = 0
alphaL = 0

# setup log file
logFile = open(logFile, 'w')
print(f'\nStart NMF for gFN using NumPy at '+time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))+'\n', file=logFile)

# Input data size
dim_time, dim_space = Data.shape

# Median number of graph neighbors
nM = np.median(np.unique(gNb[:, 0], return_counts=True)[1])

# Use Alpha and Beta to set alphaS and alphaL if they are 0
if alphaS == 0 and Alpha > 0:
    alphaS = np.round(Alpha * dim_time / K)
if alphaL == 0 and Beta > 0:
    alphaL = np.round(Beta * dim_time / K / nM)

# Setup data precision and eps
np_float, np_eps = pNet.set_data_precision(dataPrecision)

Data = Data.astype(np_float)

# Prepare and normalize scan
X = pNet.normalize_data(Data, 'vp', 'vmax')

# Initialize U and V
mean_X = np.sum(X) / (dim_time*dim_space)
rng = np.random.default_rng(1)
U = (rng.random((dim_time, K), dtype=np_float) + 1) * (np.sqrt(mean_X/K))
V = (rng.random((dim_space, K), dtype=np_float) + 1) * (np.sqrt(mean_X/K))

# To ensure the same results
U = pNet.load_matlab_single_array('gFN_init_U.mat')
V = pNet.load_matlab_single_array('gFN_init_V.mat')

# Normalize data
U, V = pNet.normalize_u_v(U, V, 1, 1)

# Construct the spatial affinity graph
L, W, D = pNet.construct_Laplacian_gNb(gNb, dim_space, vxI, X, alphaL, normW, dataPrecision)

if ard > 0:
    ard = 1
    eta = 0.1
    lambdas = np.sum(U, axis=0) / dim_time
    hyperLam = eta * np.sum(np.power(X, 2)) / (dim_time * dim_space * 2)
else:
    lambdas = 0
    hyperLam = 0

oldLogL = np.inf

# Multiplicative update of U and V
for i in range(1, 1+maxIter):
    # ===================== update V ========================
    # Eq. 8-11
    XU = X.T @ U
    UU = U.T @ U
    VUU = V @ UU

    tmpl2 = np.power(V, 2)

    if alphaS > 0:
        tmpl21 = np.sqrt(tmpl2)
        tmpl22 = np.tile(np.sqrt(np.sum(tmpl2, axis=0)), (dim_space, 1))
        tmpl21s = np.tile(np.sum(tmpl21, axis=0), (dim_space, 1))
        posTerm = V / np.maximum(tmpl21 * tmpl22, np_eps)
        negTerm = V * tmpl21s / np.maximum(np.power(tmpl22, 3), np_eps)

        VUU = VUU + 0.5 * alphaS * posTerm
        XU = XU + 0.5 * alphaS * negTerm

    if alphaL > 0:
        WV = W @ V.astype(np.float64)
        DV = D @ V.astype(np.float64)

        XU = XU + WV
        VUU = VUU + DV

    V = V * (XU / np.maximum(VUU, np_eps))

    # Prune V if empty components are found in V
    # This is almost impossible to happen without combining FNs
    prunInd = np.sum(V != 0, axis=0) == 1
    if np.any(prunInd):
        V[:, prunInd] = np.zeros((dim_space, np.sum(prunInd)))
        U[:, prunInd] = np.zeros((dim_time, np.sum(prunInd)))

    # normalize U and V
    U, V = pNet.normalize_u_v(U, V, 1, 1)

    # ===================== update U =========================
    XV = X @ V
    VV = V.T @ V
    UVV = U @ VV

    if ard > 0:  # ard term for U
        posTerm = 1 / np.maximum(np.tile(lambdas, (dim_time, 1)), np_eps)

        UVV = UVV + posTerm * hyperLam

    U = U * (XV / np.maximum(UVV, np_eps))

    # Prune U if empty components are found in U
    # This is almost impossible to happen without combining FNs
    prunInd = np.sum(U, axis=0) == 0
    if np.any(prunInd):
        V[:, prunInd] = np.zeros((dim_space, np.sum(prunInd)))
        U[:, prunInd] = np.zeros((dim_time, np.sum(prunInd)))

    # update lambda
    if ard > 0:
        lambdas = np.sum(U) / dim_time

    # ==== calculate objective function value ====
    ardU = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmpl21 = np.power(V, 2)

    if ard > 0:
        su = np.sum(U, axis=0)
        su[su == 0] = 1
        ardU = np.sum(np.log(su)) * dim_time * hyperLam

    if alphaL > 0:
        tmp2 = V.T @ L * V.T

    L21 = alphaS * np.sum(np.sum(np.sqrt(tmpl21), axis=0) / np.maximum(np.sqrt(np.sum(tmpl21, axis=0)), np_eps))
    # LDf = pNet.data_fitting_error(X, U, V, 0, 1)
    LDf = np.sum(np.power(X - U @ V.T, 2))
    LSl = np.sum(tmp2)

    # Objective function
    LogL = L21 + LDf + LSl + ardU
    print(f"    Iter = {i}: LogL: {LogL}, dataFit: {LDf}, spaLap: {LSl}, L21: {L21}, ardU: {ardU}", file=logFile)

    # The iteration needs to meet minimum iteration number and small changes of LogL
    if i > minIter and abs(oldLogL - LogL) / np.maximum(oldLogL, np_eps) < error:
        break
    oldLogL = LogL.copy()

savemat('gFN_U_end.mat', {'U': U})
savemat('gFN_V_end.mat', {'V': V})

print(f'\n Finished at '+time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))+'\n', file=logFile)

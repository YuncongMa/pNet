# Yuncong Ma, 8/29/2023
# Test the computation of gFNs with PyTorch

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
import torch

import pNet

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
logFile = 'Log_gFN_NMF_torch.log'

# Internal parameters
alphaS = 0
alphaL = 0

torch.set_num_threads(20)

# setup log file
logFile = open(logFile, 'w')
print(f'\nStart NMF for gFN using PyTorch at '+time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))+'\n', file=logFile)

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
torch_float, torch_eps = pNet.set_data_precision_torch(dataPrecision)

Data = torch.tensor(Data, dtype=torch_float)
error = torch.tensor(error, dtype=torch_float)
eta = torch.tensor(eta, dtype=torch_float)
alphaS = torch.tensor(alphaS, dtype=torch_float)

# Prepare and normalize scan
X = pNet.normalize_data_torch(Data, 'vp', 'vmax')

# Initialize U and V
mean_X = torch.divide(torch.sum(X), torch.tensor(dim_time*dim_space))
U = (torch.rand((dim_time, K), dtype=torch_float) + 1) * (torch.sqrt(torch.div(mean_X, K)))
V = (torch.rand((dim_space, K), dtype=torch_float) + 1) * (torch.sqrt(torch.div(mean_X, K)))

# To ensure the same results
# U = pNet.load_matlab_single_array('gFN_init_U.mat')
# V = pNet.load_matlab_single_array('gFN_init_V.mat')

# Disable torch gradient
torch.set_grad_enabled(False)

# Normalize data
U, V = pNet.normalize_u_v_torch(U, V, 1, 1)

# Construct the spatial affinity graph
L, W, D = pNet.construct_Laplacian_gNb_torch(gNb, dim_space, vxI, X, alphaL, normW, dataPrecision)


if ard > 0:
    ard = 1
    eta = 0.1
    lambdas = torch.sum(U, dim=0) / dim_time
    hyperLam = eta * torch.sum(torch.pow(X, 2)) / (dim_time * dim_space * 2)
else:
    lambdas = 0
    hyperLam = 0

oldLogL = torch.inf

# Multiplicative update of U and V
for i in range(1, 1+maxIter):
    # ===================== update V ========================
    # Eq. 8-11
    XU = torch.matmul(X.T, U)
    UU = torch.matmul(U.T, U)
    VUU = torch.matmul(V, UU)

    tmpl2 = torch.pow(V, 2)

    if alphaS > 0:
        tmpl21 = torch.sqrt(tmpl2)
        tmpl22 = torch.tile(torch.sqrt(torch.sum(tmpl2, dim=0)), (dim_space, 1))
        tmpl21s = torch.tile(torch.sum(tmpl21, dim=0), (dim_space, 1))
        posTerm = torch.div(V, torch.maximum(torch.mul(tmpl21, tmpl22), torch_eps))
        negTerm = torch.div(torch.mul(V, tmpl21s), torch.maximum(torch.pow(tmpl22, 3), torch_eps))

        VUU = VUU + 0.5 * alphaS * posTerm
        XU = XU + 0.5 * alphaS * negTerm

    if alphaL > 0:
        WV = torch.matmul(W, V.type(torch.float64))
        DV = torch.matmul(D, V.type(torch.float64))

        XU = torch.add(XU, WV)
        VUU = torch.add(VUU, DV)

    V = torch.mul(V, (torch.div(XU, torch.maximum(VUU, torch_eps))))

    # Prune V if empty components are found in V
    # This is almost impossible to happen without combining FNs
    prunInd = torch.sum(V != 0, dim=0) == 1
    if torch.any(prunInd):
        V[:, prunInd] = torch.zeros((dim_space, torch.sum(prunInd)))
        U[:, prunInd] = torch.zeros((dim_time, torch.sum(prunInd)))

    # normalize U and V
    U, V = pNet.normalize_u_v_torch(U, V, 1, 1)

    # ===================== update U =========================
    XV = torch.matmul(X, V)
    VV = torch.matmul(V.T, V)
    UVV = torch.matmul(U, VV)

    if ard > 0:  # ard term for U
        posTerm = torch.div(torch.tensor(1.0), torch.maximum(torch.tile(lambdas, (dim_time, 1)), torch_eps))

        UVV = torch.add(UVV, posTerm * hyperLam)

    U = torch.mul(U, torch.div(XV, torch.maximum(UVV, torch_eps)))

    # Prune U if empty components are found in U
    # This is almost impossible to happen without combining FNs
    prunInd = torch.sum(U, dim=0) == 0
    if torch.any(prunInd):
        V[:, prunInd] = torch.zeros((dim_space, torch.sum(prunInd)))
        U[:, prunInd] = torch.zeros((dim_time, torch.sum(prunInd)))

    # update lambda
    if ard > 0:
        lambdas = torch.sum(U) / dim_time

    # ==== calculate objective function value ====
    ardU = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmpl21 = torch.pow(V, 2)

    if ard > 0:
        su = torch.sum(U, dim=0)
        su[su == 0] = 1
        ardU = torch.sum(torch.log(su)) * dim_time * hyperLam

    if alphaL > 0:
        tmp2 = torch.mul(torch.matmul(V.T, L), V.T)

    L21 = torch.mul(alphaS, torch.sum(torch.div(torch.sum(torch.sqrt(tmpl21), dim=0), torch.maximum(torch.sqrt(torch.sum(tmpl21, dim=0)), torch_eps))))
    # LDf = pNet.data_fitting_error(X, U, V, 0, 1)
    LDf = torch.sum(torch.pow(torch.subtract(X, torch.matmul(U, V.T)), 2))
    LSl = torch.sum(tmp2)

    # Objective function
    LogL = L21 + LDf + LSl + ardU
    print(f"    Iter = {i}: LogL: {LogL}, dataFit: {LDf}, spaLap: {LSl}, L21: {L21}, ardU: {ardU}", file=logFile)

    # The iteration needs to meet minimum iteration number and small changes of LogL
    if i > minIter and abs(oldLogL - LogL) / torch.maximum(oldLogL, torch_eps) < error:
        break
    oldLogL = LogL.clone()

savemat('gFN_U_end_torch.mat', {'U': U.numpy()})
savemat('gFN_V_end_torch.mat', {'V': V.numpy()})

print(f'\n Finished at '+time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))+'\n', file=logFile)

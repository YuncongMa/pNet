# Yuncong Ma, 9/1/2023
# Test the computation of pFNs using precomputed group-level FNs
# Use Pytorch

#########################################
# Packages
import numpy as np
import os
import re
from scipy.io import savemat
from scipy.sparse import spdiags, lil_matrix
import time
import torch
import pNet

#########################################
# Example data from HCP surf Test_FN17
gFN = pNet.load_matlab_array(os.path.join(pNet.Example.HCP_surf.dir_pnet, 'Group_FN','FN.mat'), 'FN')
# gNb uses index starting from 1, but it starts from 0 in Python
gNb = pNet.load_matlab_single_array(os.path.join(pNet.Example.HCP_surf.dir_pnet, 'FN_Computation', 'gNb.mat'))
scan = pNet.load_matlab_single_array(os.path.join(pNet.Example.HCP_surf.dir_data, '100206', '1', 'LR', 'Image.mat'))
Setting = pNet.load_setting(os.path.join(pNet.Example.HCP_surf.dir_pnet, 'FN_Computation', 'Setting.json'))

#########################################
# Setup parameters
K = gFN.shape[1]
maxIter = 1000
minIter = 30
meanFitRatio = 0.1
error = 1e-4
normW = 1
eta = Setting['PersonalizedFN']['eta']
Alpha = 2
Beta = 30
alphaS = Setting['PersonalizedFN']['alphaS21']
alphaL = Setting['PersonalizedFN']['alphaL']
vxI = 0
initConv = 1
ard = Setting['PersonalizedFN']['ard']
dataPrecision = 'double'
logFile = 'Log_pFN_NMF_torch.log'

# setup log file
logFile = open(logFile, 'w')
print(f'\nStart NMF for pFN using PyTorch at '+time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))+'\n', file=logFile)

#########################################
# Setup data precision and eps
torch_float, torch_eps = pNet.set_data_precision_torch(dataPrecision)

scan = torch.tensor(scan, dtype=torch_float)
gFN = torch.tensor(gFN, dtype=torch_float)
error = torch.tensor(error, dtype=torch_float)
eta = torch.tensor(eta, dtype=torch_float)
alphaS = torch.tensor(alphaS, dtype=torch_float)


# initialization
initV = gFN

dim_space, dim_time = scan.shape

# Median number of graph neighbors
nM = np.median(np.unique(gNb[:, 0], return_counts=True)[1])

# Use Alpha and Beta to set alphaS and alphaL if they are 0
if alphaS == 0 and Alpha > 0:
    alphaS = torch.tensor(np.round(Alpha * dim_time / K))
if alphaL == 0 and Beta > 0:
    alphaL = np.round(Beta * dim_time / K / nM)

#########################################
# Prepare and normalize scan
X = pNet.normalize_data_torch(scan.T, 'vp', 'vmax')


# Construct the spatial affinity graph
L, W, D = pNet.construct_Laplacian_gNb_torch(gNb, dim_space, vxI, X, alphaL, normW, dataPrecision)


# Initialize V
V = initV.clone()
miv = torch.max(V, dim=0)[0]
trimInd = V / torch.maximum(torch.tile(miv, (dim_space, 1)), torch_eps) < torch.tensor(5e-2)
V[trimInd] = 0

# Initialize U
U = X @ V / torch.tile(torch.sum(V, dim=0), (dim_time, 1))

U = pNet.initialize_u_torch(X, U, V, error, maxIter, minIter, meanFitRatio, initConv)

initU = U.clone()
initV = V.clone()

# Alternative update of U and V
# Variables
U = initU.clone()
V = initV.clone()

if ard > 0:
    lambdas = torch.sum(U, dim=0) / dim_time
    hyperLam = eta * torch.sum(torch.pow(X, 2)) / (dim_time * dim_space * 2)
else:
    lambdas = 0
    hyperLam = 0

flagQC = 0
oldLogL = torch.inf
oldU = U.clone()
oldV = V.clone()

for i in range(1, 1+maxIter):
    # ===================== update V ========================
    # Eq. 8-11
    XU = X.T @ U
    UU = U.T @ U
    VUU = V @ UU

    tmpl2 = torch.pow(V, 2)

    if alphaS > 0:
        tmpl21 = torch.sqrt(tmpl2)
        tmpl22 = torch.tile(torch.sqrt(torch.sum(tmpl2, dim=0)), (dim_space, 1))
        tmpl21s = torch.tile(torch.sum(tmpl21, dim=0), (dim_space, 1))
        posTerm = V / torch.maximum(tmpl21 * tmpl22, torch_eps)
        negTerm = V * tmpl21s / torch.maximum(torch.pow(tmpl22, 3), torch_eps)

        VUU = VUU + 0.5 * alphaS * posTerm
        XU = XU + 0.5 * alphaS * negTerm

    if alphaL > 0:
        WV = W @ V
        DV = D @ V

        XU = XU + WV
        VUU = VUU + DV

    V = V * (XU / torch.maximum(VUU, torch_eps))

    # Prune V if empty components are found in V
    # This is almost impossible to happen without combining FNs
    prunInd = torch.sum(V != 0, dim=0) == 1
    if torch.any(prunInd):
        V[:, prunInd] = torch.zeros((dim_space, torch.sum(prunInd)))
        U[:, prunInd] = torch.zeros((dim_time, torch.sum(prunInd)))

    # normalize U and V
    U, V = pNet.normalize_u_v_torch(U, V, 1, 1)

    # ===================== update U =========================
    XV = X @ V
    VV = V.T @ V
    UVV = U @ VV

    if ard > 0:  # ard term for U
        posTerm = torch.tensor(1) / torch.maximum(torch.tile(lambdas, (dim_time, 1)), torch_eps)

        UVV = UVV + posTerm * hyperLam

    U = U * (XV / torch.maximum(UVV, torch_eps))

    # Prune U if empty components are found in U
    # This is almost impossible to happen without combining FNs
    prunInd = torch.sum(U, dim=0) == 0
    if torch.any(prunInd):
        V[:, prunInd] = torch.zeros((dim_space, torch.sum(prunInd)))
        U[:, prunInd] = torch.zeros((dim_time, torch.sum(prunInd)))

    # update lambda
    if ard > 0:
        lambdas = torch.sum(U, dim=0) / dim_time

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
        tmp2 = V.T @ L * V.T

    L21 = alphaS * torch.sum(torch.sum(torch.sqrt(tmpl21), dim=0) / torch.maximum(torch.sqrt(torch.sum(tmpl21, dim=0)), torch_eps))
    LDf = pNet.data_fitting_error(X, U, V, 0, 1)
    LSl = torch.sum(tmp2)

    # Objective function
    LogL = L21 + ardU + LDf + LSl
    print(f"    Iter = {i}: LogL: {LogL}, dataFit: {LDf}, spaLap: {LSl}, L21: {L21}, ardU: {ardU}", file=logFile)

    # The iteration needs to meet minimum iteration number and small changes of LogL
    if i > minIter and abs(oldLogL - LogL) / torch.maximum(oldLogL, torch_eps) < error:
        break
    oldLogL = LogL.clone()

    # QC Control
    temp = pNet.mat_corr_torch(initV, V)
    QC_Spatial_Correspondence = torch.clone(torch.diag(temp))
    temp -= torch.diag(torch.diag(temp))
    QC_Spatial_Correspondence_Control = torch.max(temp, dim=1)[0]
    QC_Delta_Sim = torch.min(QC_Spatial_Correspondence - QC_Spatial_Correspondence_Control)
    QC_Delta_Sim = QC_Delta_Sim.cpu().numpy()

    if QC_Delta_Sim <= 0:
        flagQC = 1
        U = oldU.clone()
        V = oldV.clone()
        print(f'\n  QC: Meet QC constraint: Delta sim = {QC_Delta_Sim}', file=logFile)
        print(f'    Use results from last iteration', file=logFile)
        break
    else:
        oldU = U.clone()
        oldV = V.clone()
        print(f'        QC: Delta sim = {QC_Delta_Sim}', file=logFile)

savemat('pFN_torch_U_end.mat', {'U': U.cpu().numpy()})
savemat('pFN_torch_V_end.mat', {'V': V.cpu().numpy()})

print(f'\n Finished at '+time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))+'\n', file=logFile)

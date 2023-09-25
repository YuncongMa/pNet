# Yuncong Ma, 9/5/2023
# FN Computation module of pNet
# Pytorch version

#########################################
# Packages
import scipy.io as sio
import os
import re
import time
import torch


# other functions of pNet
from Data_Input import *
from FN_Computation import construct_Laplacian_gNb, compute_gNb, bootstrap_scan, setup_NMF_setting, setup_pFN_folder


def mat_corr_torch(X, Y=None, dataPrecision='double'):
    """
    mat_corr_torch(X, Y=None, dataPrecision='double')
    Perform corr as in MATLAB, pair-wise Pearson correlation between columns in X and Y

    :param X: 1D or 2D matrix, numpy.ndarray or torch.Tensor
    :param Y: 1D or 2D matrix, or None, numpy.ndarray or torch.Tensor
    :param dataPrecision: 'double' or 'single'
    X and Y have the same number of rows
    :return: Corr

    Note: this method will use memory as it concatenates X and Y along column direction.
    By Yuncong Ma, 9/5/2023
    """

    torch_float, torch_eps = set_data_precision_torch(dataPrecision)
    if not isinstance(X, torch.Tensor):
        X = torch.tensor(X, dtype=torch_float)
    else:
        X = X.type(torch_float)
    if Y is not None:
        if not isinstance(Y, torch.Tensor):
            Y = torch.tensor(Y, dtype=torch_float)
        else:
            Y = Y.type(torch_float)

    # Check size of X and Y
    if len(X.shape) > 2 or (Y is not None and len(Y.shape) > 2):
        raise ValueError("X and Y must be 1D or 2D matrices")
    if Y is not None and X.shape[0] != Y.shape[0]:
        raise ValueError("X and Y must have the same number of columns")

    if Y is not None:
        # Subtract the mean to calculate the covariance
        X_centered = X - torch.mean(X, dim=0, keepdim=True)
        Y_centered = Y - torch.mean(Y, dim=0, keepdim=True)
        # Compute the standard deviation of the columns
        std_X = torch.std(X_centered, dim=0, keepdim=True, unbiased=True)
        std_Y = torch.std(Y_centered, dim=0, keepdim=True, unbiased=True)
        # Compute the correlation matrix
        numerator = (X_centered.T @ Y_centered)
        denominator = (std_X.T @ std_Y) * torch.tensor(X.shape[0] - 1)
        Corr = numerator / denominator
    else:
        if len(X.shape) != 2:
            raise ValueError("X must be a 2D matrix")
        X_centered = X - torch.mean(X, dim=0, keepdim=True)
        std_X = torch.std(X_centered, dim=0, keepdim=True, unbiased=True)
        numerator = (X_centered.T @ X_centered)
        denominator = (std_X.T @ std_X) * torch.tensor(X.shape[0] - 1)
        Corr = numerator / denominator

    return Corr


def normalize_data_torch(data, algorithm='vp', normalization='vmax', dataPrecision='double'):
    """
    normalize_data_torch(data, algorithm='vp', normalization='vmax', dataPrecision='double')
    Normalize data by algorithm and normalization settings

    :param data: data in 2D matrix [dim_time, dim_space], numpy.ndarray or torch.Tensor, recommend to use reference mode to save memory
    :param algorithm: 'z' 'gp' 'vp'
    :param normalization: 'n2' 'n1' 'rn1' 'g' 'vmax'
    :param dataPrecision: 'double' or 'single'
    :return: data

    Consistent to MATLAB function normalize_data(X, algorithm, normalization, dataPrecision)
    By Yuncong Ma, 9/8/2023
    """

    if len(data.shape) != 2:
        raise ValueError("data must be a 2D matrix")

    torch_float, torch_eps = set_data_precision_torch(dataPrecision)
    if not isinstance(data, torch.Tensor):
        data = torch.tensor(data, dtype=torch_float)
    else:
        data = data.type(torch_float)

    if algorithm.lower() == 'z':
        # standard score for each variable
        mVec = torch.mean(data, dim=1)
        sVec = torch.maximum(torch.std(data, dim=1), torch_eps)
        data = (data - mVec[:, torch.newaxis]) / sVec[:, torch.newaxis]
    elif algorithm.lower() == 'gp':
        # remove negative value globally
        minVal = torch.min(data)
        shiftVal = torch.abs(torch.minimum(minVal, torch.tensor(0.0)))
        data += shiftVal
    elif algorithm.lower() == 'vp':
        # remove negative value voxel-wisely
        minVal = torch.min(data, dim=0, keepdim=True)[0]
        shiftVal = torch.abs(torch.minimum(minVal, torch.tensor(0.0)))
        data += shiftVal
    else:
        # do nothing
        data = data

    if normalization.lower() == 'n2':
        # l2 normalization for each observation
        l2norm = torch.sqrt(torch.sum(data ** 2, dim=1)) + torch_eps
        data = data / l2norm[:, torch.newaxis]
    elif normalization.lower() == 'n1':
        # l1 normalization for each observation
        l1norm = torch.sum(data, dim=1) + torch_eps
        data = data / l1norm[:, torch.newaxis]
    elif normalization.lower() == 'rn1':
        # l1 normalization for each variable
        l1norm = torch.sum(data, dim=0) + torch_eps
        data = data / l1norm
    elif normalization.lower() == 'g':
        # global scale
        sVal = torch.sort(data, dim=None)
        perT = 0.001
        minVal = sVal[int(len(sVal) * perT)]
        maxVal = sVal[int(len(sVal) * (1 - perT))]
        data[data < minVal] = minVal
        data[data > maxVal] = maxVal
        data = (data - minVal) / max((maxVal - minVal), torch_eps)
    elif normalization.lower() == 'vmax':
        cmin = torch.min(data, dim=0, keepdim=True).values
        cmax = torch.max(data, dim=0, keepdim=True).values
        data = (data - cmin) / torch.maximum(cmax - cmin, torch_eps)
    else:
        # do nothing
        data = data

    if torch.isnan(data).any():
        raise ValueError('  nan exists, check the preprocessed data')

    return data


def initialize_u_torch(X, U0, V0, error=1e-4, maxIter=1000, minIter=30, meanFitRatio=0.1, initConv=1, dataPrecision='double'):
    """
    initialize_u_torch(X, U0, V0, error=1e-4, maxIter=1000, minIter=30, meanFitRatio=0.1, initConv=1, dataPrecision='double')

    :param X: data, 2D matrix [dim_time, dim_space], numpy.ndarray or torch.Tensor
    :param U0: initial temporal component, 2D matrix [dim_time, k], numpy.ndarray or torch.Tensor
    :param V0: initial spatial component, 2D matrix [dim_space, k], numpy.ndarray or torch.Tensor
    :param error: data fitting error
    :param maxIter: maximum iteration
    :param minIter: minimum iteration
    :param meanFitRatio: a 0-1 scaler, exponential moving average coefficient
    :param initConv: 0 or 1, flag for convergence
    :param dataPrecision: 'double' or 'single'
    :return: U_final: temporal components of FNs, a 2D matrix [dim_time, K]

    Consistent to MATLAB function initialize_u(X, U0, V0, error, maxIter, minIter, meanFitRatio, initConv, dataPrecision)
    By Yuncong Ma, 9/5/2023
    """

    torch_float, torch_eps = set_data_precision_torch(dataPrecision)
    if not isinstance(X, torch.Tensor):
        X = torch.tensor(X, dtype=torch_float)
    else:
        X = X.type(torch_float)
    if not isinstance(U0, torch.Tensor):
        U0 = torch.tensor(U0, dtype=torch_float)
    else:
        U0 = U0.type(torch_float)
    if not isinstance(X, torch.Tensor):
        V0 = torch.tensor(V0, dtype=torch_float)
    else:
        V0 = V0.type(torch_float)

    # Check the data size of X, U0 and V0
    if len(X.shape) != 2 or len(U0.shape) != 2 or len(V0.shape) != 2:
        raise ValueError("X, U0 and V0 must be 2D matrices")
    if X.shape[0] != U0.shape[0] or X.shape[1] != V0.shape[0] or U0.shape[1] != V0.shape[1]:
        raise ValueError("X, U0 and V0 need to have appropriate sizes")

    U = U0.clone()
    V = V0.clone()

    newFit = data_fitting_error_torch(X, U, V, 0, 1, dataPrecision)
    meanFit = newFit / meanFitRatio

    maxErr = 1
    for i in range(1, maxIter+1):
        # update U with fixed V
        XV = X @ V
        VV = V.T @ V
        UVV = U @ VV

        U = U * (XV / torch.maximum(UVV, torch_eps))

        if i > minIter:
            if initConv:
                newFit = data_fitting_error_torch(X, U, V, 0, 1, dataPrecision)
                meanFit = meanFitRatio * meanFit + (1 - meanFitRatio) * newFit
                maxErr = (meanFit - newFit) / meanFit

        if maxErr <= error:
            break

    U_final = U
    return U_final


def data_fitting_error_torch(X, U, V, deltaVU=0, dVordU=1, dataPrecision='double'):
    """
    data_fitting_error(X, U, V, deltaVU, dVordU, dataPrecision='double')
    Calculate the datat fitting of X'=UV' with terms

    :param X: 2D matrix, [Space, Time]
    :param U: 2D matrix, [Time, k]
    :param V: 2D matrix, [Space, k]
    :param deltaVU: 0
    :param dVordU: 1
    :param dataPrecision: 'double' or 'single'
    :return: Fitting_Error

    Consistent to MATLAB function fitting_initialize_u(X, U, V, deltaVU, dVordU, dataPrecision)
    By Yuncong Ma, 9/6/2023
    """

    torch_float, torch_eps = set_data_precision_torch(dataPrecision)
    if not isinstance(X, torch.Tensor):
        X = torch.tensor(X, dtype=torch_float)
    else:
        X = X.type(torch_float)
    if not isinstance(U, torch.Tensor):
        U = torch.tensor(U, dtype=torch_float)
    else:
        U = U.type(torch_float)
    if not isinstance(V, torch.Tensor):
        V = torch.tensor(V, dtype=torch_float)
    else:
        V = V.type(torch_float)

    # Check data size of X, U and V
    if len(X.shape) != 2 or len(U.shape) != 2 or len(V.shape) != 2:
        raise ValueError("X, U and V must be 2D matrices")
    if X.shape[0] != U.shape[0] or X.shape[1] != V.shape[0] or U.shape[1] != V.shape[1]:
        raise ValueError("X, U and V need to have appropriate sizes")

    dV = []
    maxM = 62500000  # To save memory
    dim_time, dim_space = X.shape
    mn = np.prod(X.shape)
    nBlock = int(np.floor(mn*3/maxM))
    if mn < maxM:
        dX = U @ V.T - X
        obj_NMF = torch.sum(torch.pow(dX, 2))
        if deltaVU:
            if dVordU:
                dV = dX.T * U
            else:
                dV = dX * V
    else:
        obj_NMF = 0
        if deltaVU:
            if dVordU:
                dV = torch.zeros_like(V)
            else:
                dV = torch.zeros_like(U)
        for i in range(int(np.ceil(dim_space/nBlock))):
            if i == int(np.ceil(dim_space/nBlock)):
                smpIdx = range(i*nBlock, dim_space)
            else:
                smpIdx = range(i*nBlock, np.minimum(dim_space, (i+1)*nBlock))
            dX = (U @ V[smpIdx, :].T) - X[:, smpIdx]
            obj_NMF += torch.sum(torch.sum(torch.pow(dX, 2)))
            if deltaVU:
                if dVordU:
                    dV[smpIdx, :] = torch.dot(dX.T, U)
                else:
                    dV += torch.dot(dX, V[smpIdx, :])
        if deltaVU:
            if dVordU:
                dV = dV

    Fitting_Error = obj_NMF
    return Fitting_Error


def normalize_u_v_torch(U, V, NormV, Norm, dataPrecision='double'):
    """
    normalize_u_v_torch(U, V, NormV, Norm, dataPrecision='double')
    Normalize U and V with terms

    :param U: 2D matrix, [Time, k]
    :param V: 2D matrix, [Space, k]
    :param NormV: 1 or 0
    :param Norm: 1 or 2
    :param dataPrecision: 'double' or 'single'
    :return: U, V

    Consistent to MATLAB function normalize_u_v(U, V, NormV, Norm, dataPrecision)
    By Yuncong Ma, 9/5/2023
    """

    torch_float, torch_eps = set_data_precision_torch(dataPrecision)
    if not isinstance(U, torch.Tensor):
        U = torch.tensor(U, dtype=torch_float)
    else:
        U = U.type(torch_float)
    if not isinstance(V, torch.Tensor):
        V = torch.tensor(V, dtype=torch_float)
    else:
        V = V.type(torch_float)

    # Check data size of U and V
    if len(U.shape) != 2 or len(V.shape) != 2:
        raise ValueError("U and V must be 2D matrices")
    if U.shape[1] != V.shape[1]:
        raise ValueError("U and V need to have appropriate sizes")

    dim_space = V.shape[0]
    dim_time = U.shape[0]

    if Norm == 2:
        norms = torch.sqrt(torch.sum(torch.pow(V, 2), dim=0))
        norms = torch.maximum(norms, torch_eps)
    else:
        norms = torch.max(V, dim=0)[0]  # torch.max return Value and Index
        norms = torch.maximum(norms, torch_eps)

    if NormV:
        U = U * torch.tile(norms, (dim_time, 1))
        V = V / torch.tile(norms, (dim_space, 1))
    else:
        U = U / torch.tile(norms, (dim_time, 1))
        V = V * torch.tile(norms, (dim_space, 1))

    return U, V


def construct_Laplacian_gNb_torch(gNb, dim_space, vxI=0, X=None, alphaL=10, normW=1, dataPrecision='double'):
    """
    construct_Laplacian_gNb_torch(gNb, dim_space, vxI=0, X=None, alphaL=10, normW=1, dataPrecision='double')
    construct Laplacian matrices for Laplacian spatial regularization term

    :param gNb: graph neighborhood, a 2D matrix [N, 2] storing rows and columns of non-zero elements
    :param dim_space: dimension of space (number of voxels or vertices)
    :param vxI: 0 or 1, flag for using the temporal correlation between nodes (vertex, voxel)
    :param X: fMRI data, a 2D matrix, [dim_time, dim_space]
    :param alphaL: internal hyper parameter for Laplacian regularization term
    :param normW: 1 or 2, normalization method for Laplacian matrix W
    :param dataPrecision: 'double' or 'single'
    :return: L, W, D: sparse 2D matrices [dim_space, dim_space]

    Yuncong Ma, 9/7/2023
    """

    torch_float, torch_eps = set_data_precision_torch(dataPrecision)
    np_float, np_eps = set_data_precision(dataPrecision)

    # Use numpy version to do
    # Current torch version does NOT support sparse matrix multiplication without CUDA
    if not isinstance(X, np.ndarray):
        X2 = np.array(X, dtype=np_float)
    else:
        X2 = X.astype(np_float)
    L, W, D = construct_Laplacian_gNb(gNb, dim_space, vxI, X2, alphaL, normW, dataPrecision)

    L = L.tocoo()
    # Create PyTorch sparse tensor using the COO format data
    indices = torch.tensor(np.array([L.row, L.col]), dtype=torch.long)
    values = torch.tensor(L.data, dtype=torch_float)
    L = torch.sparse_coo_tensor(indices, values, L.shape)

    D = D.tocoo()
    # Create PyTorch sparse tensor using the COO format data
    indices = torch.tensor(np.array([D.row, D.col]), dtype=torch.long)
    values = torch.tensor(D.data, dtype=torch_float)
    D = torch.sparse_coo_tensor(indices, values, D.shape)

    W = W.tocoo()
    # Create PyTorch sparse tensor using the COO format data
    indices = torch.tensor(np.array([W.row, W.col]), dtype=torch.long)
    values = torch.tensor(W.data, dtype=torch_float)
    W = torch.sparse_coo_tensor(indices, values, W.shape)

    return L, W, D


def pFN_NMF_torch(Data, gFN, gNb, maxIter=1000, minIter=30, meanFitRatio=0.1, error=1e-4, normW=1,
            Alpha=2, Beta=30, alphaS=0, alphaL=0, vxI=0, initConv=1, ard=0, eta=0, dataPrecision='double', logFile='Log_pFN_NMF.log'):
    """
    pFN_NMF_torch(Data, gFN, gNb, maxIter=1000, minIter=30,
            meanFitRatio=0.1, error=1e-4, normW=1,
            Alpha=2, Beta=30, alphaS=2, alphaL=10, initConv=1, ard=0, eta=0,
            dataPrecision='double', logFile='Log_pFN_NMF.log')
    Compute personalized FNs by spatially-regularized NMF method with group FNs as initialization

    :param Data: 2D matrix [dim_time, dim_space], numpy.ndarray or torch.Tensor. Data will be formatted to Tensor and normalized.
    :param gFN: group level FNs 2D matrix [dim_space, K], K is the number of functional networks, numpy.ndarray or torch.Tensor. gFN will be cloned
    :param gNb: graph neighborhood, a 2D matrix [N, 2] storing rows and columns of non-zero elements
    :param maxIter: maximum iteration number for multiplicative update
    :param minIter: minimum iteration in case fast convergence
    :param meanFitRatio: a 0-1 scaler, exponential moving average coefficient, used for the initialization of U when using group initialized V
    :param error: difference of cost function for convergence
    :param normW: 1 or 2, normalization method for W used in Laplacian regularization
    :param Alpha: hyper parameter for spatial sparsity
    :param Beta: hyper parameter for Laplacian sparsity
    :param alphaS: internally determined, the coefficient for spatial sparsity based Alpha, data size, K, and gNb
    :param alphaL: internally determined, the coefficient for Laplacian sparsity based Beta, data size, K, and gNb
    :param vxI: 0 or 1, flag for using the temporal correlation between nodes (vertex, voxel)
    :param initConv: flag for convergence of initialization of U
    :param ard: 0 or 1, flat for combining similar clusters
    :param eta: a hyper parameter for the ard regularization term
    :param dataPrecision: 'single' or 'float32', 'double' or 'float64'
    :param logFile: str, directory of a txt log file
    :return: U and V. U is the temporal components of pFNs, a 2D matrix [dim_time, K], and V is the spatial components of pFNs, a 2D matrix [dim_space, K]

    Yuncong Ma, 9/7/2023
    """

    # Setup data precision and eps
    torch_float, torch_eps = set_data_precision_torch(dataPrecision)

    # Transform data format if necessary
    if not isinstance(Data, torch.Tensor):
        Data = torch.tensor(Data, dtype=torch_float)
    else:
        Data = Data.type(torch_float)
    if not isinstance(gFN, torch.Tensor):
        gFN = torch.tensor(gFN, dtype=torch_float)
    else:
        gFN = gFN.type(torch_float)
    if not isinstance(error, torch.Tensor):
        error = torch.tensor(error, dtype=torch_float)
    else:
        error = error.type(torch_float)
    if not isinstance(eta, torch.Tensor):
        eta = torch.tensor(eta, dtype=torch_float)
    else:
        eta = eta.type(torch_float)
    if not isinstance(alphaS, torch.Tensor):
        alphaS = torch.tensor(alphaS, dtype=torch_float)
    else:
        alphaS = alphaS.type(torch_float)

    # check dimension of Data and gFN
    if Data.shape[1] != gFN.shape[0]:
        raise ValueError("The second dimension of Data should match the first dimension of gFn, as they are space dimension")

    K = gFN.shape[1]

    # setup log file
    logFile = open(logFile, 'a')
    print(f'\nStart NMF for pFN using PyTorch at '+time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))+'\n', file=logFile, flush=True)

    # initialization
    initV = gFN.clone()

    dim_time, dim_space = Data.shape

    # Median number of graph neighbors
    nM = np.median(np.unique(gNb[:, 0], return_counts=True)[1])

    # Use Alpha and Beta to set alphaS and alphaL if they are 0
    if alphaS == 0 and Alpha > 0:
        alphaS = torch.tensor(np.round(Alpha * dim_time / K))
    if alphaL == 0 and Beta > 0:
        alphaL = np.round(Beta * dim_time / K / nM)

    # Prepare and normalize scan
    Data = normalize_data_torch(Data, 'vp', 'vmax')
    X = Data    # Save memory

    # Construct the spatial affinity graph
    L, W, D = construct_Laplacian_gNb_torch(gNb, dim_space, vxI, X, alphaL, normW, dataPrecision)

    # Initialize V
    V = initV.clone()
    miv = torch.max(V, dim=0)[0]
    trimInd = V / torch.maximum(torch.tile(miv, (dim_space, 1)), torch_eps) < torch.tensor(5e-2)
    V[trimInd] = 0

    # Initialize U
    U = X @ V / torch.tile(torch.sum(V, dim=0), (dim_time, 1))

    U = initialize_u_torch(X, U, V, error, maxIter, minIter, meanFitRatio, initConv)

    initV = V.clone()

    # Alternative update of U and V
    # Variables

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
        U, V = normalize_u_v_torch(U, V, 1, 1)

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

        tmpDf = torch.pow(X - U @ V.T, 2)
        tmp1 = torch.sum(tmpDf)

        if alphaL > 0:
            tmp2 = V.T @ L * V.T

        L21 = alphaS * torch.sum(torch.sum(torch.sqrt(tmpl21), dim=0) / torch.maximum(torch.sqrt(torch.sum(tmpl21, dim=0)), torch_eps))
        LDf = tmp1
        LSl = torch.sum(tmp2)

        # Objective function
        LogL = L21 + ardU + LDf + LSl
        print(f"    Iter = {i}: LogL: {LogL}, dataFit: {LDf}, spaLap: {LSl}, L21: {L21}, ardU: {ardU}", file=logFile)

        # The iteration needs to meet minimum iteration number and small changes of LogL
        if i > minIter and abs(oldLogL - LogL) / torch.maximum(oldLogL, torch_eps) < error:
            break
        oldLogL = LogL.clone()

        # QC Control
        temp = mat_corr_torch(initV, V, dataPrecision=dataPrecision)
        QC_Spatial_Correspondence = torch.clone(torch.diag(temp))
        temp -= torch.diag(torch.diag(temp))
        QC_Spatial_Correspondence_Control = torch.max(temp, dim=1)[0]
        QC_Delta_Sim = torch.min(QC_Spatial_Correspondence - QC_Spatial_Correspondence_Control)
        QC_Delta_Sim = QC_Delta_Sim.cpu().numpy()

        if QC_Delta_Sim <= 0:
            flagQC = 1
            U = oldU.clone()
            V = oldV.clone()
            print(f'\n  QC: Meet QC constraint: Delta sim = {QC_Delta_Sim}', file=logFile, flush=True)
            print(f'    Use results from last iteration', file=logFile, flush=True)
            break
        else:
            oldU = U.clone()
            oldV = V.clone()
            print(f'        QC: Delta sim = {QC_Delta_Sim}', file=logFile, flush=True)

    print(f'\n Finished at '+time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))+'\n', file=logFile, flush=True)

    return U, V


def gFN_NMF_torch(Data, K, gNb, maxIter=1000, minIter=30, error=1e-6, normW=1,
            Alpha=2, Beta=30, alphaS=0, alphaL=0, vxI=0, ard=0, eta=0, nRepeat=5, dataPrecision='double', logFile='Log_pFN_NMF.log'):
    """
    gFN_NMF_torch(Data, K, gNb, maxIter=1000, minIter=30, error=1e-6, normW=1,
            Alpha=2, Beta=30, alphaS=0, alphaL=0, vxI=0, ard=0, eta=0, nRepeat=5, dataPrecision='double', logFile='Log_pFN_NMF.log')
    Compute group-level FNs using NMF method

    :param Data: 2D matrix [dim_time, dim_space], numpy.ndarray or torch.Tensor, recommend to normalize each fMRI scan before concatenate them along the time dimension
    :param K: number of FNs
    :param gNb: graph neighborhood, a 2D matrix [N, 2] storing rows and columns of non-zero elements
    :param maxIter: maximum iteration number for multiplicative update
    :param minIter: minimum iteration in case fast convergence
    :param error: difference of cost function for convergence
    :param normW: 1 or 2, normalization method for W used in Laplacian regularization
    :param Alpha: hyper parameter for spatial sparsity
    :param Beta: hyper parameter for Laplacian sparsity
    :param alphaS: internally determined, the coefficient for spatial sparsity based Alpha, data size, K, and gNb
    :param alphaL: internally determined, the coefficient for Laplacian sparsity based Beta, data size, K, and gNb
    :param vxI: flag for using the temporal correlation between nodes (vertex, voxel)
    :param ard: 0 or 1, flat for combining similar clusters
    :param eta: a hyper parameter for the ard regularization term
    :param nRepeat: Any positive integer, the number of repetition to avoid poor initialization
    :param dataPrecision: 'single' or 'float32', 'double' or 'float64'
    :param logFile: str, directory of a txt log file
    :return: gFN, 2D matrix [dim_space, K]

    Yuncong Ma, 9/24/2023
    """

    # setup log file
    logFile = open(logFile, 'a')
    print(f'\nStart NMF for gFN using PyTorch at '+time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))+'\n', file=logFile, flush=True)

    # Setup data precision and eps
    torch_float, torch_eps = set_data_precision_torch(dataPrecision)

    # Transform data format if necessary
    if not isinstance(Data, torch.Tensor):
        Data = torch.tensor(Data, dtype=torch_float)
    else:
        Data = Data.type(torch_float)
    if not isinstance(error, torch.Tensor):
        error = torch.tensor(error, dtype=torch_float)
    else:
        error = error.type(torch_float)
    if not isinstance(eta, torch.Tensor):
        eta = torch.tensor(eta, dtype=torch_float)
    else:
        eta = eta.type(torch_float)
    if not isinstance(alphaS, torch.Tensor):
        alphaS = torch.tensor(alphaS, dtype=torch_float)
    else:
        alphaS = alphaS.type(torch_float)

    # Input data size
    dim_time, dim_space = Data.shape

    # Median number of graph neighbors
    nM = np.median(np.unique(gNb[:, 0], return_counts=True)[1])

    # Use Alpha and Beta to set alphaS and alphaL if they are 0
    if alphaS == 0 and Alpha > 0:
        alphaS = torch.tensor(np.round(Alpha * dim_time / K))
    if alphaL == 0 and Beta > 0:
        alphaL = np.round(Beta * dim_time / K / nM)

    # Prepare and normalize scan
    Data = normalize_data_torch(Data, 'vp', 'vmax')
    X = Data  # Save memory

    # Construct the spatial affinity graph
    L, W, D = construct_Laplacian_gNb_torch(gNb, dim_space, vxI, X, alphaL, normW, dataPrecision)

    flag_Repeat = 0
    for repeat in range(1, 1 + nRepeat):
        flag_Repeat = 0
        print(f'\n Starting {repeat}-th repetition\n', file=logFile, flush=True)

        # Initialize U and V
        mean_X = torch.divide(torch.sum(X), torch.tensor(dim_time*dim_space))
        U = (torch.rand((dim_time, K), dtype=torch_float) + 1) * (torch.sqrt(torch.div(mean_X, K)))
        V = (torch.rand((dim_space, K), dtype=torch_float) + 1) * (torch.sqrt(torch.div(mean_X, K)))

        # Normalize data
        U, V = normalize_u_v_torch(U, V, 1, 1, dataPrecision)

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
            U, V = normalize_u_v_torch(U, V, 1, 1, dataPrecision)

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
            tmp2 = 0
            tmpl21 = torch.pow(V, 2)

            if ard > 0:
                su = torch.sum(U, dim=0)
                su[su == 0] = 1
                ardU = torch.sum(torch.log(su)) * dim_time * hyperLam

            if alphaL > 0:
                tmp2 = torch.mul(torch.matmul(V.T, L), V.T)

            L21 = torch.mul(alphaS, torch.sum(torch.div(torch.sum(torch.sqrt(tmpl21), dim=0), torch.maximum(torch.sqrt(torch.sum(tmpl21, dim=0)), torch_eps))))
            # LDf = data_fitting_error(X, U, V, 0, 1)
            LDf = torch.sum(torch.pow(torch.subtract(X, torch.matmul(U, V.T)), 2))
            LSl = torch.sum(tmp2)

            # Objective function
            LogL = L21 + LDf + LSl + ardU
            print(f"    Iter = {i}: LogL: {LogL}, dataFit: {LDf}, spaLap: {LSl}, L21: {L21}, ardU: {ardU}", file=logFile, flush=True)

            if i > 1 and i < minIter and abs(oldLogL - LogL) / torch.maximum(oldLogL, torch_eps) < error:
                flag_Repeat = 1
                print('\n Iteration stopped before the minimum iteration number. The results might be poor.\n', file=logFile, flush=True)
                break
            elif i > minIter and abs(oldLogL - LogL) / torch.maximum(oldLogL, torch_eps) < error:
                break
            oldLogL = LogL.clone()
        if flag_Repeat == 0:
            break

    if flag_Repeat == 1:
        print('\n All repetition stopped before the minimum iteration number. The final results might be poor\n', file=logFile, flush=True)

    gFN = V
    print(f'\nFinished at '+time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))+'\n', file=logFile, flush=True)
    return gFN


def gFN_fusion_NCut_torch(gFN_BS, K, NCut_MaxTrial=100, dataPrecision='double', logFile='Log_gFN_fusion_NCut'):
    """
    gFN_fusion_NCut_torch(gFN_BS, K, NCut_MaxTrial=100, dataPrecision='double')
    Fuses FN results to generate representative group-level FNs

    :param gFN_BS: FNs obtained from bootstrapping method, FNs are concatenated along the K dimension
    :param K: Number of FNs, not the total number of FNs obtained from bootstrapping
    :param NCut_MaxTrial: Max number trials for NCut method
    :param dataPrecision: 'double' or 'single'
    :param logFile: str, directory of a txt log file
    :return: gFNs, 2D matrix [dim_space, K]

    Yuncong Ma, 9/6/2023
    """

    # setup log file
    logFile = open(logFile, 'a')
    print(f'\nStart NCut for gFN fusion using PyTorch at '+time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))+'\n', file=logFile, flush=True)

    # Setup data precision and eps
    torch_float, torch_eps = set_data_precision_torch(dataPrecision)

    if not isinstance(gFN_BS, torch.Tensor):
        gFN_BS = torch.tensor(gFN_BS, dtype=torch_float)
    else:
        gFN_BS = gFN_BS.type(torch_float)

    # clustering by NCut

    # Get similarity between samples
    corrVal = mat_corr_torch(gFN_BS, dataPrecision=dataPrecision)  # similarity between FNs, [K * n_BS, K * n_BS]
    corrVal[torch.isnan(corrVal)] = -1
    nDis = 1 - corrVal  # Transform Pearson correlation to non-negative values similar to distance
    triuInd = torch.triu(torch.ones(nDis.shape), 1)  # Index of upper triangle
    nDisVec = nDis[triuInd == 1]  # Take the values in upper triangle
    # Make all distance non-negative and normalize their distribution
    nW = torch.exp(-torch.pow(nDis, 2) / (torch.pow(torch.median(nDisVec), 2)))  # Transform distance values using exp(-X/std^2) with std as the median value
    nW[torch.isnan(nW)] = 0
    sumW = torch.sum(nW, dim=1)  # total distance for each FN
    sumW[sumW == 0] = 1  # In case two FNs are the same
    # Construct Laplacian matrix
    D = torch.diag(sumW)
    L = torch.sqrt(torch.linalg.inv(D)) @ nW @ torch.linalg.inv(torch.sqrt(D))  # A way to normalize nW based on the total distance of each FN
    L = (L + L.T) / 2  # Ensure L is symmetric. Computation error may result in asymmetry

    # Get first K eigenvectors, sign of vectors may be different to MATLAB results
    eigenvalues, eigenvectors = torch.linalg.eigh(L.type(torch_float))
    # Sort by eigenvalues and select the K largest
    sorted_indices = torch.argsort(eigenvalues, descending=True)[:K]
    eVal = eigenvalues[sorted_indices]
    Ev = eigenvectors[:, sorted_indices]
    Ev = torch.real(Ev)
    # Correct the sign of eigenvectors to make them same as derived from MATLAB
    temp = torch.sign(torch.sum(Ev, dim=0))  # Use the total value of each eigenvector to reset its sign
    temp[temp == 0.0] = 1.0
    Ev = Ev * torch.tile(temp, (Ev.shape[0], 1))  # Reset the sign of each eigenvector
    normvect = torch.sqrt(torch.diag(Ev @ Ev.T))  # Get the norm of each eigenvector
    normvect[normvect == 0.0] = 1  # Incase all 0 eigenvector
    Ev = torch.linalg.solve(torch.diag(normvect), Ev)  # Use linear solution to normalize Ev satisfying normvect * Ev_new = Ev_old

    # Multiple trials to get reproducible results
    Best_C = []
    Best_NCutValue = torch.inf
    for i in range(1,NCut_MaxTrial+1):
        print(f'    Iter = ' + str(i), file=logFile)
        EigenVectors = Ev
        n, k = EigenVectors.shape  # n is K * n_BS, k is K

        vm = torch.sqrt(torch.sum(EigenVectors**2, dim=1, keepdims=True))  # norm of each row
        EigenVectors = EigenVectors / torch.tile(vm, (1, k))  # normalize eigenvectors to ensure each FN vector's norm = 1

        R = torch.zeros((k, k), dtype=torch_float)
        ps = torch.randint(0, n, (1,))  # Choose a random row in eigenvectors
        R[:, 0] = EigenVectors[ps, :]  # This randomly selected row in eigenvectors is used as an initial center

        c_index = torch.zeros(n)
        c = torch.zeros(n, dtype=torch_float)  # Total distance to different rows in R [K * n_BS, 1]
        c_index[0] = ps  # Store the index of selected samples
        c[ps] = torch.inf

        for j in range(2, k+1):  # Find another K-1 rows in eigenvectors which have the minimum similarity to previous selected rows, similar to initialization in k++
            c += torch.abs(EigenVectors @ R[:, j-2])
            ps = torch.argmin(c)
            c_index[j-1] = ps
            c[ps] = torch.inf
            R[:, j-1] = EigenVectors[ps, :]

        lastObjectiveValue = 0
        exitLoop = 0
        nbIterationsDiscretisation = 0
        nbIterationsDiscretisationMax = 20

        while exitLoop == 0:
            nbIterationsDiscretisation += 1

            EigenVectorsR = EigenVectors @ R
            n, k = EigenVectorsR.shape
            J = torch.argmax(EigenVectorsR, dim=1)  # Assign each sample to K centers of R based on highest similarity
            Indice = torch.stack((torch.arange(n).type(torch.int32), J.type(torch.int32)))
            EigenvectorsDiscrete = torch.sparse_coo_tensor(Indice, torch.ones(n, dtype=torch_float), (n, k))  # Generate a 0-1 matrix with each row containing only one 1

            U, S, Vh = torch.linalg.svd(EigenvectorsDiscrete.T @ EigenVectors, full_matrices=False)  # Economy-size decomposition

            S = torch.diag(S)  # To match MATLAB svd results
            V = Vh.T  # To match MATLAB svd results
            NcutValue = 2 * (n - torch.trace(S))

            # escape the loop when converged or meet max iteration
            if abs(NcutValue - lastObjectiveValue) < torch_eps or nbIterationsDiscretisation > nbIterationsDiscretisationMax:
                exitLoop = 1
                print(f'    Reach stop criterion of NCut, NcutValue = '+str(NcutValue.numpy())+'\n', file=logFile, flush=True)
            else:
                print(f'    NcutValue = '+str(NcutValue.numpy()), file=logFile)
                lastObjectiveValue = NcutValue
                R = V @ U.T  # Update R which stores the new centers

        C = torch.argmax(EigenvectorsDiscrete.to_dense(), dim=1)  # Assign each sample to K centers in R

        if len(torch.unique(C)) < K:  # Check whether there are empty results
            print(f'    Found empty results in iteration '+str(i)+'\n', file=logFile, flush=True)
        else:  # Update the best result
            if NcutValue < Best_NCutValue:
                Best_NCutValue = NcutValue
                Best_C = C

    if len(set(Best_C)) < K:  # In case even the last trial has empty results
        raise ValueError('  Cannot generate non-empty gFNs\n')
        Flag = 1
        Message = "Cannot generate non-empty FN"

    print(f'Best NCut value = '+str(Best_NCutValue.numpy())+'\n', file=logFile, flush=True)

    # Get centroid
    C = Best_C
    gFN = torch.zeros((gFN_BS.shape[0], K))
    for ki in range(K):
        if torch.sum(C == ki) > 1:
            candSet = gFN_BS[:, C == ki]  # Get the candidate set of FNs assigned to cluster ki
            corrW = torch.abs(mat_corr_torch(candSet, dataPrecision=dataPrecision))  # Get the similarity between candidate FNs
            corrW[torch.isnan(corrW)] = 0
            mInd = torch.argmax(torch.sum(corrW, dim=0), dim=0)  # Find the FN with the highest total similarity to all other FNs
            gFN[:, ki] = candSet[:, mInd]
        elif torch.sum(C == ki) == 1:
            mInd = int((C == ki).nonzero(as_tuple=True)[0])
            gFN[:, ki] = gFN_BS[:, mInd]

    gFN = gFN / torch.maximum(torch.tile(torch.max(gFN, dim=0)[0], (gFN.shape[0], 1)), torch_eps)  # Normalize each FN by its max value
    print(f'\nFinished at '+time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))+'\n', file=logFile, flush=True)

    return gFN


def compute_gNb_torch(Brain_Template, logFile=None):
    """
    compute_gNb_torch(Brain_Template, logFile=None)
    Prepare a graph neighborhood variable, using indices as its sparse representation

    :param Brain_Template: a structure variable with keys 'Data_Type', 'Data_Format', 'Shape', 'Mask'.
        If Brain_Template.Data_Type is 'Surface', Shape contains L and R, with vertices and faces as sub keys. Mask contains L and R.
        If Brain_Template.Data_Type is 'Volume', Shape is None, Mask is a 3D 0-1 matrix, Overlay_Image is a 3D matrix
    :param logFile:
    :return: gNb: a 2D matrix [N, 2], which labels the non-zero elements in a graph. Index starts from 1

    Yuncong Ma, 9/13/2023
    """

    gNb = compute_gNb(Brain_Template, logFile)

    return gNb


def bootstrap_scan_torch(dir_output: str, file_scan: str, file_subject_ID: str, file_subject_folder: str, file_group_ID=None, combineScan=0,
                         samplingMethod='Subject', sampleSize=10, nBS=50, logFile=None):
    """
    bootstrap_scan_torch(dir_output: str, file_scan: str, file_subject_ID: str, file_subject_folder: str, file_group=None, BS=10, N_BS=50, combineFLag=0, samplingMethod='Subject', logFile=None)
    prepare bootstrapped scan file lists

    :param dir_output: directory of a folder to store bootstrapped files
    :param file_scan: a txt file that stores directories of all fMRI scans
    :param file_subject_ID: a txt file that store subject ID information corresponding to fMRI scan in file_scan
    :param file_subject_folder: a txt file that store subject folder names corresponding to fMRI scans in file_scan
    :param file_group_ID: a txt file that store group information corresponding to fMRI scan in file_scan
    :param combineScan: 0 or 1, whether to combine multiple fMRI scans for each subject
    :param samplingMethod: 'Subject' or 'Group_Subject'. Uniform sampling based subject ID, or group and then subject ID
    :param sampleSize: number of subjects selected for each bootstrapping run
    :param nBS: number of runs for bootstrap
    :param logFile: directory of a txt file
    :return: None

    Yuncong Ma, 9/21/2023
    """

    bootstrap_scan(dir_output, file_scan, file_subject_ID=file_subject_ID, file_subject_folder=file_subject_folder,
                   file_group_ID=file_group_ID, combineScan=combineScan, samplingMethod=samplingMethod,
                   sampleSize=sampleSize, nBS=nBS, logFile=logFile)


def setup_NMF_setting_torch(dir_pnet_result: str, K=17, Combine_Scan=False, Compute_gFN=True, samplingMethod='Subject', sampleSize=10, nBS=50, maxIter=1000, minIter=30, meanFitRatio=0.1, error=1e-6,
                      normW=1, Alpha=2, Beta=30, alphaS=0, alphaL=0, vxI=0, ard=0, eta=0, nRepeat=5, Parallel=False, Computation_Mode='CPU', N_Thread=1, dataPrecision='double'):
    """
    setup_NMF_setting_torch(dir_pnet_result: str, K=17, Combine_Scan=False, Compute_gFN=True, samplingMethod='Subject', sampleSize=10, nBS=50, maxIter=1000, minIter=30, meanFitRatio=0.1, error=1e-6,
                      normW=1, Alpha=2, Beta=30, alphaS=0, alphaL=0, vxI=0, ard=0, eta=0, nRepeat=5, Parallel=False, Computation_Mode='CPU', N_Thread=1, dataPrecision='double')
    Setup the setting for NMF-based method to compute gFNs and pFNs

    :param dir_pnet_result: directory of the pNet result folder
    :param K: number of FNs
    :param Combine_Scan: False or True, whether to combine multiple scans for the same subject
    :param Compute_gFN: True or False, whether to compute gFNs from the provided data or load a precomputed gFN set
    :param samplingMethod: 'Subject' or 'Group_Subject'. Uniform sampling based subject ID, or group and then subject ID
    :param sampleSize: number of subjects selected for each bootstrapping run
    :param nBS: number of runs for bootstrap
    :param maxIter: maximum iteration number for multiplicative update
    :param minIter: minimum iteration in case fast convergence
    :param meanFitRatio: a 0-1 scaler, exponential moving average coefficient, used for the initialization of U when using group initialized V
    :param error: difference of cost function for convergence
    :param normW: 1 or 2, normalization method for W used in Laplacian regularization
    :param Alpha: hyper parameter for spatial sparsity
    :param Beta: hyper parameter for Laplacian sparsity
    :param alphaS: internally determined, the coefficient for spatial sparsity based Alpha, data size, K, and gNb
    :param alphaL: internally determined, the coefficient for Laplacian sparsity based Beta, data size, K, and gNb
    :param vxI: flag for using the temporal correlation between nodes (vertex, voxel)
    :param ard: 0 or 1, flat for combining similar clusters
    :param eta: a hyper parameter for the ard regularization term
    :param nRepeat: Any positive integer, the number of repetition to avoid poor initialization
    :param Parallel: False or True, whether to enable parallel computation
    :param Computation_Mode: 'CPU'
    :param N_Thread: positive integers, used for parallel computation
    :param dataPrecision: 'double' or 'single'
    :return: setting: a structure

    Yuncong Ma, 9/21/2023
    """

    setting = setup_NMF_setting(dir_pnet_result, K=K, Combine_Scan=Combine_Scan, Compute_gFN=Compute_gFN,
                                samplingMethod=samplingMethod, sampleSize=sampleSize, nBS=nBS,
                                maxIter=maxIter, minIter=minIter, meanFitRatio=meanFitRatio, error=error, normW=normW,
                                Alpha=Alpha, Beta=Beta, alphaS=alphaS, alphaL=alphaL, vxI=vxI, ard=ard, eta=eta,
                                nRepeat=nRepeat,
                                Parallel=Parallel, Computation_Mode=Computation_Mode, N_Thread=N_Thread,
                                dataPrecision=dataPrecision)

    return setting


def setup_pFN_folder_torch(dir_pnet_result: str):
    """
    setup_pFN_folder_torch(dir_pnet_result: str)
    Setup sub-folders in Personalized_FN to

    :param dir_pnet_result: directory of the pNet result folder
    :return: list_subject_folder_unique: unique subject folder array for getting sub-folders in Personalized_FN

    Yuncong Ma, 9/21/2023
    """

    list_subject_folder_unique = setup_pFN_folder(dir_pnet_result)
    return list_subject_folder_unique


def run_FN_Computation_torch(dir_pnet_result: str):
    """
    run_FN_Computation_torch(dir_pnet_result: str)
    run the FN Computation module with settings ready in Data_Input and FN_Computation

    :param dir_pnet_result: directory of pNet result folder

    Yuncong Ma, 9/24/2023
    """

    # get directories of sub-folders
    dir_pnet_dataInput, dir_pnet_FNC, dir_pnet_gFN, dir_pnet_pFN, _, _ = setup_result_folder(dir_pnet_result)

    # load settings for data input and FN computation
    if not os.path.isfile(os.path.join(dir_pnet_dataInput, 'Setting.json')):
        raise ValueError('Cannot find the setting json file in folder Data_Input')
    if not os.path.isfile(os.path.join(dir_pnet_FNC, 'Setting.json')):
        raise ValueError('Cannot find the setting json file in folder FN_Computation')
    settingDataInput = load_json_setting(os.path.join(dir_pnet_dataInput, 'Setting.json'))
    settingFNC = load_json_setting(os.path.join(dir_pnet_FNC, 'Setting.json'))
    setting = {'Data_Input': settingDataInput, 'FN_Computation': settingFNC}

    # load basic settings
    dataType = setting['Data_Input']['Data_Type']
    dataFormat = setting['Data_Input']['Data_Format']

    # load Brain Template
    Brain_Template = load_brain_template(os.path.join(dir_pnet_dataInput, 'Brain_Template.json'))
    if dataType == 'Volume':
        Brain_Mask = Brain_Template['Brain_Mask']
    else:
        Brain_Mask = None

    # ============== gFN Computation ============== #
    # Start computation using SP-NMF
    if setting['FN_Computation']['Method'] == 'SR-NMF':

        if setting['FN_Computation']['Group_FN']['Compute_gFN']:
            # 2 steps
            # step 1 ============== bootstrap
            # sub-folder in FNC for storing bootstrapped results
            dir_pnet_BS = os.path.join(dir_pnet_FNC, 'BootStrapping')
            if not os.path.exists(dir_pnet_BS):
                os.makedirs(dir_pnet_BS)
            # Log
            logFile = os.path.join(dir_pnet_BS, 'Log.log')

            # Generate additional parameters
            gNb = compute_gNb_torch(Brain_Template)
            scipy.io.savemat(os.path.join(dir_pnet_FNC, 'gNb.mat'), {'gNb': gNb})
            # Input files
            file_scan = os.path.join(dir_pnet_dataInput, 'Scan_List.txt')
            file_subject_ID = os.path.join(dir_pnet_dataInput, 'Subject_ID.txt')
            file_subject_folder = os.path.join(dir_pnet_dataInput, 'Subject_Folder.txt')
            file_group_ID = os.path.join(dir_pnet_dataInput, 'Group_ID.txt')
            if not os.path.exists(file_group_ID):
                file_group = None
            # Parameters
            combineScan = setting['FN_Computation']['Combine_Scan']
            samplingMethod = setting['FN_Computation']['Group_FN']['BootStrap']['samplingMethod']
            sampleSize = setting['FN_Computation']['Group_FN']['BootStrap']['sampleSize']
            nBS = setting['FN_Computation']['Group_FN']['BootStrap']['nBS']

            # create scan lists for bootstrap
            bootstrap_scan_torch(dir_pnet_BS, file_scan, file_subject_ID, file_subject_folder,
                                 file_group_ID=file_group_ID, combineScan=combineScan,
                                 samplingMethod=samplingMethod, sampleSize=sampleSize, nBS=nBS, logFile=logFile)

            # Parameters
            K = setting['FN_Computation']['K']
            maxIter = setting['FN_Computation']['Group_FN']['maxIter']
            minIter = setting['FN_Computation']['Group_FN']['minIter']
            error = setting['FN_Computation']['Group_FN']['error']
            normW = setting['FN_Computation']['Group_FN']['normW']
            Alpha = setting['FN_Computation']['Group_FN']['Alpha']
            Beta = setting['FN_Computation']['Group_FN']['Beta']
            alphaS = setting['FN_Computation']['Group_FN']['alphaS']
            alphaL = setting['FN_Computation']['Group_FN']['alphaL']
            vxI = setting['FN_Computation']['Group_FN']['vxI']
            ard = setting['FN_Computation']['Group_FN']['ard']
            eta = setting['FN_Computation']['Group_FN']['eta']
            nRepeat = setting['FN_Computation']['Group_FN']['nRepeat']
            dataPrecision = setting['FN_Computation']['Computation']['dataPrecision']

            # NMF on bootstrapped subsets
            for rep in range(1, 1+nBS):
                # log file
                logFile = os.path.join(dir_pnet_BS, str(rep), 'Log.log')
                # load data
                file_scan_list = os.path.join(dir_pnet_BS, str(rep), 'Scan_List.txt')
                Data = load_fmri_scan(file_scan_list, dataType=dataType, dataFormat=dataFormat, Reshape=True, Brain_Mask=Brain_Mask,
                                      Normalization='vp-vmax', logFile=logFile)
                # perform NMF
                FN_BS = gFN_NMF_torch(Data, K, gNb, maxIter=maxIter, minIter=minIter, error=error, normW=normW,
                                      Alpha=Alpha, Beta=Beta, alphaS=alphaS, alphaL=alphaL, vxI=vxI, ard=ard, eta=eta,
                                      nRepeat=nRepeat, dataPrecision=dataPrecision, logFile=logFile)
                # save results
                FN_BS = reshape_FN(FN_BS.numpy(), dataType=dataType, Brain_Mask=Brain_Mask)
                sio.savemat(os.path.join(dir_pnet_BS, str(rep), 'FN.mat'), {"FN": FN_BS})

            # step 2 ============== fuse results
            # Generate gFNs
            FN_BS = np.empty(nBS, dtype=np.ndarray)
            # load bootstrapped results
            for rep in range(1, nBS+1):
                FN_BS[rep-1] = np.array(reshape_FN(load_matlab_single_array(os.path.join(dir_pnet_BS, str(rep), 'FN.mat')), dataType=dataType, Brain_Mask=Brain_Mask))
            gFN_BS = np.concatenate(FN_BS, axis=1)
            # log
            logFile = os.path.join(dir_pnet_gFN, 'Log.log')
            # Fuse bootstrapped results
            gFN = gFN_fusion_NCut_torch(gFN_BS, K, logFile=logFile)
            # output
            gFN = reshape_FN(gFN.numpy(), dataType=dataType, Brain_Mask=Brain_Mask)
            sio.savemat(os.path.join(dir_pnet_gFN, 'FN.mat'), {"FN": gFN})
        # ============================================= #

        # ============== pFN Computation ============== #
        # load precomputed gFNs
        gFN = load_matlab_single_array(os.path.join(dir_pnet_gFN, 'FN.mat'))
        # additional parameter
        gNb = load_matlab_single_array(os.path.join(dir_pnet_FNC, 'gNb.mat'))
        # reshape to 2D if required
        gFN = reshape_FN(gFN, dataType=dataType, Brain_Mask=Brain_Mask)
        # setup folders in Personalized_FN
        list_subject_folder = setup_pFN_folder(dir_pnet_result)
        N_Scan = len(list_subject_folder)
        for i in range(1, N_Scan+1):
            # print(' Processing ' + str(i))
            dir_pnet_pFN_indv = os.path.join(dir_pnet_pFN, list_subject_folder[i-1])
            # parameter
            maxIter = setting['FN_Computation']['Personalized_FN']['maxIter']
            minIter = setting['FN_Computation']['Personalized_FN']['minIter']
            meanFitRatio = setting['FN_Computation']['Personalized_FN']['meanFitRatio']
            error = setting['FN_Computation']['Personalized_FN']['error']
            normW = setting['FN_Computation']['Personalized_FN']['normW']
            Alpha = setting['FN_Computation']['Personalized_FN']['Alpha']
            Beta = setting['FN_Computation']['Personalized_FN']['Beta']
            alphaS = setting['FN_Computation']['Personalized_FN']['alphaS']
            alphaL = setting['FN_Computation']['Personalized_FN']['alphaL']
            vxI = setting['FN_Computation']['Personalized_FN']['vxI']
            ard = setting['FN_Computation']['Personalized_FN']['ard']
            eta = setting['FN_Computation']['Personalized_FN']['eta']
            dataPrecision = setting['FN_Computation']['Computation']['dataPrecision']
            # log file
            logFile = os.path.join(dir_pnet_pFN_indv, 'Log.log')
            # load data
            Data = load_fmri_scan(os.path.join(dir_pnet_pFN_indv, 'Scan_List.txt'),
                                  dataType=dataType, dataFormat=dataFormat,
                                  Reshape=True, Brain_Mask=Brain_Mask, logFile=logFile)
            # perform NMF
            TC, pFN = pFN_NMF_torch(Data, gFN, gNb, maxIter=maxIter, minIter=minIter, meanFitRatio=meanFitRatio,
                                    error=error, normW=normW,
                                    Alpha=Alpha, Beta=Beta, alphaS=alphaS, alphaL=alphaL,
                                    vxI=vxI,  ard=ard, eta=eta,
                                    dataPrecision=dataPrecision, logFile=logFile)
            pFN = pFN.numpy()
            TC = TC.numpy()

            # output
            pFN = reshape_FN(pFN, dataType=dataType, Brain_Mask=Brain_Mask)
            sio.savemat(os.path.join(dir_pnet_pFN_indv, 'FN.mat'), {"FN": pFN})
            sio.savemat(os.path.join(dir_pnet_pFN_indv, 'TC.mat'), {"TC": TC})
        # ============================================= #

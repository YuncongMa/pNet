#########################################
# Packages
import os
import numpy
import pNet
import torch


#########################################
# Example data from HCP surf Test_FN17
gFN = pNet.load_matlab_array(os.path.join(pNet.Example.HCP_surf.dir_pnet, 'Group_FN', 'FN.mat'), 'FN')
# gNb uses index starting from 1, but it starts from 0 in Python
gNb = pNet.load_matlab_single_array(os.path.join(pNet.Example.HCP_surf.dir_pnet, 'FN_Computation', 'gNb.mat'))

Setting = pNet.load_json_setting(os.path.join(pNet.Example.HCP_surf.dir_pnet, 'FN_Computation', 'Setting.json'))

#########################################
# Setup parameters
K = gFN.shape[1]
maxIter = 1000
minIter = 30
meanFitRatio = 0.1
error = 1e-4
normW = 1
eta = Setting['PersonalizedFN']['eta']
alphaS21 = Setting['PersonalizedFN']['alphaS21']
alphaS1 = 0
alphaL = Setting['PersonalizedFN']['alphaL']
initConv = 1
ard = Setting['PersonalizedFN']['ard']
dataPrecision = 'double'
file_log = 'Log_pFN_NMF_torch.log'


#########################################
# Setup data precision and eps
torch_float, torch_eps = pNet.set_data_precision_torch(dataPrecision)

gFN = torch.tensor(gFN, dtype=torch_float)
error = torch.tensor(error, dtype=torch_float)
eta = torch.tensor(eta, dtype=torch_float)
alphaS21 = torch.tensor(alphaS21, dtype=torch_float)
alphaS1 = torch.tensor(alphaS1, dtype=torch_float)

dim_space = 59412
dim_time = 400

n_indice = 0
for vi in range(dim_space):
    n_indice += gNb[vi][0].shape[1]
print(n_indice)

indices = torch.zeros((2, n_indice*2), dtype=torch_float)
count = 0
for vi in range(dim_space):
    for ni in range(gNb[vi][0].shape[1]):
        nei = gNb[vi][0][0, ni] - 1  # Python index
        indices[0, count] = vi
        indices[1, count] = nei
        count += 1
        indices[1, count] = vi
        indices[0, count] = nei
        count += 1
W = torch.sparse_coo_tensor(indices, torch.ones((count,)), (dim_space, dim_space), dtype=torch_float)

# Defining other matrices
DCol = W.sum(dim=1).to_dense()
indices = torch.cat(((torch.arange(0, dim_space)).reshape(-1, 1), (torch.arange(0, dim_space)).reshape(-1, 1)), dim=1)
D = torch.sparse_coo_tensor(indices.T, DCol, (dim_space, dim_space), dtype=torch_float)
L = D - W
print(L.shape)
print(D.shape)
print(W.shape)

if isinstance(alphaL, torch.Tensor):
    alphaL = torch.tensor(alphaL)

if normW > 0:
    D_mhalf = torch.sparse_coo_tensor(indices.T, torch.pow(DCol, -0.5), (dim_space, dim_space), dtype=torch_float)
    L = D_mhalf @ L @ D_mhalf * alphaL
    W = D_mhalf @ W @ D_mhalf * alphaL
    D = D_mhalf @ D @ D_mhalf * alphaL

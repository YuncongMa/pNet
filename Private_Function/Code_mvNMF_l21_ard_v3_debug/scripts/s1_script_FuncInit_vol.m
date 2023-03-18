clc;
clear;

%addpath(genpath('C:\Users\LiHon\Google Drive\Code\Inhouse\ongoing\Code_mvNMF_l21_ard_v3_release\Release'));
%addpath(genpath('D:\Google_drive\Code\Inhouse\ongoing\Code_mvNMF_l21_ard_v3_release\Release'));

sbjListFile = '/cbica/home/zhouz/projects/istaging/LiHM_NMF/sublis_test.txt';
maskFile = '/cbica/home/zhouz/projects/istaging/LiHM_NMF/atlas/MNI-maxprob-thr50-2mm-mask.nii.gz';
prepDataFile = '/cbica/home/zhouz/projects/istaging/LiHM_NMF/test_CreatePrepData.mat';
outDir = ['/cbica/home/zhouz/projects/istaging/LiHM_NMF/BLSA_test', filesep, 'init_r2'];
spaR = 1;
vxI = 0;
ard = 1;
iterNum = 1000;
K = 17;
tNum = 180;
alpha = 2;
beta = 10;
resId = 'blsa_vol';

deployFuncInit_vol(sbjListFile,maskFile,prepDataFile,outDir,spaR,vxI,ard,iterNum,K,tNum,alpha,beta,resId);

clc;
clear;

%addpath(genpath('C:\Users\LiHon\Google Drive\Code\Inhouse\ongoing\Code_mvNMF_l21_ard_v3_release\Release'));

sbjListFile = '/cbica/home/zhouz/projects/istaging/LiHM_NMF/sublist_test.txt';
maskFile = '/cbica/home/zhouz/projects/istaging/LiHM_NMF/atlas/MNI-maxprob-thr50-2mm-mask.nii.gz';
prepDataFile = '/cbica/home/zhouz/projects/istaging/LiHM_NMF/test_CreatePrepData.mat';
outDir = ['/cbica/home/zhouz/projects/istaging/LiHM_NMF/BLSA_test', filesep, 'res_test'];

resId = 'blsa_vol';
%initName = 'F:\data\test_brain_decomp_release\res\pnc_vol\init\pnc_vol_num5_comp17_S1_69_L_17_spaR_1_vxInfo_0_ard_1\init.mat';
initName = '/cbica/home/zhouz/projects/istaging/LiHM_NMF/BLSA_test/robust_init/init.mat';
K = 17;
alphaS21 = 2;
alphaL = 10;
vxI = 0;
spaR = 1;
ard = 1;
eta = 1;
iterNum = 30;
calcGrp = 1;
parforOn = 0;

deployFuncMvnmfL21p1_func_vol(sbjListFile,maskFile,prepDataFile,outDir,resId,initName,K,alphaS21,alphaL,vxI,spaR,ard,eta,iterNum,calcGrp,parforOn);

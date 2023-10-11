clc;
clear;

addpath(genpath('C:\Users\LiHon\Google Drive\Code\Inhouse\ongoing\Code_mvNMF_l21_ard_v3_release\Release'));

sbjListFile = 'F:\data\test_brain_decomp_release\data\hcp\hcp_sbjLst.txt';
wbPath = 'C:\work_code\download\workbench\bin_windows64\wb_command.exe';
surfL = 'F:\data\HCP\HCP-20-subjects\100307\MNINonLinear\fsaverage_LR32k\100307.L.inflated.32k_fs_LR.surf.gii';
surfR = 'F:\data\HCP\HCP-20-subjects\100307\MNINonLinear\fsaverage_LR32k\100307.R.inflated.32k_fs_LR.surf.gii';
surfML = 'F:\data\HCP\HCP-20-subjects\100307\MNINonLinear\fsaverage_LR32k\100307.L.atlasroi.32k_fs_LR.shape.gii';
surfMR = 'F:\data\HCP\HCP-20-subjects\100307\MNINonLinear\fsaverage_LR32k\100307.R.atlasroi.32k_fs_LR.shape.gii';
prepDataFile = 'F:\data\test_brain_decomp_release\res\hcp\prepData.mat';
outDir = ['F:\data\test_brain_decomp_release\res\hcp', filesep, 'init'];
spaR = 1;
vxI = 0;
ard = 0;
iterNum = 1000;
K = 17;
tNum = 1200;
alpha = 2;
beta = 10;
resId = 'hcp';

deployFuncInit_surf_hcp(sbjListFile,wbPath,surfL,surfR,surfML,surfMR,prepDataFile,outDir,spaR,vxI,ard,iterNum,K,tNum,alpha,beta,resId); %,svFileL,svFileR);

clc;
clear;

addpath(genpath('D:\Google_drive\Code\Inhouse\ongoing\Code_mvNMF_l21_ard_v3_release\Release'));
%addpath(genpath('C:\Users\LiHon\Google Drive\Code\Inhouse\ongoing\Code_mvNMF_l21_ard_v3_release\Release'));

sbjListFile = 'F:\data\test_brain_decomp_release\data\hcp\hcp_sbjLst.txt';
%wbPath = 'C:\work_code\download\workbench\bin_windows64\wb_command.exe';
wbPath = 'D:\Code\Download\workbench\bin_windows64\wb_command.exe';
prepDataFile = 'F:\data\test_brain_decomp_release\res\hcp\prepData.mat';
outDir = ['F:\data\test_brain_decomp_release\res\hcp', filesep, 'res'];

resId = 'hcp';
initName = 'F:\data\test_brain_decomp_release\res\hcp\init\hcp_num5_comp17_S1_706_L_294_spaR_1_vxInfo_0_ard_0\init.mat';
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

deployFuncMvnmfL21p1_func_surf_hcp(sbjListFile,wbPath,prepDataFile,outDir,resId,initName,K,alphaS21,alphaL,vxI,spaR,ard,eta,iterNum,calcGrp,parforOn);

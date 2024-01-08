clc;
clear;

addpath(genpath('C:\Users\LiHon\Google Drive\Code\Inhouse\ongoing\Code_mvNMF_l21_ard_v3_release\Release'));

sbjListFile = 'F:\data\test_brain_decomp_release\data\pnc_fs\pnc_fs_sbjLst.txt';
medialWallFileL = 'C:\Users\LiHon\Google Drive\Code\Inhouse\ongoing\Code_mvNMF_l21_ard_v3_release\Release\lib\freesurfer\subjects\fsaverage5\label\lh.Medial_wall.label';
medialWallFileR = 'C:\Users\LiHon\Google Drive\Code\Inhouse\ongoing\Code_mvNMF_l21_ard_v3_release\Release\lib\freesurfer\subjects\fsaverage5\label\rh.Medial_wall.label';
prepDataFile = 'F:\data\test_brain_decomp_release\res\pnc_fs\prepData.mat';
outDir = ['F:\data\test_brain_decomp_release\res\pnc_fs', filesep, 'res'];

resId = 'pnc_fs';
initName = 'F:\data\test_brain_decomp_release\res\pnc_fs\init\pnc_fs_num5_comp17_S1_71_L_29_spaR_1_vxInfo_0_ard_0\init.mat';
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

deployFuncMvnmfL21p1_func_surf_fs(sbjListFile,medialWallFileL,medialWallFileR,prepDataFile,outDir,resId,initName,K,alphaS21,alphaL,vxI,spaR,ard,eta,iterNum,calcGrp,parforOn);

clc;
clear;

addpath(genpath('C:\Users\LiHon\Google Drive\Code\Inhouse\ongoing\Code_mvNMF_l21_ard_v3_release\Release'));

sbjListFile = 'F:\data\test_brain_decomp_release\data\pnc_fs\pnc_fs_sbjLst.txt';
surfL = 'C:\Users\LiHon\Google Drive\Code\Inhouse\ongoing\Code_mvNMF_l21_ard_v3_release\Release\lib\freesurfer\subjects\fsaverage5\surf\lh.pial';
surfR = 'C:\Users\LiHon\Google Drive\Code\Inhouse\ongoing\Code_mvNMF_l21_ard_v3_release\Release\lib\freesurfer\subjects\fsaverage5\surf\rh.pial';
surfML = 'C:\Users\LiHon\Google Drive\Code\Inhouse\ongoing\Code_mvNMF_l21_ard_v3_release\Release\lib\freesurfer\subjects\fsaverage5\label\lh.Medial_wall.label';
surfMR = 'C:\Users\LiHon\Google Drive\Code\Inhouse\ongoing\Code_mvNMF_l21_ard_v3_release\Release\lib\freesurfer\subjects\fsaverage5\label\rh.Medial_wall.label';
prepDataFile = 'F:\data\test_brain_decomp_release\res\pnc_fs\prepData.mat';
outDir = ['F:\data\test_brain_decomp_release\res\pnc_fs', filesep, 'init'];
spaR = 1;
vxI = 0;
ard = 0;
iterNum = 1000;
K = 17;
tNum = 120;
alpha = 2;
beta = 10;
resId = 'pnc_fs';

deployFuncInit_surf_fs(sbjListFile,surfL,surfR,surfML,surfMR,prepDataFile,outDir,spaR,vxI,ard,iterNum,K,tNum,alpha,beta,resId);

clc;
clear;

%addpath(genpath(pwd));
%addpath(genpath('D:\Google_drive\Code\Inhouse\ongoing\Code_mvNMF_l21_ard_v3_release\Release'));


% for volumetric data
maskFile = '/cbica/home/zhouz/projects/istaging/LiHM_NMF/atlas/mask_thr0p5_wmparc.2_cc.nii.gz';
maskNii = load_untouch_nii(maskFile);

gNb = createPrepData('volumetric', maskNii.img, 1);

% for surface data
% surfL = 'C:\Users\LiHon\Google Drive\Code\Inhouse\ongoing\Code_mvNMF_l21_ard_v3_release\Release\lib\freesurfer\subjects\fsaverage5\surf\lh.pial';
% surfR = 'C:\Users\LiHon\Google Drive\Code\Inhouse\ongoing\Code_mvNMF_l21_ard_v3_release\Release\lib\freesurfer\subjects\fsaverage5\surf\rh.pial';
% surfML = 'C:\Users\LiHon\Google Drive\Code\Inhouse\ongoing\Code_mvNMF_l21_ard_v3_release\Release\lib\freesurfer\subjects\fsaverage5\label\lh.Medial_wall.label';
% surfMR = 'C:\Users\LiHon\Google Drive\Code\Inhouse\ongoing\Code_mvNMF_l21_ard_v3_release\Release\lib\freesurfer\subjects\fsaverage5\label\rh.Medial_wall.label';
% 
% [surfStru, surfMask] = getFsSurf(surfL, surfR, surfML, surfMR);

%gNb = createPrepData('surface', surfStru, 1, surfMask);

% save gNb into file for later use
prepDataName = 'new_CreatePrepData.mat';
save(prepDataName, 'gNb');

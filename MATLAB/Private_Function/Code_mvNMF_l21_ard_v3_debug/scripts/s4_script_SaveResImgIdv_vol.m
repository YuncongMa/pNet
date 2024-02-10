clc;
clear;

%addpath(genpath('C:\Users\LiHon\Google Drive\Code\Inhouse\ongoing\Code_mvNMF_l21_ard_v3_release\Release'));
%addpath(genpath('D:\Google_drive\Code\Inhouse\ongoing\Code_mvNMF_l21_ard_v3_release\Release'));

resFileName = '/cbica/home/zhouz/projects/istaging/LiHM_NMF/BLSA_test/res_part1/blsa_vol_sbj20_comp17_alphaS21_40_alphaL10_vxInfo0_ard1_eta1/final_UV.mat';
maskName = '/cbica/home/zhouz/projects/istaging/LiHM_NMF/atlas/MNI-maxprob-thr50-2mm-mask.nii.gz';
outDir = '/cbica/home/zhouz/projects/istaging/LiHM_NMF/BLSA_test/res_part1/blsa_vol_sbj20_comp17_alphaS21_40_alphaL10_vxInfo0_ard1_eta1/fig_idv';
saveFig = 1;
refNiiName = '/cbica/home/zhouz/projects/istaging/LiHM_NMF/atlas/MNI-maxprob-thr50-2mm-mask.nii.gz';
sbjId = 1;

func_saveVolRes2Nii_idv(resFileName,maskName,sbjId,outDir,saveFig,refNiiName);

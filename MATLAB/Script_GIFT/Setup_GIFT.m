% This script is setup GIFT for obtaining group ICA results
% It uses a customized function called gica_cmd_Yuncong other than the
% original version to setup the number of PCA for data dimension reduction


%%%%%% Files and directories
% directory to the gift toolbox
dir_gift='./gift';
addpath(genpath(dir_gift));
% root directory for the fMRI dataset
dir_fMRI_dataset='.';
% find all fMRI scan files that end with '.nii.gz'
list_fMRI_scan=dir(fullfile(dir_fMRI_dataset, ['**',filesep,'*.nii.gz']));
% mask file
file_mask='./Brain_Mask.nii.gz';
% output directory for GIFT
dir_output='./group_ICA/Result';

%%%%%% Parameters for gICA
% Number of PCA components for dimension reduction
N_PCA=200;
% Number of ICA components
N_ICA=100;
% Repetitions to get robust group-ICA results via the ICASSO method
N_ICASSO=50;

%%%%%% prepare files
% create a scan list file
cd(dir_output);
fileID = fopen('Scan_List.txt','w');
for i=1:length(list_fMRI_scan)
    fprintf(fileID,'%s\n',fullfile(list_fMRI_scan(i).folder,list_fMRI_scan(i).name));
end
fclose(fileID);

% Steps need to do in group-ICA
GICA_allSteps = {'all', 'parameter_initialization', 'group_pca', 'calculate_ica', 'back_reconstruct', 'scale_components', 'group_stats', 'resume'};
Setup_ICA.groupICAStep=GICA_allSteps([2,3,4]);

% preprocessing steps
%       1. Remove Mean Per Timepoint
%       2. Remove Mean Per Voxel
%       3. Intensity Normalization
%       4. Variance Normalization


%%%%%% Run group ICA
% parallel computation is enabled to use 8 cores
gica_cmd_Yuncong(Setup_ICA,...
    '--data',fullfile(dir_output,'Scan_List.txt'),...
    '--mask',file_mask,...
    '--sess','1',...
    '--output',Output_Path,...
    '--preproc','1',...
    '--n',num2str(N_ICA),...
    '--algorithm','infomax',...
    '--icasso',num2str(N_ICASSO),...
    '--reductions','2',...
    '--df',num2str(N_PCA),...
    '--parallel','8');


%%%%%% Read the group ICA results and output results for pNet

% output brain mask file into Brain_Mask.mat
NII=load_nii(file_mask);
Brain_Mask = NII.img;
save(fullfile(dir_output, 'Brain_Mask.mat'), 'Brain_Mask');
% output group ICA results into gFN_gICA.mat
FN=load(fullfile(dir_output, '_icat_.mat'),'sc');
save(fullfile(dir_output, ['gFN_',num2str(N_ICA),'_gICA.mat']), 'FN');




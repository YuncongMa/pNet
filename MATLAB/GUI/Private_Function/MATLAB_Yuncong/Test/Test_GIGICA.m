% Yuncong Ma, 1/29/2024

%% Load data

%Mask=fLoad_MATLAB_Single_Variable('/Volumes/Scratch_0/pNet/Example/UKBB_Volume/Data/Brain_Mask.mat');

temp=load_nii('/Users/yuncongma/Desktop/UKBiobank_BrainImaging_GroupMeanTemplates/rfMRI_ICA_d25.nii.gz');
gFN = temp.img;
Mask=sum(abs(gFN),4)>0;
Brain_Mask=Mask;
save('/Volumes/Scratch_0/pNet/Example/UKBB_Volume/Test_FN17_GIGICA/Data_Input/Brain_Mask.mat','Brain_Mask');
save_nii(make_nii(double(Brain_Mask),[2,2,2]),'/Volumes/Scratch_0/pNet/Example/UKBB_Volume/Test_FN17_GIGICA/Data_Input/Brain_Mask.nii')

temp = load_nii('/Volumes/Scratch_0/pNet/Example/UKBB_Volume/Data/1054816/filtered_func_data_clean.nii.gz');
Data=temp.img;
Data=fApply_Mask(Mask, Data, -1)';
%Data = zscore(Data, 1);



FID=fopen('/Users/yuncongma/Desktop/UKBiobank_BrainImaging_GroupMeanTemplates/rfMRI_GoodComponents_d25_v1.txt');
List_goodICA=textscan(FID, '%d');
List_goodICA=List_goodICA{1};
fclose(FID);

gFN=gFN(:,:,:,List_goodICA);

dir_pNet_result='/Volumes/Scratch_0/pNet/Example/UKBB_Volume/Test_FN17_GIGICA';

%% Prepare

[pFN, TC] = pFN_GIGICA(zscore(Data,[],1), fApply_Mask(Mask, gFN, -1));

pFN = fInverse_Mask(Mask, pFN, -1);

%% output
temp.img=pFN;
save_nii(temp, '/Users/yuncongma/Desktop/GIG-ICA toolbox20220506/pFN.nii.gz')

%% spatial similarity

SS = corr(fApply_Mask(Mask, gFN, -1), fApply_Mask(Mask, pFN, -1));

mean(diag(SS))

Delta_SS = diag(SS) - max(SS-diag(diag(SS)),[],2);

% FN=pFN;
% save(fullfile(dir_pNet_result,'Personalized_FN','1054816','FN.mat'),'FN');

%% Batch processing

dir_pNet_result='/Volumes/Scratch_0/pNet/Example/UKBB_Volume/Test_FN17_GIGICA';

FID=fopen(fullfile(dir_pNet_result,'Data_Input/Subject_Folder.txt'));
Subject_Folder=textscan(FID,'%s\n');
Subject_Folder=Subject_Folder{1};
fclose(FID);

for i=1:length(Subject_Folder)
    fprintf('%s\n',Subject_Folder{i})
    % get pFNs
    temp = load_nii(fullfile('/Volumes/Scratch_0/pNet/Example/UKBB_Volume/Data',Subject_Folder{i},'/filtered_func_data_clean.nii.gz'));
    Data=temp.img;
    Data=fApply_Mask(Mask, Data, -1)';
    [pFN, TC] = GIGICA(zscore(Data,[],1), fApply_Mask(Mask, gFN, -1));
    pFN = fInverse_Mask(Mask, pFN, -1);
    FN=pFN;
    save(fullfile(dir_pNet_result,'Personalized_FN',Subject_Folder{i},'FN.mat'),'FN');
    %continue
    % visualize
    File_FN=fullfile(dir_pNet_result,'Personalized_FN',Subject_Folder{i},'FN.mat');
    File_Brain_Mask=fullfile(dir_pNet_result,'Data_Input/Brain_Mask.mat');
    File_Overlay_Image=fullfile(dir_pNet_result,'Data_Input/Overlay_Image.mat');
    Dir_Figure=fullfile(dir_pNet_result,'Personalized_FN',Subject_Folder{i});
    fVisualize_FN_Volume(File_FN,File_Brain_Mask,File_Overlay_Image,Dir_Figure,{'Color_Range_Style','GIG-ICA'});
    return
end

%% visualize gFN
File_Brain_Mask=fullfile(dir_pNet_result,'Data_Input/Brain_Mask.mat');
File_Overlay_Image=fullfile(dir_pNet_result,'Data_Input/Overlay_Image.mat');
File_FN=fullfile(dir_pNet_result,'Group_FN','FN.mat');
Dir_Figure=fullfile(dir_pNet_result,'Group_FN');
fVisualize_FN_Volume(File_FN,File_Brain_Mask,File_Overlay_Image,Dir_Figure,{'Color_Range_Style','GIG-ICA'},{'Output_Setting',1});

%% visualize pFN
dir_pNet_result='/Volumes/Scratch_0/pNet/Example/UKBB_Volume/Test_FN17_GIGICA';

FID=fopen(fullfile(dir_pNet_result,'Data_Input/Subject_Folder.txt'));
Subject_Folder=textscan(FID,'%s\n');
Subject_Folder=Subject_Folder{1};
fclose(FID);

i=1;
File_FN=fullfile(dir_pNet_result,'Personalized_FN',Subject_Folder{i},'FN.mat');
File_Brain_Mask=fullfile(dir_pNet_result,'Data_Input/Brain_Mask.mat');
File_Overlay_Image=fullfile(dir_pNet_result,'Data_Input/Overlay_Image.mat');
Dir_Figure=fullfile(dir_pNet_result,'Personalized_FN',Subject_Folder{i});
fVisualize_FN_Volume(File_FN,File_Brain_Mask,File_Overlay_Image,Dir_Figure, ...
    {'Color_Range_Style','GIG-ICA'},{'Setting_Folder',fullfile(dir_pNet_result,'Group_FN','Figure_Setting')});

%% Prepare groupFN
temp=load_nii('/Users/yuncongma/Desktop/UKBiobank_BrainImaging_GroupMeanTemplates/rfMRI_ICA_d25.nii.gz');
gFN = temp.img;
FN=gFN;
save(fullfile(dir_pNet_result,'Group_FN/FN.mat'),'FN');



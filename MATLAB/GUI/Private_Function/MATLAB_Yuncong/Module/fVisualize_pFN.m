function [Flag,Message]=fVisualize_pFN(App_Dir,Work_Dir)
% Yuncong Ma, 2/1/2024
% Make visualization for personalized FN
% [Flag,Message]=fVisualize_pFN(Work_Dir)

Flag=0;
Message='';

Setting.Load_Data=fLoad_MATLAB_Single_Variable(fullfile(Work_Dir,'Data_Input','Setting.mat'));
Setting.Compute_FN=fLoad_MATLAB_Single_Variable(fullfile(Work_Dir,'FN_Computation','Setting.mat'));

% set the color range style
switch Setting.FN_Computation.Method
    case 'SR-NMF'
        Color_Range_Style='SR-NMF';
    case 'GIG-ICA'
        Color_Range_Style='GIG-ICA';
    otherwise
        Flag=1;
        Message=['Unknow pFN method: ',Setting.FN_Computation.Method];
        return
end


FID=fopen(fullfile(Work_Dir,'Data_Input','Subject_Folder.txt'));
Subject_Folder=textscan(FID,'%s\n');
Subject_Folder=Subject_Folder{1};
Subject_Folder_Unique=fUnique_Cell_String(Subject_Folder);
fclose(FID);

% Start parallel
if Setting.Compute_FN.Parallel.Flag
    CPU=gcp('nocreate');
    if isempty(CPU)
        parpool(Setting.Compute_FN.Parallel.N_Thread);
    elseif CPU.NumWorkers~=Setting.Compute_FN.Parallel.N_Thread
        delete(CPU);
        parpool(Setting.Compute_FN.Parallel.N_Thread);
    end
end

N_Subject_Folder_Unique=length(Subject_Folder_Unique);
switch Setting.Load_Data.Data_Format
    case 'HCP Surface (*.cifti, *.mat)'
        File_Brain_Surface=fullfile(Work_Dir,'Data_Input','Brain_Surface.mat');
        File_FN=cell(N_Subject_Folder_Unique,1);
        for i=1:N_Subject_Folder_Unique
            File_FN{i}=fullfile(Work_Dir,'Personalized_FN',Subject_Folder_Unique{i},'FN.mat');
            if ~exist(File_FN{i},'file')
                Flag=1;
                Message=['Cannot find FN.mat in folder: ',File_FN{i}];
                return
            end
        end
        if Setting.Compute_FN.Parallel.Flag==0
            for i=1:N_Subject_Folder_Unique
                fVisualize_FN_HCP_Surface(File_FN{i},File_Brain_Surface,fullfile(Work_Dir,'Personalized_FN',Subject_Folder{i}), ...
                    {'Color_Range_Style',Color_Range_Style});
            end
        else
            parfor i=1:N_Subject_Folder_Unique
                fVisualize_FN_HCP_Surface(File_FN{i},File_Brain_Surface,fullfile(Work_Dir,'Personalized_FN',Subject_Folder{i}), ...
                    {'Color_Range_Style',Color_Range_Style});
            end
        end

    case {'MGH Surface (*.mgh)','Surface (*.mgz)'}
        File_Brain_Surface=fullfile(Work_Dir,'Data_Input','Brain_Surface.mat');
        File_FN=cell(N_Subject_Folder_Unique,1);
        for i=1:N_Subject_Folder_Unique
            File_FN{i}=fullfile(Work_Dir,'Personalized_FN',Subject_Folder_Unique{i},'FN.mat');
            if ~exist(File_FN{i},'file')
                Flag=1;
                Message=['Cannot find FN.mat in folder: ',File_FN{i}];
                return
            end
        end
        if Setting.Compute_FN.Parallel.Flag==0
            for i=1:N_Subject_Folder_Unique
                fVisualize_FN_FreeSurfer(File_FN{i},File_Brain_Surface,fullfile(Work_Dir,'Personalized_FN',Subject_Folder{i}), ...
                    {'Color_Range_Style',Color_Range_Style});
            end
        else
            parfor i=1:N_Subject_Folder_Unique
                fVisualize_FN_FreeSurfer(File_FN{i},File_Brain_Surface,fullfile(Work_Dir,'Personalized_FN',Subject_Folder{i}), ...
                    {'Color_Range_Style',Color_Range_Style});
            end
        end

    case 'Volume (*.nii, *.nii.gz, *.mat)'
        Setting=struct('Synchronized_View',1, 'Synchronized_Colorbar',0);
        File_Setting=fullfile(Work_Dir,'Personalized_FN','Setting.json');
        write_json_setting(Setting, File_Setting);

        for i=1:N_Subject_Folder_Unique
            File_FN{i}=fullfile(Work_Dir,'Personalized_FN',Subject_Folder_Unique{i},'FN.mat');
            if ~exist(File_FN{i},'file')
                Flag=1;
                Message=['Cannot find FN.mat in folder: ',File_FN{i}];
                return
            end
        end

        File_Brain_Mask=fullfile(Work_Dir,'Data_Input','Brain_Mask.mat');
        File_Overlay_Image=fullfile(Work_Dir,'Data_Input','Overlay_Image.mat');
        
        if Setting.Compute_FN.Parallel.Flag==0
            for i=1:N_Subject_Folder_Unique
                fVisualize_FN_Volume(File_FN{i},File_Brain_Mask,File_Overlay_Image,fullfile(Work_Dir,'Personalized_FN',Subject_Folder{i}), ...
                    {'Color_Range_Style',Color_Range_Style}, {'Setting_Folder',fullfile(Work_Dir,'Group_FN','Figure_Setting')});
            end
        else
            parfor i=1:N_Subject_Folder_Unique
                fVisualize_FN_Volume(File_FN{i},File_Brain_Mask,File_Overlay_Image,fullfile(Work_Dir,'Personalized_FN',Subject_Folder{i}), ...
                    {'Color_Range_Style',Color_Range_Style}, {'Setting_Folder',fullfile(Work_Dir,'Group_FN','Figure_Setting')});
            end
        end

    otherwise
        Flag=1;
        Message=['Unknown data format : ',Setting.Load_Data.Data_Format];
end


end



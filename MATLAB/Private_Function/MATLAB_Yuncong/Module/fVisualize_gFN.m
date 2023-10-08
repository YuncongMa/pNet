function [Flag,Message]=fVisualize_gFN(App_Dir,Work_Dir)
% Yuncong Ma, 2/27/2023
% Make visualization for group-level FN
% [Flag,Message]=fVisualize_gFN(Work_Dir)

Flag=0;
Message='';

Setting.Load_Data=fLoad_MATLAB_Single_Variable(fullfile(Work_Dir,'Load_Data','Setting.mat'));

switch Setting.Load_Data.Data_Format
    case 'HCP Surface (*.cifti, *.mat)'
        File_Brain_Surface=fullfile(Work_Dir,'Load_Data','Brain_Surface.mat');
        if ~exist(fullfile(Work_Dir,'Group_FN','FN.mat'),'file')
            Flag=1;
            Message='Cannot find the FN.mat in folder Group_FN';
        end
        fVisualize_FN_HCP_Surface(fullfile(Work_Dir,'Group_FN','FN.mat'),File_Brain_Surface,fullfile(Work_Dir,'Group_FN'));

    case {'MGH Surface (*.mgh)','Surface (*.mgz)'}
        File_Brain_Surface=fullfile(Work_Dir,'Load_Data','Brain_Surface.mat');
        if ~exist(fullfile(Work_Dir,'Group_FN','FN.mat'),'file')
            Flag=1;
            Message='Cannot find the FN.mat in folder Group_FN';
        end
        fVisualize_FN_FreeSurfer(fullfile(Work_Dir,'Group_FN','FN.mat'),File_Brain_Surface,fullfile(Work_Dir,'Group_FN'));
        
    case 'Volume (*.nii, *.nii.gz, *.mat)'
        File_Brain_Mask=fullfile(Work_Dir,'Load_Data','Brain_Mask.mat');
        File_Overlay_Image=fullfile(Work_Dir,'Load_Data','Overlay_Image.mat');
        if ~exist(fullfile(Work_Dir,'Group_FN','FN.mat'),'file')
            Flag=1;
            Message='Cannot find the FN.mat in folder Group_FN';
        end
        fVisualize_FN_Volume(fullfile(Work_Dir,'Group_FN','FN.mat'),File_Brain_Mask,File_Overlay_Image,fullfile(Work_Dir,'Group_FN'));
        
    otherwise
        Flag=1;
        Message=['Unknown data format : ',Setting.Load_Data.Data_Format];
end

end



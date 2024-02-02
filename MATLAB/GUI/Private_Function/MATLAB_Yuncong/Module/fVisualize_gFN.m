function [Flag,Message]=fVisualize_gFN(App_Dir,Work_Dir)
% Yuncong Ma, 2/1/2024
% Make visualization for group-level FN
% [Flag,Message]=fVisualize_gFN(Work_Dir)

Flag=0;
Message='';

Setting.Load_Data=fLoad_MATLAB_Single_Variable(fullfile(Work_Dir,'Data_Input','Setting.mat'));
Setting.FN_Computation=fLoad_MATLAB_Single_Variable(fullfile(Work_Dir,'FN_Computation','Setting.mat'));

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

switch Setting.Load_Data.Data_Format
    case 'HCP Surface (*.cifti, *.mat)'
        File_Brain_Surface=fullfile(Work_Dir,'Data_Input','Brain_Surface.mat');
        if ~exist(fullfile(Work_Dir,'Group_FN','FN.mat'),'file')
            Flag=1;
            Message='Cannot find the FN.mat in folder Group_FN';
        end
        fVisualize_FN_HCP_Surface(fullfile(Work_Dir,'Group_FN','FN.mat'),File_Brain_Surface,fullfile(Work_Dir,'Group_FN'), ...
            {'Color_Range_Style',Color_Range_Style});

    case {'MGH Surface (*.mgh)','Surface (*.mgz)'}
        File_Brain_Surface=fullfile(Work_Dir,'Data_Input','Brain_Surface.mat');
        if ~exist(fullfile(Work_Dir,'Group_FN','FN.mat'),'file')
            Flag=1;
            Message='Cannot find the FN.mat in folder Group_FN';
        end
        fVisualize_FN_FreeSurfer(fullfile(Work_Dir,'Group_FN','FN.mat'),File_Brain_Surface,fullfile(Work_Dir,'Group_FN'), ...
            {'Color_Range_Style',Color_Range_Style});
        
    case 'Volume (*.nii, *.nii.gz, *.mat)'
        File_Brain_Mask=fullfile(Work_Dir,'Data_Input','Brain_Mask.mat');
        File_Overlay_Image=fullfile(Work_Dir,'Data_Input','Overlay_Image.mat');
        if ~exist(fullfile(Work_Dir,'Group_FN','FN.mat'),'file')
            Flag=1;
            Message='Cannot find the FN.mat in folder Group_FN';
        end
        
        fVisualize_FN_Volume(fullfile(Work_Dir,'Group_FN','FN.mat'),File_Brain_Mask,File_Overlay_Image,fullfile(Work_Dir,'Group_FN'), ...
            {'Color_Range_Style',Color_Range_Style}, {'Output_Setting',1});
        
    otherwise
        Flag=1;
        Message=['Unknown data format : ',Setting.Load_Data.Data_Format];
end

end



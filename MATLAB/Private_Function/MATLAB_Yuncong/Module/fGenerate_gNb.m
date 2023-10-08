function [Flag,Message]=fGenerate_gNb(Data_Format,Dir)
% Yuncong Ma, 3/2/2023
% Use files in Load_Data to generate gNb file in Compute_FN
% [Flag,Message]=fGenerate_gNb(Data_Format,Dir)

Flag=0;
Message='';

nb=1;

switch Data_Format
    case {'HCP Surface (*.cifti, *.mat)','MGH Surface (*.mgh)','MGZ Surface (*.mgz)'}
        Brain_Surface=fLoad_MATLAB_Single_Variable(fullfile(Dir,'Load_Data','Brain_Surface.mat'));
        surfStru.faces_l = Brain_Surface.Shape.L.faces;
        surfStru.faces_r = Brain_Surface.Shape.R.faces;
        surfStru.vx_l = Brain_Surface.Shape.L.vertices;
        surfStru.vx_r = Brain_Surface.Shape.R.vertices;
        surfMask.l = Brain_Surface.MW.L==0;
        surfMask.r = Brain_Surface.MW.R==0;
        
        gNb = constructW_surf(surfStru,nb,surfMask);
        save(fullfile(Dir,'Compute_FN','gNb.mat'),'gNb');

    case 'Volume (*.nii, *.nii.gz, *.mat)'
        Brain_Mask=fLoad_MATLAB_Single_Variable(fullfile(Dir,'Load_Data','Brain_Mask.mat'));
        gNb = constructW_vol(Brain_Mask,nb);
        save(fullfile(Dir,'Compute_FN','gNb.mat'),'gNb');
    otherwise
        Flag=1;
        Message=['Unsupported data format: ',Data_Format];

end
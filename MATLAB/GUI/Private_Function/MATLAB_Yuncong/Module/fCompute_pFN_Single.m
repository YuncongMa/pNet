function [Flag,Message]=fCompute_pFN_Single(App_Dir,Work_Dir,Subject_Folder)
% Yuncong Ma, 2/1/2024
% Compute Personalized FN in each subject folder
% [Flag,Message]=fCompute_pFN_Single(App_Dir,Work_Dir,Subject_Folder)
% It uses either SR-NMF or GIG-ICA to compute pFNs

Flag=0;
Message='';

FN_Computation=fLoad_MATLAB_Single_Variable(fullfile(Work_Dir,'FN_Computation','Setting.mat'));

FN_Model=FN_Computation.Method;

switch FN_Model
    case 'SR-NMF'
        pFN_SR_NMF(App_Dir,Work_Dir,Subject_Folder);
    case 'GIG-ICA'
        pFN_GIG_ICA(App_Dir,Work_Dir,Subject_Folder);
    otherwise
        Flag=1;
        Message=['Unknown pFN model ', FN_Model];
end

end






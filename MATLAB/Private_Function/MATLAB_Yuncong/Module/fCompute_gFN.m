function [Flag,Message]=fCompute_gFN(App_Dir,Work_Dir)
% Yuncong Ma, 2/1/2024
% Perform group-level FN computation
% [Flag,Message]=fCompute_gFN(App_Dir,Work_Dir)

Flag=0;
Message='';

FN_Computation=fLoad_MATLAB_Single_Variable(fullfile(Work_Dir,'FN_Computation','Setting.mat'));

FN_Model=FN_Computation.Method;

switch FN_Model
    case 'SR-NMF'
        gFN_SR_NMF(App_Dir,Work_Dir);
    otherwise
        Flag=1;
        Message=['Unknown pFN model ', FN_Model, ' for obtaining group-level FNs'];
end

end






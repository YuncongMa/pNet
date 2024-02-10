function [Flag, Message]=check_gFN(Work_Dir)
% Yuncong Ma, 2/1/2024
% 

% load gFN
[gFN,Flag,Message]=fReformat_FN(App_Dir,Work_Dir,fullfile(Work_Dir,'Group_FN','FN.mat'));
if Flag
    return;
end

FN_Computation=fLoad_MATLAB_Single_Variable(fullfile(Work_Dir,'FN_Computation','Setting.mat'));
FN_Model=FN_Computation.Method;


switch FN_Model
    case 'SR-NMF'
        if sum(gFN(:)<0)>0
            Flag=1;
            Message='SR-NMF requires gFNs to be non-negative';
        end
    case 'GIG-ICA'
        K=size(gFN,2);
        if sum(gFN<0,1)+sum(gFN>0,1)<2*K
            Flag=1;
            Message='GIG-ICA requires inputs have both positive and negative values in each group-level FN';
        end
    otherwise
        Flag=1;
        Message=['Unknown pFN model ', FN_Model];
end

end
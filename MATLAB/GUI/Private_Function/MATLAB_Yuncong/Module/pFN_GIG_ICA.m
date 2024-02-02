function [Flag,Message]=pFN_GIG_ICA(App_Dir,Work_Dir,Subject_Folder)
% Yuncong Ma, 2/1/204
% Compute Personalized FN in each subject folder using GIG-ICA
% [Flag,Message]=pFN_GIG_ICA(App_Dir,Work_Dir,Subject_Folder)

Flag=0;
Message='';

outDir=fullfile(Work_Dir,'Personalized_FN',Subject_Folder);

% Parameter
Scan_File=fullfile(outDir,'Scan_List.txt');
FID=fopen(Scan_File,'r');
Scan_List=textscan(FID,'%s\n');
Scan_List=Scan_List{1};
fclose(FID);


FN_Computation=fLoad_MATLAB_Single_Variable(fullfile(Work_Dir,'FN_Computation','Setting.mat'));
K=FN_Computation.K;

Options.threshold_eign=FN_Computation.PersonalizedFN.threshold_eign;
Options.maxIter=FN_Computation.PersonalizedFN.maxIter;
Options.a=FN_Computation.PersonalizedFN.a;
Options.EGv=FN_Computation.PersonalizedFN.EGv;
Options.ErChuPai=FN_Computation.PersonalizedFN.ErChuPai;
Options.ftol=FN_Computation.PersonalizedFN.ftol;
Options.error=FN_Computation.PersonalizedFN.error;
Options.Nembda=FN_Computation.PersonalizedFN.Nembda;

% Log file
FID=fopen(fullfile(outDir,'Log.txt'),'w');
fprintf(FID,'pFN_GIG_ICA\n Start at %s\n',char(datetime('now')));

% Load data
[Data,Flag,Message]=fLoad_Scan(App_Dir,Work_Dir,Scan_List);
% Concatenate Data on time dimension after normalization
temp=[];
for i=1:length(Data)
    temp(end+1:end+size(Data{i},1),:)=zscore(Data{i},[],1);
end
Data=temp;
clear temp

if Flag
    return
end
[gFN,Flag,Message]=fReformat_FN(App_Dir,Work_Dir,fullfile(Work_Dir,'Group_FN','FN.mat'));
if Flag
    return
end

% Do GIG-ICA
[FN, TC] = GIGICA(Data, gFN, Options);

% finalize the results
[FN,Flag,Message]=fReformat_FN(App_Dir,Work_Dir,FN);
if Flag
    return
end
save(fullfile(outDir,'FN.mat'),'FN');

save(fullfile(outDir,'TC.mat'),'TC');

fprintf(FID,'\nDone');
fclose(FID);


end
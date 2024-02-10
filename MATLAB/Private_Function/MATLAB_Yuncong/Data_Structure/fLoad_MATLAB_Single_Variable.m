function [Data,Flag,Message] = fLoad_MATLAB_Single_Variable(File)
% Yuncong Ma, 3/17/2023
% Load matlab file with a single variable
% [Data,Flag,Message] = fLoad_MATLAB_Single_Variable(File)

Data=[];
Flag=0;
Message='';

if ~exist(File,"file")
    Flag=1;
    Message=['Cannot find the file ',File];
    Data=[];
    return;
end
[~,~,EXT]=fileparts(File);
if ~strcmp(EXT,'.mat')
    Flag=1;
    Message='This file needs to be in MATLAB data format (.mat)';
    return
end
temp=load(File);
FN=fieldnames(temp);
if length(FN)>1
    Flag=1;
    Message='This matlab file contains more than a single variable';
else
    Data=temp.(FN{1});
end

end
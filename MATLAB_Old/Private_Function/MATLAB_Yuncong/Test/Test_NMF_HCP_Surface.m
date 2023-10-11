% Yuncong Ma, 2/22/2023
% Test computation parts in NMF APP

%% Use HCP data in mat format to test
App_Dir='/Volumes/Data/Users/yuncongma/Documents/Document/fMRI/Myworks/MATLAB_APP/NMF';
Work_Dir='/Volumes/Data/Users/yuncongma/Documents/Document/fMRI/Myworks/NMF/Result/HCP/Surface/Test_17';

[Flag,Message]=fCompute_gFN(App_Dir,Work_Dir);

[Flag,Message]=fVisualize_gFN(App_Dir,Work_Dir);

[Flag,Message]=fCompute_iFN(App_Dir,Work_Dir);

[Flag,Message]=fVisualize_iFN(App_Dir,Work_Dir);


% [FN,Flag,Message]=fReformat_FN(App_Dir,Work_Dir,zeros(59412,10));


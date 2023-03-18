clc;
clear;

%addpath(genpath('C:\Users\LiHon\Google Drive\Code\Inhouse\ongoing\Code_mvNMF_l21_ard_v3_release\Release'));
%addpath(genpath('D:\Google_drive\Code\Inhouse\ongoing\Code_mvNMF_l21_ard_v3_release\Release'));


dataset = 'All';
%K = [175,200];
K = [300,400];

for i=1:length(K)
    k=K(i);
    candidateLstFile = strcat('/cbica/home/zhouz/projects/istaging/LiHM_NMF/result/',dataset,'/',num2str(k),'_network/init_file_network.txt');
    outDir = strcat('/cbica/home/zhouz/projects/istaging/LiHM_NMF/result/',dataset,'/',num2str(k),'_network/robust_init');
    initV = selRobustInit(candidateLstFile,k,outDir);
end
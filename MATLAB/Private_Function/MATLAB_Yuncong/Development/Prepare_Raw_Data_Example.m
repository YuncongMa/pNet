% Yuncong Ma, 2/20/2023
% Prepare the short-time HCP surface data in mat format
% Prepare multi site short-time volume data in mat format
% ~10 subject data for each format
% 
%% Prepare HCP surface data in mat format
% 20 scans for HCP in surface format

N_Subject=10;
N_Scan_Per_Subject=2;
N_t=400;

Raw_Data=fLoad_SSD_Structure('/Volumes/Scratch_0/Scratch/Raw_Data_HCP_FS_REST12_FIX');

Path_Test_Data='/Volumes/Data/Users/yuncongma/Documents/Document/fMRI/Myworks/NMF/Result/HCP/Surface/APP_Data';
fMake_Folder(Path_Test_Data);

for scan=1:N_Subject*N_Scan_Per_Subject
    Data=fSSD_Structure_Operation(Raw_Data,'Get_Value',scan,'Image_Raw');
    Data=Data(:,1:N_t);
    Data=fTemporal_Smoothing(ones(59412,1),Data,0.72,0.01,.08,{'Parallel',0},{'Order',6});
    Subject_ID=num2str(fSSD_Structure_Operation(Raw_Data,'Get_Value',scan,'Subject'));
    REST=num2str(fSSD_Structure_Operation(Raw_Data,'Get_Value',scan,'REST'));
    Type=fSSD_Structure_Operation(Raw_Data,'Get_Value',scan,'Type');
    fMake_Folder(fullfile(Path_Test_Data,Subject_ID,REST,Type));
    save(fullfile(Path_Test_Data,Subject_ID,REST,Type,'Image.mat'),'Data');
end

%% Zhen's multi site volume data in mat format









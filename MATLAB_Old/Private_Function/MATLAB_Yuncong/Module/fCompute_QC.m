function [Flag,Message]=fCompute_QC(App_Dir,Work_Dir)
% Yuncong Ma, 2/21/2023
% Peform quality control of individualized FN computation
% [Flag,Message]=fCompute_QC(App_Dir,Work_Dir)

Flag=0;
Message='';


Setting.Load_Data=fLoad_MATLAB_Single_Variable(fullfile(Work_Dir,'Load_Data','Setting.mat'));
Setting.Compute_FN=fLoad_MATLAB_Single_Variable(fullfile(Work_Dir,'Compute_FN','Setting.mat'));

FID=fopen(fullfile(Work_Dir,'Load_Data','Scan_List.txt'));
Scan_File=textscan(FID,'%s\n');
Scan_File=Scan_File{1};
fclose(FID);

FID=fopen(fullfile(Work_Dir,'Load_Data','Subject_Folder.txt'));
Subject_Folder=textscan(FID,'%s\n');
Subject_Folder=Subject_Folder{1};
fclose(FID);

Group_FN=fLoad_MATLAB_Single_Variable(fullfile(Work_Dir,'Group_FN','init.mat'));

K=Setting.Compute_FN.K;
N_Scan=length(Subject_Folder);

Result.Spatial_Correspondence=zeros(K,N_Scan);
Result.Spatial_Correspondence_Control=zeros(K,N_Scan);
Result.Spatial_Correspondence_Control_Index=zeros(K,N_Scan);
for i=1:length(Subject_Folder)
    
    temp=dir(fullfile(Work_Dir,'Individualized_FN',Subject_Folder{i},'**','final_UV.mat'));
    if isempty(temp)
        Flag=1;
        Message=['Cannot find the individualized FN in folder: ',fullfile(Work_Dir,'Individualized_FN',Subject_Folder{i})];
        return
    end

    load(fullfile(temp(1).folder,temp(1).name),'V');

    temp=corr(Group_FN,V{1});
    Result.Spatial_Correspondence(:,i)=diag(temp);
    temp=temp-diag(diag(temp));
    [Result.Spatial_Correspondence_Control(:,i),Result.Spatial_Correspondence_Index(:,i)]=max(temp,[],2);
end


QC.Spatial_Correspondence=Result.Spatial_Correspondence;
QC.Spatial_Correspondence_Control=Result.Spatial_Correspondence_Control;
QC.Delta_Similarity=min(Result.Spatial_Correspondence-Result.Spatial_Correspondence_Control,[],1);

[temp,ps]=min(Result.Spatial_Correspondence-Result.Spatial_Correspondence_Control,[],1);
temp=sortrows([temp;1:N_Scan;ps]');
QC.Individual(1:N_Scan)=struct('Miss_Match',[]);

FID=fopen(fullfile(Work_Dir,'Quality_Control','Failed_Scan.txt'),'w');
for i=1:sum(temp(:,1)<0)
    fprintf(FID,'%s ',Result.Individualized_FN(temp(i,2)).SubjectID);
    ps2=find(Result.Spatial_Correspondence(:,temp(i,2))<=Result.Spatial_Correspondence_Control(:,temp(i,2)));
    for j=1:length(ps2)
        fprintf(FID,'%d-%d ',ps2(j),Result.Spatial_Correspondence_Index(ps2(j),temp(i,2)));
        QC.Individual(temp(i,2)).Miss_Match(j,:)=[ps2(j),Result.Spatial_Correspondence_Index(ps2(j),temp(i,2))];
    end
    fprintf(FID,'\n');
end
fclose(FID);

save(fullfile(Work_Dir,'Quality_Control','QC.mat'),'QC');


end









function [Flag,Message]=fCompute_gFN(App_Dir,Work_Dir)
% Yuncong Ma, 3/2/2023
% Perform group-level FN computation
% [Flag,Message]=fCompute_gFN(App_Dir,Work_Dir)

Flag=0;
Message='';

File_Group_FN=fullfile(Work_Dir,'Group_FN','init.mat');

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

% Reorganzie files by subjects
fMake_Folder(fullfile(Work_Dir,'Group_FN','BootStrapping'));
Scan_Organization=struct('Subject_Folder',[],'Scan_List',{},'N_Scan',0);
Bootstrapping(1:Setting.Compute_FN.GroupFN.BootStrap.Repetition)=struct('File',[],'Out_Dir',[]);

% Generate Scan_List in each bootstrapping folder
Subject_Folder_Unique=fUnique_Cell_String(Subject_Folder);
N_Subject_Folder=length(Subject_Folder_Unique);
for i=1:N_Subject_Folder
    ps=find(strcmp(Subject_Folder,Subject_Folder_Unique{i}));
    Scan_Organization(i).Subject_Folder=Subject_Folder_Unique{i};
    Scan_Organization(i).N_Scan=length(ps);
    Scan_Organization(i).Scan_List=Scan_File(ps);
end
for i=1:Setting.Compute_FN.GroupFN.BootStrap.Repetition
    ps=randperm(N_Subject_Folder,Setting.Compute_FN.GroupFN.BootStrap.File_Selection);

    Bootstrapping(i).Out_Dir=fullfile(Work_Dir,'Group_FN','BootStrapping',num2str(i));
    fMake_Folder(Bootstrapping(i).Out_Dir);
    Bootstrapping(i).File=fullfile(Bootstrapping(i).Out_Dir,'Scan_List.txt');

    FID=fopen(Bootstrapping(i).File,'w');
    for j=1:Setting.Compute_FN.GroupFN.BootStrap.File_Selection
        fprintf(FID,'%s\n',Scan_Organization(ps(j)).Scan_List{randi(Scan_Organization(ps(j)).N_Scan)});
    end
    fclose(FID);
end

% Start parallel
if Setting.Compute_FN.Parallel.Flag
    CPU=gcp('nocreate');
    if isempty(CPU)
        parpool(Setting.Compute_FN.Parallel.N_Thread);
    elseif CPU.NumWorkers~=Setting.Compute_FN.Parallel.N_Thread
        delete(CPU);
        parpool(Setting.Compute_FN.Parallel.N_Thread);
    end
end

% Compute each bootstrapping
switch Setting.Load_Data.Data_Format
    case {'HCP Surface (*.cifti, *.mat)','MGH Surface (*.mgh)','MGZ Surface (*.mgz)','Volume (*.nii, *.nii.gz, *.mat)'}
        if Setting.Compute_FN.IndividualizedFN.Combine_Flag==0
            if Setting.Compute_FN.Parallel.Flag==0
                for i=1:Setting.Compute_FN.GroupFN.BootStrap.Repetition
                    [Flag,Message]=fCompute_gFN_Bootstrapping(App_Dir,Work_Dir,Setting.Compute_FN,Bootstrapping(i));
                    if Flag==1
%                         Message=['Error in running fCompute_gFN_Bootstrapping for ',num2str(i),'-th bootstrapping'];
                        return
                    end
                end
            else
                Parallel(1:Setting.Compute_FN.GroupFN.BootStrap.Repetition)=struct('Flag',0,'Message','');
                parfor i=1:Setting.Compute_FN.GroupFN.BootStrap.Repetition
                    [Parallel(i).Flag,Parallel(i).Message]=fCompute_gFN_Bootstrapping(App_Dir,Work_Dir,Setting.Compute_FN,Bootstrapping(i));
                end
                for i=1:Setting.Compute_FN.GroupFN.BootStrap.Repetition
                    if Parallel(i).Flag==1
                        Flag=1;
                        Message=Parallel(i).Message;
%                         Message=['Error in running fCompute_gFN_Bootstrapping for ',num2str(i),'-th bootstrapping'];
                        return
                    end
                end
            end
        end

    otherwise
        Flag=1;
        Message=['Unknown data format : ',Setting.Load_Data.Data_Format];
end

[Flag,Message]=fCompute_gFN_Fusion(App_Dir,Work_Dir,Setting.Compute_FN,Bootstrapping);
if Flag==1
    return
end


end






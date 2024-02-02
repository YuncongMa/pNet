function [Flag,Message]=fCompute_pFN(App_Dir,Work_Dir)
% Yuncong Ma, 2/1/2024
% Perform personalized FN computation
% [Flag,Message]=fCompute_pFN(App_Dir,Work_Dir)

Flag=0;
Message='';

Setting.Data_Input=fLoad_MATLAB_Single_Variable(fullfile(Work_Dir,'Data_Input','Setting.mat'));
Setting.FN_Computation=fLoad_MATLAB_Single_Variable(fullfile(Work_Dir,'FN_Computation','Setting.mat'));

FID=fopen(fullfile(Work_Dir,'Data_Input','Scan_List.txt'));
Scan_File=textscan(FID,'%s\n');
Scan_File=Scan_File{1};
fclose(FID);

FID=fopen(fullfile(Work_Dir,'Data_Input','Subject_ID.txt'));
Subject_Folder=textscan(FID,'%s\n');
Subject_Folder=Subject_Folder{1};
fclose(FID);

% Make Scan_List in each subject folder
Scan_Organization=struct('Subject_Folder',[],'Scan_List',{},'N_Scan',0);
Subject_Folder_Unique=fUnique_Cell_String(Subject_Folder);
N_Subject_Folder=length(Subject_Folder_Unique);
for i=1:N_Subject_Folder
    ps=find(strcmp(Subject_Folder,Subject_Folder_Unique{i}));
    Scan_Organization(i).Subject_Folder=Subject_Folder_Unique{i};
    Scan_Organization(i).N_Scan=length(ps);
    Scan_Organization(i).Scan_List=Scan_File(ps);
    fMake_Folder(fullfile(Work_Dir,'Personalized_FN',Subject_Folder_Unique{i}));
    FID=fopen(fullfile(Work_Dir,'Personalized_FN',Subject_Folder_Unique{i},'Scan_List.txt'),'w');
    for j=1:length(ps)
        fprintf(FID,'%s\n',Scan_Organization(i).Scan_List{j});
    end
    fclose(FID);
end

% Start parallel
if Setting.FN_Computation.Parallel.Flag
    CPU=gcp('nocreate');
    if isempty(CPU)
        parpool('Processes',Setting.FN_Computation.Parallel.N_Thread);
    elseif CPU.NumWorkers~=Setting.FN_Computation.Parallel.N_Thread
        delete(CPU);
        parpool('Processes',Setting.FN_Computation.Parallel.N_Thread);
    end
end

switch Setting.Data_Input.Data_Format
    case {'HCP Surface (*.cifti, *.mat)','MGH Surface (*.mgh)','MGZ Surface (*.mgz)','Volume (*.nii, *.nii.gz, *.mat)'}
        if Setting.FN_Computation.Parallel.Flag==0
            for i=1:N_Subject_Folder
                [Flag,Message]=fCompute_pFN_Single(App_Dir,Work_Dir,Setting.FN_Computation,Subject_Folder{i});
                if Flag==1
                    %                         Message=['Error in running fCompute_pFN_Single in ',Subject_Folder{i}];
                    return
                end
            end
        else
            Parallel(1:N_Subject_Folder)=struct('Flag',0,'Message','');
            parfor i=1:N_Subject_Folder
                [Parallel(i).Flag,Parallel(i).Message]=fCompute_pFN_Single(App_Dir,Work_Dir,Setting.FN_Computation,Subject_Folder{i});
            end
            for i=1:N_Subject_Folder
                if Parallel(i).Flag==1
                    Flag=1;
                    Message=Parallel(i).Message;
                    %                         Message=['Error in running fCompute_pFN_Single in ',Subject_Folder{i}];
                    return
                end
            end
        end


    otherwise
        Flag=1;
        Message=['Unknown data format : ',Setting.Data_Input.Data_Format];
end

end


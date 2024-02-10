function [Data,Flag,Message]=fLoad_Scan(App_Dir,Work_Dir,Scan_List)
% Yuncong Ma, 4/24/2023
% Load fMRI scans to NMF APP
% [Data,Flag,Message]=fLoad_Scan(App_Dir,Work_Dir,Scan_List)
% Works for different formats of fMRI scans
% HCP cifti; FreeSurfer gii, mgz; NIFTI nii, nii.gz; MATLAB .mat
% Support both surface and volume data formats in MATLAB files
% Can concatenate multiple scans on the time dimension
% Return 2D matrix for both surface and volume image data
% Time is the first dimension for NMF
% Scan_List could be a single string or a cell string
% Data is a matrix when Scan_List is a string
% Data is a cell matrix when Scan_List is a cell string

Data=[];
Flag=0;
Message='';

Setting.Load_Data=fLoad_MATLAB_Single_Variable(fullfile(Work_Dir,'Load_Data','Setting.mat'));

if ischar(Scan_List)
    [Data,Flag,Message]=fLoad_Scan_2(Setting, App_Dir, Work_Dir,Scan_List);
elseif iscell(Scan_List)
    Data=cell(length(Scan_List),1);
    for i=1:length(Scan_List)
        [temp,Flag,Message]=fLoad_Scan_2(Setting,App_Dir,Work_Dir,Scan_List{i});
        temp=temp';
        if Flag==1
            return
        else
            if i==1
                Data{1}=temp;
            else
                if size(temp,2)==size(Data{1},2)
                    Data{i}=temp;
                else
                    Flag=1;
                    Message=['Unequal number of vertex/voxel in scan : ',Scan_List{i}];
                    return
                end
            end
        end
    end
else
    Flag=1;
    Message='Only support a single string or a cell string for the input Scan_List';
end


end


function [Data,Flag,Message]=fLoad_Scan_2(Setting,App_Dir,Work_Dir,Scan_File)
Data=[];
Flag=0;
Message='';

if ismac
    OS='macOS';
elseif ispc
    OS='Windows';
elseif isunix
    OS='Linux';
else
    Flag=1;
    Message='Can only work on macOS, Windows and Linux';
    return
end

switch Setting.Load_Data.Data_Format

    case 'HCP Surface (*.cifti, *.mat)'
        if ~exist(Scan_File,'file')
            Flag=1;
            Message=['Cannot find file ',Scan_File];
            return
        end

        [~,~,EXT]=fileparts(Scan_File);
        switch EXT
            case '.nii'
                Struct=cifti_read(Scan_File);
                if isempty(Struct) || ~isfield(Struct,'cdata') || isempty(Struct.cdata)
                    Flag=1;
                    Message='Failed to load the CIFTI file, or the CIFTI file contains an empty matrix';
                    return
                end
                Data=Struct.cdata;
                N_Vertex=59412; % For surface mesh
                if size(Data,1)>N_Vertex
                    Data(N_Vertex+1:end,:)=[];
                elseif size(Data,1)<N_Vertex
                    Flag=1;
                    Message=['For data format HCP Surface, it requires data contains at least ',num2str(N_Vertex),' nodes\nBut the scan ',...
                        Scan_File, ' contains ',nu2mstr(size(Data,1)),' nodes'];
                    return
                end
            case '.mat'
                [Data,Flag,Message]=fLoad_MATLAB_Single_Variable(Scan_File);
                if Flag==1
                    return
                end
                N_Vertex=59412; % For surface mesh
                if size(Data,1)>N_Vertex
                    Data(N_Vertex+1:end,:)=[];
                elseif size(Data,1)<N_Vertex
                    Flag=1;
                    Message=['For data format HCP Surface, it requires data contains at least ',num2str(N_Vertex),' nodes\nBut the scan ',...
                        Scan_File, ' contains ',nu2mstr(size(Data,1)),' nodes'];
                    return
                end
                
            otherwise
                Flag=1;
                Message=['Only support .cifti and .mat for HCP Surface, does not support : ',EXT];
                return
        end

    case 'MGH Surface (*.mgh)'
        Brain_Surface=fLoad_MATLAB_Single_Variable(fullfile(Work_Dir,'Load_Data','Brain_Surface.mat'));
        [~,~,EXT]=fileparts(Scan_File);
        switch EXT
            case '.mgh'
                Scan_File_LR=split(Scan_File,';');
                if ~exist(Scan_File_LR{1},'file')
                    Flag=1;
                    Message=['Cannot find file ',Scan_File_LR{1}];
                    return
                elseif ~exist(Scan_File_LR{2},'file')
                    Flag=1;
                    Message=['Cannot find file ',Scan_File_LR{2}];
                    return
                end
                Data_L=MRIread(Scan_File_LR{1});
                Data_R=MRIread(Scan_File_LR{2});
                Data=[squeeze(Data_L.vol);squeeze(Data_R.vol)];
                Data(Brain_Surface.MW.LR>0,:)=[];

            otherwise
                Flag=1;
                Message=['Only support .mgh for MGH Surface, does not support : ',EXT];
                return
        end

    case 'MGZ Surface (*.mgz)'
        if ~exist(Scan_File_LR{1},'file')
            Flag=1;
            Message=['Cannot find file ',Scan_File_LR{1}];
            return
        end
        [~,~,EXT]=fileparts(Scan_File);
        switch EXT
            case '.mgz'
                temp=MRIread(Scan_File);
                Data=squeeze(temp.vol);

            otherwise
                Flag=1;
                Message=['Only support .mgz for Surface, does not support : ',EXT];
                return
        end

    case 'Volume (*.nii, *.nii.gz, *.mat)'
        if ~exist(Scan_File,'file')
            Flag=1;
            Message=['Cannot find file ',Scan_File];
            return
        end

        Brain_Mask=fLoad_MATLAB_Single_Variable(fullfile(Work_Dir,'Load_Data','Brain_Mask.mat'));

        [~,~,EXT]=fileparts(Scan_File);
        switch EXT
            case '.mat'
                [Data,Flag,Message]=fLoad_MATLAB_Single_Variable(Scan_File);
                if Flag==1
                    return
                end
                if ~isequal(size(Data(:,:,:,1)),size(Brain_Mask))
                    Flag=1;
                    Message=['For volume data type, it requires data size matches to the brain mask file \nBut the scan ',...
                        Scan_File, ' has a different image size ',num2str(size(Data))];
                    return
                else
                    Data=fApply_Mask(Brain_Mask>0,Data,-1);
                end
            case '.gz'
                if contains(Scan_File,'.nii.gz')==0
                    Flag=1;
                    Message=['This is not a .nii.gz file: ',Scan_File];
                    return
                end
                NII=load_untouch_nii(Scan_File);
                Data=NII.img;
                if ~isequal(size(Data(:,:,:,1)),size(Brain_Mask))
                    Flag=1;
                    Message=['For volume data type, it requires data size matches to the brain mask file \nBut the scan ',...
                        Scan_File, ' has a different image size ',num2str(size(Data))];
                    return
                else
                    Data=fApply_Mask(Brain_Mask>0,Data,-1);
                end
            case '.nii'
                NII=load_untouch_nii(Scan_File);
                Data=NII.img;
                if ~isequal(size(Data(:,:,:,1)),size(Brain_Mask))
                    Flag=1;
                    Message=['For volume data type, it requires data size matches to the brain mask file \nBut the scan ',...
                        Scan_File, ' has a different image size ',num2str(size(Data))];
                    return
                else
                    Data=fApply_Mask(Brain_Mask>0,Data,-1);
                end
                
            otherwise
                Flag=1;
                Message=['Only support .nii, .nii.gz and .mat for volume data format, does not support : ',EXT];
                return
        end


    otherwise
        Flag=1;
        Message=['Unsupported data format: ',Setting.Load_Data.Data_Format];
        return
end

if sum(isnan(Data(:)))>0
    Flag=1;
    Message=['NaN is found in scan %s',Scan_File];
    return
end

end





function Result=fAPP_Statistics(File_List,varargin)
% Yuncong Ma, 11/21/2022
% Result=fAPP_Statistics(File_List)
% Perform statistic analysis for NMF APP
% Options support Data_Type, Data_Format, Method, FDR, P_Value, Behavior

Options.Data_Type='Surface';
Options.Data_Format='HCP_Surface';
Options.Method='1 sample ttest';
Options.FDR=0;
Options.P_Value=0.05;
Options.Behavior=[];
Options.File_Brain_Mask='';
Options=fOption('fAPP_Statistics',Options,varargin);
if isempty(Options)
    return;
end

if Options.P_Value==0 || Options.P_Value==1
    error('Error in fAPP_Statistics: pvalue needs to be within 0-1 %f',Options.P_Value);
    return
end


Result=[];

N_File=length(File_List);
% load data
Data=[];
switch Options.Data_Type
    case 'Surface'
        for i=1:N_File
            load(File_List{i},'V');
            V=V{1,1};
            Data(i,:,:)=V'; %[k,N]
        end
    case 'Volume'
        if ~exist(Options.File_Brain_Mask,'file')
            error('Error in fAPP_Statistics: cannot find the brain mask file for volume-based results');
        else
            load(Options.File_Brain_Mask,'Brain_Mask');
        end
        for i=1:N_File
            load(File_List{i},'V');
            V=V{1,1};
            Data(i,:,:)=fApply_Mask(Brain_Mask,V,-1)';
        end
    otherwise
        error('Error in fAPP_Statistics: unknown Data_Type: %s',Options.Data_Type);
end

% Statistics
N_Dim=size(Data,3);
K=size(Data,2);
P_Value=zeros(K,N_Dim);
T_Value=zeros(K,N_Dim);
switch Options.Method
    case '1 sample ttest'
        fig=waitbar(0,'Computing ...');
        for i=1:K
            for j=1:N_Dim
                [~,P_Value(i,j),~,stats]=ttest(squeeze(Data(:,i,j)));
                T_Value(i,j)=stats.tstat;
            end
            if ~isvalid(fig)
                return
            end
            waitbar(i/K,fig,'Computing ...');
        end
        if isvalid(fig)
            close(fig)
        end
        Result.P_Value=P_Value';
        Result.T_Value=T_Value';
    case '2 sample ttest'
        Unique=unique(Options.Behavior);
        Data1=Data(Options.Behavior==Unique(1),:,:);
        Data2=Data(Options.Behavior==Unique(1),:,:);
        clear Data
        fig=waitbar(0,'Computing ...');
        for i=1:K
            for j=1:N_Dim
                [~,P_Value(i,j),~,stats]=ttest2(squeeze(Data1(:,i,j)),squeeze(Data2(:,i,j)));
                T_Value(i,j)=stats.tstat;
            end
            if ~isvalid(fig)
                return
            end
            waitbar(i/K,fig,'Computing ...');
        end
        if isvalid(fig)
            close(fig)
        end
        Result.P_Value=P_Value';
        Result.T_Value=T_Value';
    otherwise
        error('Error in fAPP_Statistics: unknown statistic Method: %s',Options.Method);
end


if Options.FDR
    for k=1:K
        [~,~,P_Value(:,k)]=fdr_bh(P_Value(:,k),Options.P_Value);
    end
end


% output
switch Options.Data_Type
    case 'Volume'
        Result.P_Value=fInverse_Mask(Brain_Mask,Result.P_Value,-1);
        Result.T_Value=fInverse_Mask(Brain_Mask,Result.T_Value,-1);
end

end



















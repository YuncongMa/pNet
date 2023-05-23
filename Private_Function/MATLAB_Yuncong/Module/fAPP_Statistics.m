function Result=fAPP_Statistics(File_List,varargin)
% Yuncong Ma, 5/23/2023
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

if Options.P_Value<=0 || Options.P_Value>=1
    error('Error in fAPP_Statistics: pvalue needs to be within 0-1 %f',Options.P_Value);
    Result=[];
    return
end


Result=[];

N_File=length(File_List);
% load data
Data=[];
switch Options.Data_Type
    case 'Surface'
        fig=waitbar(0,'Loading FN ...');
        for i=1:N_File
            FN=fLoad_MATLAB_Single_Variable(fullfile(File_List{i},'FN.mat'));
            Data(i,:,:)=FN'; %[k,N]
            waitbar(i/N_File,fig,'Loading FN ...');
        end
        if isvalid(fig)
            close(fig)
        end
    case 'Volume'
        if ~exist(Options.File_Brain_Mask,'file')
            error('Error in fAPP_Statistics: cannot find the brain mask file for volume-based results');
        else
            Brain_Mask=fLoad_MATLAB_Single_Variable(Options.File_Brain_Mask,'Brain_Mask');
        end
        fig=waitbar(0,'Loading FN ...');
        for i=1:N_File
            FN=fLoad_MATLAB_Single_Variable(fullfile(File_List{i},'FN.mat'));
            Data(i,:,:)=fApply_Mask(Brain_Mask,FN,-1)';
            waitbar(i/N_File,fig,'Loading FN ...');
        end
        if isvalid(fig)
            close(fig)
        end
    otherwise
        error('Error in fAPP_Statistics: unknown Data_Type: %s',Options.Data_Type);
end
Data(isnan(Data))=0;

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
    case 'Wilcoxon singed rank test'
        fig=waitbar(0,'Computing ...');
        Z_Value=zeros(K,N_Dim);
        for i=1:K
            for j=1:N_Dim
                [P_Value(i,j),~,stats]=signrank(squeeze(Data(:,i,j)));
                Z_Value(i,j)=stats.zval;
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
        Result.Z_Value=Z_Value';
    case '2 sample ttest'
        Unique=unique(Options.Behavior);
        Data1=Data(Options.Behavior==Unique(1),:,:);
        Data2=Data(Options.Behavior==Unique(2),:,:);
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
    case 'Wilcoxon rank sum test'
        Unique=unique(Options.Behavior);
        Data1=Data(Options.Behavior==Unique(1),:,:);
        Data2=Data(Options.Behavior==Unique(2),:,:);
        clear Data
        fig=waitbar(0,'Computing ...');
        Z_Value=zeros(K,N_Dim);
        for i=1:K
            for j=1:N_Dim
                [P_Value(i,j),~,stats]=ranksum(squeeze(Data1(:,i,j)),squeeze(Data2(:,i,j)));
                Z_Value(i,j)=stats.zval;
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
        Result.Z_Value=Z_Value';
    otherwise
        error('Error in fAPP_Statistics: unknown statistic Method: %s',Options.Method);
end

% Remove NaN
% FDR
Result.P_Value(isnan(Result.P_Value))=1;
if Options.FDR
    for k=1:K
        [~,~,Result.P_Value(:,k)]=fdr_bh(Result.P_Value(:,k),Options.P_Value);
    end
end

% Remove NaN
if isfield(Result,'T_Value')
    Result.T_Value(isnan(Result.T_Value))=0;
elseif isfield(Result,'Z_Value')
    Result.Z_Value(isnan(Result.Z_Value))=0;
end


% output
switch Options.Data_Type
    case 'Volume'
        Result.P_Value=fInverse_Mask(Brain_Mask,Result.P_Value,-1);
        if isfield(Result,'T_Value')
            Result.T_Value=fInverse_Mask(Brain_Mask,Result.T_Value,-1);
        end
        if isfield(Result,'Z_Value')
            Result.Z_Value=fInverse_Mask(Brain_Mask,Result.Z_Value,-1);
        end
end

end



















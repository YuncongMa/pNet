% Yuncong Ma, 8/21/2023
% This script is to test the corresponding python file, in order to ensure
% the same result of gFN computation
% Test_gFN_Computation.py

%% Add to path
addpath(genpath('/Users/yuncongma/Documents/Document/fMRI/Myworks/pNet'));
%%
maxNumCompThreads(10);

% load file
Data = fLoad_MATLAB_Single_Variable('/Users/yuncongma/Documents/Document/fMRI/Myworks/pNet/Example/HCP_Surface/Data/100206/1/LR/Image.mat');
gNb = fLoad_MATLAB_Single_Variable('/Users/yuncongma/Documents/Document/fMRI/Myworks/pNet/Example/HCP_Surface/Test_FN17/FN_Computation/gNb.mat');
Setting = fLoad_MATLAB_Single_Variable('/Users/yuncongma/Documents/Document/fMRI/Myworks/pNet/Example/HCP_Surface/Test_FN17/FN_Computation/Setting.mat');

Data=Data';
%%
tic
% Setup parameter
K = 17;
maxIter = 1000;
minIter = 30;
error = 1e-6;
normW = 1;
eta = 0;
Alpha = 2;
Beta = 30;
vxI = 0;
ard = 0;
dataPrecision = 'double';
logFile = 'Log_gFN_NMF_MATLAB.log';

% Internal parameter
alphaS = 0;
alphaL = 0;

% setup log file
% fileID = fopen(logFile, 'w');
fprintf(['\nStart NMF for gFN using MATLAB at ',char(datetime),'\n'])


% [mat_float,mat_eps] = set_data_precision(dataPrecision);
mat_eps = eps(dataPrecision);
if strcmp(dataPrecision,'single')
    Data=single(Data);
else
    Data=double(Data);
end

dim_space = size(Data,2);
dim_time = size(Data,1);

nmVec = zeros(length(gNb),1);
for gni=1:length(gNb)
    nmVec(gni) = length(gNb{gni});
end
nM = median(nmVec);

alphaS = round((Alpha*dim_time)/K);
alphaL = round((Beta*dim_time)/(K*nM));

% normalize data
X = normalize_data(Data,'vp','vmax');

% Initialize U and V
mean_X = sum(X(:))/(dim_time*dim_space);
rng('default')
U = (rand(dim_time,K)+1)*(sqrt(mean_X/K));
V = (rand(dim_space,K)+1)*(sqrt(mean_X/K));

% save('./gFN_init_U.mat','U');
% save('./gFN_init_V.mat','V');

[U,V] = normalize_u_v(U, V, 1, 1);

% construct the spatial affinity graph
[L,W,D] = construct_Laplacian_gNb(gNb, X, vxI, dim_space, alphaL, normW, dataPrecision);


if ard>0
    ard = 1;
    eta = 0.1;
    lambdas = sum(U,1) / dim_time;
    hyperLam = eta * sum(sum(X.^2,1),2) / (dim_time*dim_space*2);
else
    lambdas=0;
    hyperLam=0;
end

oldLogL=inf;

% Multiplicative update of U and V
for i=1:maxIter

    % ===================== update V ========================
    XU = X'*U;  % mnk or pk (p<<mn)
    UU = U'*U;  % mk^2
    VUU = V*UU; % nk^2

    if alphaS>0
        %             % L1 sparsity
        %             F = ones(size(V));
        %             VUU = VUU + 0.5*alphaS * F;

        % scale-invariant sparsity
        tmpNorm2 = sqrt(sum(V.^2,1));
        posTerm = 1 ./ max(repmat(tmpNorm2,dim_space,1),mat_eps);
        tmpNorm1 = sum(V,1);
        negTerm = V .* repmat(tmpNorm1,dim_space,1) ./ max(repmat(tmpNorm2.^3,dim_space,1),mat_eps);

        XU = XU + 0.5*alphaS * negTerm;
        VUU = VUU + 0.5*alphaS * posTerm;
    end

    if alphaL>0
        V = double(V);
        WV = W*V;
        DV = D*V;

        XU = XU + WV;
        VUU = VUU + DV;
    end

    V = V.*(XU./max(VUU,mat_eps));

    prunInd = sum(V~=0)==1;
    if any(prunInd)
        V(:,prunInd) = zeros(dim_space,sum(prunInd));
        U(:,prunInd) = zeros(dim_time,sum(prunInd));
    end

    [U, V] = normalize_u_v(U, V, 1, 1);

    % ===================== update U ========================
    XV = X*V;   % mnk or pk (p<<mn)
    VV = V'*V;  % nk^2
    UVV = U*VV; % mk^2

    if ard>0
        posTerm = 1./max(repmat(lambdas,dim_time,1),mat_eps);
        UVV = UVV + posTerm*hyperLam;
    end

    U = U.*(XV./max(UVV,eps));

    prunInd = sum(U)==0;
    if any(prunInd)
        V(:,prunInd) = zeros(dim_space,sum(prunInd));
        U(:,prunInd) = zeros(dim_time,sum(prunInd));
    end

    if ard>0
        lambdas = sum(U) / dim_time;
    end
    
    % ==== calculate partial objective function value ====
    ardU = 0;
    tmp1 = 0;
    tmp2 = 0;
    tmp3 = 0;
    tmpl21 = V.^2;

    if ard>0
        su=sum(U,1);
        su(su==0)=1;
        ardU=sum(log(su))*dim_time*hyperLam;
    end

    tmpDf=(X-U*V').^2;
    tmp1=sum(tmpDf(:));

    if alphaL>0
        tmp2=V'*L.*V';
    end

    L21=alphaS*sum(sum(sqrt(tmpl21),1)./max(sqrt(sum(tmpl21,1)),mat_eps));
    LDf = tmp1;
    LSl = sum(tmp2(:));

    % Objective function
    LogL = L21 + ardU + LDf + LSl;

    fprintf(['  Iter = ',num2str(i),': LogL:',num2str(LogL),',dataFit:',num2str(LDf)...
        ',spaLap:',num2str(LSl),...
        ',L21:',num2str(L21),',ardU:',num2str(ardU),'\n']);
    
    if i>minIter && abs(oldLogL-LogL)/max(oldLogL,mat_eps)<error
        break;
    end
    oldLogL = LogL;
end


toc

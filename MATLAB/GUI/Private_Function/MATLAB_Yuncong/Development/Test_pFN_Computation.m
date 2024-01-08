% Yuncong Ma, 8/14/2023
% This script is to test the corresponding python file, in order to ensure
% the same result of pFN computation
% Test_pFN_Computation.py

%% Add to path
addpath(genpath('/Users/yuncongma/Documents/Document/fMRI/Myworks/pNet'));
%%
maxNumCompThreads(10);
tic
% load file
scan = fLoad_MATLAB_Single_Variable('/Users/yuncongma/Documents/Document/fMRI/Myworks/pNet/Example/HCP_Surface/Data/100206/1/LR/Image.mat');
gNb = fLoad_MATLAB_Single_Variable('/Users/yuncongma/Documents/Document/fMRI/Myworks/pNet/Example/HCP_Surface/Test_FN17/FN_Computation/gNb.mat');
gFN = fLoad_MATLAB_Single_Variable('/Users/yuncongma/Documents/Document/fMRI/Myworks/pNet/Example/HCP_Surface/Test_FN17/Group_FN/FN.mat');
Setting = fLoad_MATLAB_Single_Variable('/Users/yuncongma/Documents/Document/fMRI/Myworks/pNet/Example/HCP_Surface/Test_FN17/FN_Computation/Setting.mat');

scan=single(scan);
gFN=single(gFN);
%%

% parameter
K = 17;
maxIter = 1000;
minIter = 30;
meanFitRatio = 0.1;
error = 1e-4;
nRepeat = 1;
rounds = 30;
normW = 1;
eta = 0;
alphaS = 2;
alphaL = 10;
vxI = 0;
initConv = 1;
ard = 0;
dataPrecision = 'double';
logFile = 'Log_pFN_NMF_MATLAB.log';

% setup log file
% fileID = fopen(logFile, 'w');
fprintf(['\nStart NMF for pFN using MATLAB at ',char(datetime),'\n'])


% [mat_float,mat_eps] = set_data_precision(dataPrecision);
mat_eps = eps(dataPrecision);
if strcmp(dataPrecision,'single')
    gFN=single(gFN);
    scan=single(scan);
else
    gFN=double(gFN);
    scan=double(scan);
end

% initialization
initV=gFN;

dim_space = size(scan,1);
dim_time = size(scan,2);

% Setup alphaS and alphaL
nmVec = zeros(length(gNb),1);
for gni=1:length(gNb)
    nmVec(gni) = length(gNb{gni});
end
nM = median(nmVec);

% Use Alpha and Beta to set alphaS and alphaL if they are 0
if alphaS==0 && Alpha>0
    alphaS = round((Alpha*dim_time)/K);
end
if alphaL==0 && Beta>0
    alphaL = round((Beta*dim_time)/(K*nM));
end

% normalize data
X = normalize_data(scan','vp','vmax');


% construct the spatial affinity graph
tmpW = sparse(dim_space,dim_space);
for vi=1:dim_space
    nei = gNb{vi};
    if vxI>0
        corrVal = (1+corr(X(:,vi),X(:,nei)))/2;
    else
        corrVal = 1;
    end
    tmpW(vi,nei) = corrVal;
    tmpW(nei,vi) = corrVal;
end
W = tmpW;


DCol = full(sum(W,2));
D = spdiags(DCol,0,dim_space,dim_space);
L = D - W;
if normW>0
    D_mhalf = spdiags(DCol.^-0.5,0,dim_space,dim_space);
    L = D_mhalf * L * D_mhalf * alphaL;
    W = D_mhalf * W * D_mhalf * alphaL;
    D = D_mhalf * D * D_mhalf * alphaL;
end

% initialize V
V = initV;
miv = max(V,[],1);
trimInd = V ./ max(repmat(miv,dim_space,1),mat_eps) < 5e-2;
V(trimInd)=0;


% initialize U
U = X * V ./ repmat(sum(V,1),[dim_time,1]);

U = initialize_u(X, U, V, error, maxIter, minIter, meanFitRatio, initConv);

initU=U;
initV=V;

%
% Alternative update of U and V

oldL = Inf;
j = 0;
iterLog = [];

% No reuse of initU and initV, can be reference
U = initU;
V = initV;

[dim_time,dim_space] = size(X);

if ard>0
    lambdas = sum(U,1) / dim_time;
    hyperLam = eta * sum(sum(X.^2,1),2) / (dim_time*dim_space*2);
else
    lambdas=0;
    hyperLam=0;
end

flagQC=0;
oldLogL=inf;
oldU=U;
oldV=V;
% Multiplicative update of U and V
for i=1:maxIter

    % ===================== update V ========================
    % Eq. 8-11
    XU = X'*U;
    UU = U'*U;
    VUU = V*UU;

    tmpl2 = V.^2;

    if alphaS>0
        tmpl21 = sqrt(tmpl2);
        tmpl22 = repmat(sqrt(sum(tmpl2,1)),dim_space,1);
        tmpl21s = repmat(sum(tmpl21,1),dim_space,1);
        posTerm = V ./ max(tmpl21.*tmpl22,mat_eps);
        negTerm = V .* tmpl21s ./ max(tmpl22.^3,mat_eps);

        VUU = VUU + 0.5 * alphaS * posTerm;
        XU = XU + 0.5 * alphaS * negTerm;
    end

    if alphaL>0
        WV = W * double(V);
        DV = D * double(V);

        XU = XU + WV;
        VUU = VUU + DV;
    end

    V = V.*(XU./max(VUU,mat_eps));

    prunInd = sum(V~=0,1)==1;
    if any(prunInd)
        V(:,prunInd) = zeros(dim_space,sum(prunInd));
        U(:,prunInd) = zeros(dim_time,sum(prunInd));
    end

    % ==== normalize U and V ====
    [U,V] = normalize_u_v(U, V, 1, 1);

    % ===================== update U =========================
    XV = X*V;
    VV = V'*V;
    UVV = U*VV;

    if ard>0 % ard term for U
        posTerm = 1./max(repmat(lambdas,dim_time,1),mat_eps);
        UVV = UVV + posTerm*hyperLam;
    end

    U = U.*(XV./max(UVV,mat_eps));

    prunInd = sum(U)==0;
    if any(prunInd)
        V(:,prunInd) = zeros(dim_space,sum(prunInd,1));
        U(:,prunInd) = zeros(dim_time,sum(prunInd,1));
    end

    % update lambda
    if ard>0
        lambdas = sum(U,1) / dim_time;
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

    % QC Control
    temp=corr(initV,V);
    QC_Spatial_Correspondence=diag(temp);
    temp=temp-diag(diag(temp));
    QC_Spatial_Correspondence_Control=max(temp,[],2);
    QC_Delta_Sim=min(QC_Spatial_Correspondence-QC_Spatial_Correspondence_Control);
    if QC_Delta_Sim<=0
        U=oldU;
        V=oldV;
        flagQC=1;
        fprintf('\n  Meet QC constraint: Delta sim = %f\n',QC_Delta_Sim);
        fprintf('    Use results from last iteration\n');
        break;
    else
        oldU=U;
        oldV=V;
        fprintf('\n    Delta sim = %f\n',QC_Delta_Sim);
    end

end

fprintf(['\n Finished at ',char(datetime),'\n'])
% save();
toc



function [Flag,Message]=gFN_SR_NMF_Bootstrapping(App_Dir,Work_Dir,FN_Computation,Bootstrapping)
% Yuncong Ma, 2/1/2024
% Perform bootstrapping computation for group-level FN
% [Flag,Message]=gFN_SR_NMF_Bootstrapping(App_Dir,Work_Dir,FN_Computation,Bootstrapping)
%

Flag=0;
Message='';

FID=fopen(Bootstrapping.File,'r');
Scan_List=textscan(FID,'%s\n');
Scan_List=Scan_List{1};
fclose(FID);

% Load data
[Data,Flag,Message]=fLoad_Scan(App_Dir,Work_Dir,Scan_List);
if Flag==1
    return
end
[gNb,Flag,Message]=fLoad_MATLAB_Single_Variable(fullfile(Work_Dir,'FN_Computation','gNb.mat'));
if Flag
    return
end

% Parameter
K=FN_Computation.K;
spaR=FN_Computation.GroupFN.spaR;
vxI=FN_Computation.GroupFN.vxI;
ard=0;
iterNum=FN_Computation.GroupFN.iterNum;
alpha=FN_Computation.GroupFN.Alpha;
beta=FN_Computation.GroupFN.Beta;

nmVec = zeros(length(gNb),1);
for gni=1:length(gNb)
    nmVec(gni) = length(gNb{gni});
end
nM = median(nmVec);

numUsed = sum(cellfun(@(x)size(x,1),Data));
pS = round((alpha*numUsed)/K);
pL = round((beta*numUsed)/(K*nM));


options = [];
options.maxIter = iterNum;
options.error = 1e-6;
options.nRepeat = 1;
options.minIter = 100; %iterNum;
options.Converge = 1;
options.meanFitRatio = 0.1;

options.S1 = pS;
options.L = pL;
options.robust = 0; % robust NMF


neiR = spaR;
vxlInfo = vxI;

if ard~=0
    options.ard = 1;
end

% Start computation

nSmp = size(Data{1},2);
totFea = 0;
mFea = zeros(numUsed,1);
for si=1:length(Data)
    mFea(si) = size(Data{si},1);
    totFea = totFea + mFea(si);
end
catX = zeros(totFea,nSmp,'single');

sInd = 1;
for si=1:length(Data)

    origSbjData = dataPrepro(Data{si},'vp','vmax');

    eInd = sInd + mFea(si) - 1;
    catX(sInd:eInd,:) = origSbjData;
    sInd = eInd + 1;
end
clear Data;


% construct the affinity graph here
vNum = size(gNb,1);
catW = sparse(vNum,vNum);
for vi=1:vNum
    for ni=1:length(gNb{vi})
        nei = gNb{vi}(ni);
        if vxlInfo~=0
            if vi<nei
                corrVal = (1+corr(catX(:,vi),catX(:,nei)))/2;
            else
                continue;
            end
        else
            corrVal = 1;
        end
        if isnan(corrVal)
            corrVal = 0;
        end
        catW(vi,nei) = corrVal;
        catW(nei,vi) = corrVal;
    end
end


[~, FN] = fmNMF_sp(catX, K, catW, options, Bootstrapping.Out_Dir);
% initU and initV in old variable naming

%Output
[FN,Flag,Message]=fReformat_FN(App_Dir,Work_Dir,FN);
if Flag
    return
end
save(fullfile(Bootstrapping.Out_Dir,'FN.mat'),'FN','-v7.3');

end


function [U_final, V_final, nIter_final, elapse_final, bSuccess, objhistory_final] = fmNMF_sp(X, k, W, options, Out_Dir)
% Non-negative Matrix Factorization (NMF) with multiplicative update
% with sparse and graph constraints
% Notation:
% X ... (mFea x nSmp) data matrix
%       mFea  ... number of time points
%       nSmp  ... number of voxels
% k ... number of ICN
% W ... weight matrix of the sample affinity graph
% options ... Structure holding all settings
%
% U_ ... initialization for time series
% V_ ... initialization for spatial maps
%
%   Written by Deng Cai (dengcai AT gmail.com)
%   Modified by Jialu Liu (jliu64 AT illinois.edu)
%   Modified by hmli
%  Modified by Yuncong Ma, 2/22/2023, add option to use litekmeans to initialize U and V

Debug=1;
FID=fopen(fullfile(Out_Dir,'Log.txt'),'w');

fprintf(FID,'\n--mNMF_sp--\n');
fprintf(FID,char(datetime('now')));
fprintf(FID,'\n');

% Parameter
differror = options.error;
maxIter = options.maxIter;
nRepeat = options.nRepeat;
minIterOrig = options.minIter;
minIter = minIterOrig-1;
meanFitRatio = options.meanFitRatio;

alphaL = options.L;
alphaS = options.S1;

Norm = 1;
NormV = 1;

[mFea,nSmp] = size(X);

if alphaL > 0
    DCol = full(sum(W,2));
    D = spdiags(DCol,0,nSmp,nSmp);
    L = D - W;

    D_mhalf = spdiags(DCol.^-0.5,0,nSmp,nSmp) ;
    L = D_mhalf*L*D_mhalf * alphaL;
    W = D_mhalf*W*D_mhalf * alphaL;
    D = D_mhalf*D*D_mhalf * alphaL;
else
    L = [];
end

bSuccess.bSuccess = 1;

fprintf(FID,'\nStart initialization\n');

Method_Choice='Original';
switch Method_Choice
    case 'Original'
        U_=[];
        V_=[];
    case 'kmeans_Least_Square'
        % Yuncong Ma, 12/1/2022
        % X is [T,F]
        T_Reduction=10;
        T_ps=randi(mFea,[round(mFea/T_Reduction),1]);
        [label, ~] = litekmeans(X(T_ps,:)',k,'Replicates',20);
        U_=zeros(mFea,k);
        for i=1:k
            U_(:,i)=mean(X(:,label==i),2);
        end
        %             V_=zeros(nSmp,k);
        %             for i=1:k
        %                 V_(:,i)=(label==i)+0.1;
        %             end
        V_=X'*U_*pinv(U_'*U_);
        V_ = dataPrepro(V_,'vp','vmax');
        V_=max(V_,eps);
end
%

selectInit = 1;
mean_X = sum(X(:))/(mFea*nSmp);
if isempty(U_)
    U = (rand(mFea,k)+1)*(sqrt(mean_X/k));
    if isempty(V_)
        V = (rand(nSmp,k)+1)*(sqrt(mean_X/k));
    else
        V = V_;
    end
else
    U = U_;
    if isempty(V_)
        V = (rand(nSmp,k)+1)*(sqrt(mean_X/k));
    else
        V = V_;
    end
end

[U,V] = NormalizeUV(U, V, NormV, Norm);

if Debug==1
    fprintf(FID,'size U is %d %d\n',size(U));
    fprintf(FID,'size V is %d %d\n',size(V));
    fprintf(FID,'\n');
end

if isfield(options,'ard') && options.ard==1
    ard = 1;
    eta = 0.1;
    lambdas = sum(U) / mFea;
    hyperLam = eta * sum(sum(X.^2)) / (mFea*nSmp*2);
else
    ard = 0;
    hyperLam = 0;
end

if nRepeat == 1
    selectInit = 0;
    minIterOrig = 0;
    %minIter = 0;
    if isempty(maxIter)
        objhistory = CalculateObj(X, U, V, L, alphaS, ard, hyperLam);
        meanFit = objhistory*10;
    else
        if isfield(options,'Converge') && options.Converge
            objhistory = CalculateObj(X, U, V, L, alphaS, ard, hyperLam);
            meanFit = objhistory*10;
        end
    end
else
    if isfield(options,'Converge') && options.Converge
        error('Not implemented!');
    end
end

fprintf(FID,'\nStart iteration\n');

tryNo = 0;
while tryNo < nRepeat
    tmp_T = cputime;
    tryNo = tryNo+1;
    nIter = 0;
    maxErr = 1;
    nStepTrial = 0;
    while(maxErr > differror)
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
            posTerm = 1 ./ max(repmat(tmpNorm2,nSmp,1),eps);
            tmpNorm1 = sum(V,1);
            negTerm = V .* repmat(tmpNorm1,nSmp,1) ./ max(repmat(tmpNorm2.^3,nSmp,1),eps);

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

        V = V.*(XU./max(VUU,eps));

        prunInd = sum(V~=0)==1;
        if any(prunInd)
            V(:,prunInd) = zeros(nSmp,sum(prunInd));
            U(:,prunInd) = zeros(mFea,sum(prunInd));
        end

        % erease very small numbers
        %ereEle = V./max(eps,repmat(max(V),size(V,1),1));
        %V(ereEle<1e-6) = 0;

        [U, V] = NormalizeUV(U, V, NormV, Norm);

        % ===================== update U ========================
        XV = X*V;   % mnk or pk (p<<mn)
        VV = V'*V;  % nk^2
        UVV = U*VV; % mk^2

        % needed if Ui is restricted to norm 1
        %         if alphaS>0
        %             % L1 sparsity
        %             F = repmat(sum(V), mFea, 1);
        %             UVV = UVV + 0.5*alphaS * F;
        %         end
        %         if alphaL>0
        %             V = double(V);
        %             VLV = repmat(diag(V'*L*V)' .* sum(U,1), mFea, 1);
        %             UVV = UVV + VLV;
        %         end

        if isfield(options,'ard') && options.ard==1
            posTerm = 1./max(repmat(lambdas,mFea,1),eps);
            UVV = UVV + posTerm*hyperLam;
        end

        U = U.*(XV./max(UVV,eps));

        prunInd = sum(U)==0;
        if any(prunInd)
            V(:,prunInd) = zeros(nSmp,sum(prunInd));
            U(:,prunInd) = zeros(mFea,sum(prunInd));
        end

        if isfield(options,'ard') && options.ard==1
            lambdas = sum(U) / mFea;
        end

        nIter = nIter + 1;
        %disp(['   iteration:',num2str(nIter)]);
        if nIter > minIter
            if selectInit
                objhistory = CalculateObj(X, U, V, L, alphaS, ard, hyperLam);
                maxErr = 0;
            else
                if isempty(maxIter)
                    newobj = CalculateObj(X, U, V, L, alphaS, ard, hyperLam);
                    objhistory = [objhistory newobj];
                    meanFit = meanFitRatio*meanFit + (1-meanFitRatio)*newobj;
                    maxErr = (meanFit-newobj)/meanFit;
                else
                    if isfield(options,'Converge') && options.Converge
                        [newobj, newObjStr] = CalculateObj(X, U, V, L, alphaS, ard, hyperLam);
                        objhistory = [objhistory newobj];
                        if mod(nIter,10)==0
                            fprintf(FID,['\n  iter ',num2str(nIter), ' in ',num2str(maxIter),'\n']);
                            fprintf(FID,[newObjStr,'\n\n']);
                        end

                        %meanFit = meanFitRatio*meanFit + (1-meanFitRatio)*newobj;
                        %maxErr = (meanFit-newobj)/meanFit;
                        maxErr = abs(objhistory(end)-objhistory(end-1))/objhistory(end);
                    else
                        maxErr = 1;
                    end

                    if nIter >= maxIter
                        maxErr = 0;
                        if isfield(options,'Converge') && options.Converge
                        else
                            objhistory = 0;
                        end
                    end
                end
            end
        end
    end

    elapse = cputime - tmp_T;

    if tryNo == 1
        U_final = U;
        V_final = V;
        nIter_final = nIter;
        elapse_final = elapse;
        objhistory_final = objhistory;
        bSuccess.nStepTrial = nStepTrial;
    else
        if objhistory(end) < objhistory_final(end)
            U_final = U;
            V_final = V;
            nIter_final = nIter;
            objhistory_final = objhistory;
            bSuccess.nStepTrial = nStepTrial;
            if selectInit
                elapse_final = elapse;
            else
                elapse_final = elapse_final+elapse;
            end
        end
    end

    if selectInit
        if tryNo < nRepeat
            %re-start
            if isempty(U_)
                U = abs(rand(mFea,k));
                norms = sqrt(sum(U.^2,1));
                norms = max(norms,eps);
                U = U./repmat(norms,mFea,1);
                if isempty(V_)
                    V = abs(rand(nSmp,k));
                    V = V/sum(sum(V));
                else
                    V = V_;
                end
            else
                U = U_;
                if isempty(V_)
                    V = abs(rand(nSmp,k));
                    V = V/sum(sum(V));
                else
                    V = V_;
                end
            end

            [U,V] = NormalizeUV(U, V, NormV, Norm);
        else
            tryNo = tryNo - 1;
            minIter = 0;
            selectInit = 0;
            U = U_final;
            V = V_final;
            objhistory = objhistory_final;
            meanFit = objhistory*10;
        end
    end

    %display(['  tryNo: ',num2str(tryNo),' ','obj: ',num2str(objhistory(end))]);% added by hmli
end

nIter_final = nIter_final + minIterOrig;

[U_final, V_final] = NormalizeUV(U_final, V_final, NormV, Norm);
objhistory_final = CalculateObj(X, U_final, V_final, L, alphaS, ard, hyperLam); % added by hmli

fprintf(FID,['\nDone at ',char(datetime('now'))]);
fclose(FID);

end



function [obj, objStr] = CalculateObj(X, U, V, L, alphaS, ard, hyperLam, deltaVU, dVordU)
if ~exist('deltaVU','var')
    deltaVU = 0;
end
if ~exist('dVordU','var')
    dVordU = 1;
end
dV = [];
maxM = 62500000;
[mFea, nSmp] = size(X);
mn = numel(X);
nBlock = floor(mn*3/maxM);

if mn < maxM
    dX = U*V'-X;
    obj_NMF = sum(sum(dX.^2));
    if deltaVU
        if dVordU
            dV = dX'*U;
        else
            dV = dX*V;
        end
    end
else
    obj_NMF = 0;
    if deltaVU
        if dVordU
            dV = zeros(size(V));
        else
            dV = zeros(size(U));
        end
    end
    for i = 1:ceil(nSmp/nBlock)
        if i == ceil(nSmp/nBlock)
            smpIdx = (i-1)*nBlock+1:nSmp;
        else
            smpIdx = (i-1)*nBlock+1:i*nBlock;
        end
        dX = U*V(smpIdx,:)'-X(:,smpIdx);
        obj_NMF = obj_NMF + sum(sum(dX.^2));
        if deltaVU
            if dVordU
                dV(smpIdx,:) = dX'*U;
            else
                dV = dU+dX*V(smpIdx,:);
            end
        end
    end
    if deltaVU
        if dVordU
            dV = dV ;
        end
    end
end
if isempty(L)
    obj_Lap = 0;
else
    V = double(V);
    obj_Lap = sum(sum((L*V).*V));
end

%obj_Spa = alphaS * sum(sum(V));
tmpNorm1 = sum(V,1);
tmpNorm2 = sqrt(sum(V.^2,1)) + eps;
obj_Spa = alphaS * sum(tmpNorm1./tmpNorm2);

if ard>0
    su = sum(U);
    su(su==0) = 1;
    obj_ard = sum(log(su))*mFea*hyperLam;
else
    obj_ard = 0;
end

obj = obj_NMF + obj_Lap + obj_Spa;
objStr = ['    totObj:',num2str(obj),',NMF:',num2str(obj_NMF),',Lap:',num2str(obj_Lap),',Spa:',num2str(obj_Spa),',Ard:',num2str(obj_ard)];

end

% function [U, V] = Normalize(U, V)
%     [U,V] = NormalizeUV(U, V, 1, 1);
function [U, V] = NormalizeUV(U, V, NormV, Norm)
nSmp = size(V,1);
mFea = size(U,1);
if Norm == 2
    if NormV
        norms = sqrt(sum(V.^2,1));
        norms = max(norms,eps);
        V = V./repmat(norms,nSmp,1);
        U = U.*repmat(norms,mFea,1);
    else
        norms = sqrt(sum(U.^2,1));
        norms = max(norms,eps);
        U = U./repmat(norms,mFea,1);
        V = V.*repmat(norms,nSmp,1);
    end
else
    if NormV
        %norms = sum(abs(V),1);
        norms = max(V);
        norms = max(norms,eps);
        V = V./repmat(norms,nSmp,1);
        U = U.*repmat(norms,mFea,1);
    else
        %norms = sum(abs(U),1);
        norms = max(U);
        norms = max(norms,eps);
        U = U./repmat(norms,mFea,1);
        V = V.*repmat(norms,nSmp,1);
    end
end
end







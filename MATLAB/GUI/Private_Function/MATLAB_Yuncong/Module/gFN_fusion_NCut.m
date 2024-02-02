function [Flag,Message]=gFN_fusion_NCut(App_Dir,Work_Dir,Compute_FN,Bootstrapping)
% Yuncong Ma, 2/1/2024
% Fuse the bootstrapping results for group-level gFN
% [Flag,Message]=gFN_fusion_NCut(Work_Dir,Compute_FN,Bootstrapping)

Flag=0;
Message='';

K=Compute_FN.K;
repNum = length(Bootstrapping);

% load results
resSet = [];
for ri=1:repNum
    resName = fullfile(Bootstrapping(ri).Out_Dir,'FN.mat');
    [FN,Flag,Message]=fReformat_FN(App_Dir,Work_Dir,resName);
    if Flag
        return
    end

    resSet = [resSet, FN];
    clear FN;
end
resSet = resSet';

% clustering by ncut
corrVal = corr(resSet');
corrVal(isnan(corrVal)) = -1;
nDis = 1 - corrVal;
triuInd = triu(ones(size(nDis)),1);
nDisVec = nDis(triuInd==1);

nW = exp(-nDis.^2 ./ (median(nDisVec).^2));
nW(isnan(nW)) = 0;

sumW = sum(nW,1);
sumW(sumW==0) = 1;
D = diag(sumW);
L = sqrt(inv(D))*nW*sqrt(inv(D));
L = (L+L')/2;
opts.disp = 0;
[Ev,~] = eigs(double(L),K,'LA',opts);
normvect = sqrt(diag(Ev*Ev'));
normvect(normvect==0.0) = 1;
Ev = diag(normvect) \ Ev;

% Repeat clustering if empty results are found
N_Max_Test=50;
C=[];
Count=0;
while length(unique(C))<K && Count<N_Max_Test
    [EvDiscrete,~] = discretisation(Ev);
    EvDiscrete = full(EvDiscrete);
    [~, C] = max(EvDiscrete,[],2);

%     if length(unique(C))<K && Debug==1
%         fprintf('Find %d empty FN, redo discretisation\n',K-length(unique(C)));
%     end
    Count=Count+1;
end
if length(unique(C))<K
    Flag=1;
    Message=['Error in fCompute_gFN_Fusion: cannot generate ',num2str(K),' FN'];
end

% get centroid
initV = zeros(size(resSet,2),K);
for ki=1:K
    % % typical point
    if sum(C==ki)>1
        candSet = resSet(C==ki,:)';
        corrW = abs(corr(candSet));
        corrW(isnan(corrW)) = 0;
        [mVal, mInd] = max(sum(corrW));
        initV(:,ki) = candSet(:,mInd);
    elseif sum(C==ki)==1
        initV(:,ki) = resSet(C==ki,:);
    end
end
initV = initV ./ max(eps,repmat(max(initV),size(initV,1),1));

% Output
FN=initV;
[FN,Flag,Message]=fReformat_FN(App_Dir,Work_Dir,FN);
if Flag
    return
end
save(fullfile(Work_Dir,'Group_FN','FN.mat'),'FN');

end



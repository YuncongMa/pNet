% Yuncong Ma, 8/26/2023
% Test the fusion of gFN results

%% Add to path
addpath(genpath('/Users/yuncongma/Documents/Document/fMRI/Myworks/pNet'));

%%

% Setup parameters
dir_pnet='/Users/yuncongma/Documents/Document/fMRI/Myworks/pNet';
folder_gFN = fullfile(dir_pnet,'Example/HCP_Surface/Test_FN17_200', 'Group_FN', 'BootStrapping');
file_output = fullfile(dir_pnet,'Example/HCP_Surface/Test_FN17_200', 'Group_FN', 'FN.mat');
file_setting = fullfile(dir_pnet,'Example/HCP_Surface/Test_FN17_200', 'FN_Computation', 'Setting.json');

% Load setting parameters
fid = fopen(file_setting);
raw = fread(fid);
str = char(raw');
fclose(fid);
setting = jsondecode(str);
K = setting.K;
BS_Repetition = 10; % Number of bootstrapping
NCut_MaxTrial=1; % Max number of trials
dataPrecision = 'double';

mat_eps = eps(dataPrecision);

% Concatenate FN along the FN dimension, gFN_BS is [dim_space, K * n_BS]
gFN_BS=[];
for i=1:BS_Repetition
    FN = fLoad_MATLAB_Single_Variable(fullfile(folder_gFN,num2str(i),'FN.mat'));
    if strcmp(dataPrecision,'single')
        FN=single(FN);
    else
        FN=double(FN);
    end
    gFN_BS=[gFN_BS, FN];
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clustering by NCut

% Get similarity between samples
corrVal = corr(gFN_BS); % similarity between FNs, [K * n_BS, K * n_BS]
corrVal(isnan(corrVal)) = -1;
nDis = 1 - corrVal; % Transform Pearson correlation to non-negative values similar to distance
triuInd = triu(ones(size(nDis)),1); % Index of upper triangle
nDisVec = nDis(triuInd==1); % Take the values in upper triangle
% Make all distance non-negative and normalize their distribution
nW = exp(-nDis.^2 ./ (median(nDisVec)^2)); % Transform distance values using exp(-X/std^2) with std as the median value
nW(isnan(nW)) = 0;
sumW = sum(nW,1); % total distance for each FN
sumW(sumW==0) = 1; % In case two FNs are the same
% Construct Laplacian matrix
D = diag(sumW);
L = sqrt(inv(D))*nW*sqrt(inv(D)); % A way to normalize nW based on the total distance of each FN
L = (L+L')/2; % Ensure L is symmetric. Computation error may result in asymmetry
opts.disp = 0;
[Ev,eVal] = eigs(double(L),K,'largestreal',opts); % Get first K eigenvectors
% Correct the sign of eigenvectors to make them same as derived from Python 
temp=sign(sum(Ev,1)); % Use the total value of each eigenvector to reset its sign
temp(temp==0)=1;
Ev=Ev.*repmat(temp,[size(Ev,1),1]); % Reset the sign of each eigenvector
normvect = sqrt(diag(Ev*Ev')); % Get the norm of each eigenvector
normvect(normvect==0.0) = 1; % Incase all 0 eigenvector
Ev = diag(normvect) \ Ev; % Use linear solution to normalize Ev satisfying normvect * Ev_new = Ev_old

% Multiple trials to get reproducible results
Best_C=[];
Best_NcutValue=inf;
for i=1:NCut_MaxTrial

    EigenVectors=Ev;
    [n,k]=size(EigenVectors); % n is K * n_BS, k is K

    vm = sqrt(sum(EigenVectors.^2,2)); % norm of each row
    EigenVectors = EigenVectors./repmat(vm,1,k); % normalize eigenvectors to ensure each FN vector's norm = 1

    R=zeros(k);
    ps=1+round(rand(1)*(n-1)); % Choose a random row in eigenvectors
    R(:,1)=EigenVectors(ps,:)'; % This randomly selected row in eigenvectors is used as an initial center
    c=zeros(n,1); % Total distance to different rows in R [K * n_BS, 1]
    c_index=zeros(n,1); % Store the index of selected samples
    c_index(1)=ps;
    c(ps)=inf;

    for j=2:k % Find another K-1 rows in eigenvectors which have the minimum similarity to previous selected rows, similar to initialization in k++ 
        c=c+abs(EigenVectors*R(:,j-1));
        [~,ps]=min(c);
        ps=ps(1); % in case more than one sample has the minimum distance
        c_index(j)=ps;
        c(ps(1))=inf;
        R(:,j)=EigenVectors(ps(1),:)';
    end
    
    lastObjectiveValue=0;
    exitLoop=0;
    nbIterationsDiscretisation = 0;
    nbIterationsDiscretisationMax = 20;
    while exitLoop == 0
        nbIterationsDiscretisation = nbIterationsDiscretisation + 1 ;
        
        EigenVectorsR=EigenVectors*R;
        [n,k]=size(EigenVectorsR);
        [~,J]=max(EigenVectorsR,[],2); % Assign each sample to K centers of R based on highest similarity

        EigenvectorsDiscrete=sparse(1:n,J',1,n,k); % Generate a 0-1 matrix with each row containing only one 1

        [U,S,V] = svd(EigenvectorsDiscrete'*EigenVectors,0); % Economy-size decomposition
        NcutValue=2*(n-trace(S));
        
        % escape the loop when converged or meet max iteration
        if abs(NcutValue-lastObjectiveValue) < mat_eps || nbIterationsDiscretisation > nbIterationsDiscretisationMax
            exitLoop=1;
            fprintf('\nReach stop criterion of NCut, NcutValue = %f\n',NcutValue);
        else
            fprintf('   NcutValue = %f\n',NcutValue);
            lastObjectiveValue = NcutValue;
            R=V*U'; % Update R which stores the new centers
        end
    end

    [~, C] = max(full(EigenvectorsDiscrete),[],2); % Assign each sample to K centers in R

    if length(unique(C))<K % Check whether there are empty results
        fprintf('Find empty results in iteration %d\n',i);
    else  % Update the best result
        if NcutValue<Best_NcutValue
            Best_NcutValue=NcutValue;
            Best_C=C;
        end
    end
end
if isempty(Best_C) % In case even the last trial has empty results
    Flag=1;
    Message='Cannot generate non-empty FN';
end
fprintf('Best NCut Value = %f\n',Best_NcutValue);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get centroid
C = Best_C;
gFN = zeros(size(gFN_BS,1),K);
for ki=1:K
    % typical point
    if sum(C==ki)>1
        candSet = gFN_BS(:, C==ki); % Get the candidate set of FNs assigned to cluster ki
        corrW = abs(corr(candSet)); % Get the similarity between candidate FNs
        corrW(isnan(corrW)) = 0;
        [mVal, mInd] = max(sum(corrW,1),[],2); % Find the FN with the highest total similarity to all other FNs
        gFN(:,ki) = candSet(:,mInd);
    elseif sum(C==ki)==1
        gFN(:,ki) = gFN_BS(:, C==ki);
    end
end

gFN = gFN ./ max(mat_eps,repmat(max(gFN,[],1),size(gFN,1),1)); % Normalize each FN by its max value

FN=gFN;

%% Test reproducibility

clear Result
Result(1:100)=struct('gFN',[]);
for i=1:100
    % NCut
    Result(i).gFN = fTest_gFN_Fusion(gFN_BS, K, 100, dataPrecision);
    % K-means
%     [~, Result(i).gFN] = kmeans(gFN_BS',K,'Replicates',100);
%     Result(i).gFN=Result(i).gFN';
    % SC
end

Reproducibility.NCut=zeros(100);
for i=1:100
    for j=i+1:100
        Reproducibility.NCut(i,j)=mean(max(corr(Result(i).gFN,Result(j).gFN)));
    end
end

Result(1:100)=struct('gFN',[]);
for i=1:100
    % NCut
%     Result(i).gFN = fTest_gFN_Fusion(gFN_BS, K, 100, dataPrecision);
    % K-means
    [~, Result(i).gFN] = kmeans(gFN_BS',K,'Replicates',100);
    Result(i).gFN=Result(i).gFN';
end

Reproducibility.kmeans=zeros(100);
for i=1:100
    for j=i+1:100
        Reproducibility.kmeans(i,j)=mean(max(corr(Result(i).gFN,Result(j).gFN)));
    end
end

Result(1:100)=struct('gFN',[]);
for i=1:100
    % Spectral clustering
    Result(i).gFN = fTest_gFN_Fusion_SC(gFN_BS,K);
end

Reproducibility.SC=zeros(100);
for i=1:100
    for j=i+1:100
        Reproducibility.SC(i,j)=mean(max(corr(Result(i).gFN,Result(j).gFN)));
    end
end

fFigure(1,1,1,'',[500,400]);
fAxes(1,1,1,[.7,.6],[0,.1]);
temp=[fApply_Mask('Upper',Reproducibility.NCut)',fApply_Mask('Upper',Reproducibility.kmeans)',fApply_Mask('Upper',Reproducibility.SC)'];
p=signrank(temp(:,1),temp(:,2))
violinplot(temp,1,'ShowData',false);
xticks([1:3]);
xticklabels({'NCut','k-means','Spectral Clustering'})
xlim([0,3]+0.5);
ylabel('Reproducibility');
set(gca,'fontsize',16,'LineWidth',2);
set(gcf,'color','w');



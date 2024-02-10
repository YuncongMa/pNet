function gFN = fTest_gFN_Fusion_SC(gFN_BS, K)
% Yuncong Ma, 8/30/2023
% A function to test the robustness of gFN fusion
% Spectral clustering

mat_eps = eps('double');

% clustering by SC



% Get centroid
C = spectralcluster(gFN_BS',K);
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

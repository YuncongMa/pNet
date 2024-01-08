function pX = normalize_data(data, algorithum, normalization)
% normalize_data(data, algorithm, normalization)
% Normalize data by algorithm and normalization settings
% :param data: data in 2D matrix [dim_time, dim_space]
% :param algorithm: 'z' 'gp' 'vp'
% :param normalization: 'n2' 'n1' 'rn1' 'g' 'vmax'
% :return: data
% consistent to python function pNet.normalize_data(X, algorithm, normalization)
% By Yuncong Ma, 8/11/2023

if length(size(data))~=2
    error('data must be a 2D matrix');
end

X=data;
np_eps = eps(class(X));

switch lower(algorithum)
    case 'z'
        % standard score for each variable
        mVec = mean(X,2);
        sVec = max(std(X,0,2),np_eps);
        pX = (X-repmat(mVec,1,size(X,2)))./repmat(sVec,1,size(X,2));
    case 'gp'
        % remove negative value globally
        minVal = min(X(:));
        shiftVal = abs(min(minVal,0));
        pX = X + shiftVal;
    case 'vp'
        % remove negative value voxel-wisely
        minVal = min(X,[],1);
        shiftVal = abs(min(minVal,0));
        pX = X + repmat(shiftVal,size(X,1),1);
    otherwise
        % do nothing
        disp('  unknown preprocess parameters, no preprocess applied');
        pX = X;
end

% normalization
switch lower(normalization)
    case 'n2'
        % l2 normalization for each observation
        l2norm = sqrt(sum(pX.^2)) + np_eps;
        pX = pX ./ repmat(l2norm,size(pX,1),1);
    case 'n1'
        % l1 normalization for each observation
        l1norm = sum(pX) + np_eps;
        pX = pX ./ repmat(l1norm,size(pX,1),1);
    case 'rn1'
        % l1 normalization for each variable
        l1norm = sum(pX,2) + np_eps;
        pX = pX ./ repmat(l1norm,1,size(pX,2));
    case 'g'
        % global scale
        [sVal,sInd] = sort(pX(:));
        perT = 0.001;
        minVal = sVal(round(length(sInd)*perT));
        maxVal = sVal(round(length(sInd)*(1-perT)));
        pX(pX<minVal) = minVal;
        pX(pX>maxVal) = maxVal;
        pX = (pX-minVal)/max((maxVal-minVal),np_eps);
    case 'vmax'
        cmin = repmat(min(pX),size(pX,1),1);
        cmax = repmat(max(pX),size(pX,1),1);
        pX = (pX-cmin)./max(np_eps,cmax-cmin);
    otherwise
        % do nothing
        disp('  unknown normalization parameters, no normalization applied');
end
if any(isnan(pX))
    error('  nan exists, check the preprocessed data');
end
end



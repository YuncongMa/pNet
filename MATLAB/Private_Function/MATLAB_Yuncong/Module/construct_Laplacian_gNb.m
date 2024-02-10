function [L,W,D]=construct_Laplacian_gNb(gNb, X, vxI, dim_space, alphaL, normW, dataPrecision)
% construct_Laplacian_gNb(gNb, X, vxI, dim_space, alphaL=10, normW=1)
% construct Laplacian matrices for Laplacian spatial regularization term
% :param gNb:
% :param X:
% :param vxI:
% :param dim_space:
% :param alphaL:
% :param normW:
% :param dataPrecision:
% :return: L, W, D: sparse 2D matrices [dim_space, dim_space]
% Yuncong Ma, 8/17/203

% [mat_float, mat_eps] = set_data_precision(dataPrecision);


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

end
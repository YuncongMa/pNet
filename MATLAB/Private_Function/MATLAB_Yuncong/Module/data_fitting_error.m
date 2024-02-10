function Fitting_Error = data_fitting_error(X, U, V, deltaVU, dVordU)
% data_fitting_error(X, U, V, deltaVU=None, dVordU=None)
% Calculate the datat fitting of X'=UV' with terms
% :param X: 2D matrix, [Space, Time]
% :param U: 2D matrix, [Time, k]
% :param V: 2D matrix, [Space, k]
% :param deltaVU: 0
% :param dVordU: 1
% :return: Fitting_Error
% Consistent to python function pNet.data_fitting_error(X, U, V, deltaVU, dVordU)
% By Yuncong Ma, 8/21/2023

if length(size(X))~=2 || length(size(U))~=2 || length(size(V))~=2
    error('X, U and V must be 2D matrices');
end
if size(X,1)~=size(U,1) || size(X,2)~=size(V,1) || size(U,2)~=size(V,2)
    error('X, U and V need to have apprpriate sizes');
end

dV = [];
maxM = 62500000; % To save memory
[dim_time, dim_space] = size(X);
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
    for i = 1:ceil(dim_space/nBlock)
        if i == ceil(dim_space/nBlock)
            smpIdx = (i-1)*nBlock+1:dim_space;
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

Fitting_Error = obj_NMF;

end
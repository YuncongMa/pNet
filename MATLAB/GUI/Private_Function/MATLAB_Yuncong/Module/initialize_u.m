function U_final = initialize_u(X, U0, V0, error, maxIter, minIter, meanFitRatio, initConv)
% initialize_u(X, U0, V0, error, maxIter, minIter, meanFitRatio, initConv)
% :param X: data in 2D matrix [dim_time, dim_space]
% :param U0: initial temporal component, 2D matrix [dim_time, k]
% :param V0: initial spatial component, 2D matrix [dim_space, k]
% :param error: fitting error
% :param maxIter: maximum iteration
% :param minIter: minimum iteration
% :param meanFitRatio: average fitting ratio
% :param initConv: flag for convergence
% :return: U_final
% Consistent to python function initialize_u(X, U0, V0, error, maxIter, minIter, meanFitRatio, initConv)
% By Yuncong Ma, 8/14/2023

if length(size(X))~=2 || length(size(U0))~=2 || length(size(V0))~=2
    error('X, U0 and V0 must be 2D matrices');
end
if size(X,1)~=size(U0,1) || size(X,2)~=size(V0,1) || size(U0,2)~=size(V0,2)
    error('X, U0 and V0 need to have apprpriate sizes');
end

np_eps = eps(class(X));

U = U0;
V = V0;

newFit = data_fitting_error(X, U, V, 0, 1);
meanFit = newFit/meanFitRatio;

maxErr = 1;
for i=1:maxIter
    % ===================== update V ========================
    %        XU = X'*U;  % mnk or pk (p<<mn)
    %        UU = U'*U;  % mk^2
    %        VUU = V*UU; % nk^2

    %        V = V.*(XU./max(VUU,eps));

    % ===================== update U ========================
    XV = X * V;   % mnk or pk (p<<mn)
    VV = V' * V;  % nk^2
    UVV = U * VV; % mk^2

    U = U .* (XV./max(UVV,np_eps)); % 3mk

    if i > minIter
        if initConv
            newFit = fitting_initialize_u(X, U, V, 0, 1);
            meanFit = meanFitRatio*meanFit + (1-meanFitRatio)*newFit;
            maxErr = (meanFit-newFit)/meanFit;
        end
    end
    if maxErr <= error
        break
    end
end

U_final = U;
end


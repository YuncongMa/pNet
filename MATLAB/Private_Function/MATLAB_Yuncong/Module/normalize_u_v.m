function [U, V] = normalize_u_v(U, V, NormV, Norm)
% normalize_u_v(U, V, NormV, Norm)
% Normalize U and V with terms
% :param U: 2D matrix, [Time, k]
% :param V: 2D matrix, [Space, k]
% :param NormV: 1 or 0
% :param Norm: 1 or 2
% :return: U, V
% Consistent to python function pNet.normalize_u_v(U, V, NormV, Norm)
% By Yuncong Ma, 8/14/2023

if length(size(U))~=2 || length(size(V))~=2
    error('U and V must be 2D matrices');
end
if size(U,2)~=size(V,2)
    error('U and V need to have apprpriate sizes');
end

np_eps = eps(class(U));

dim_space = size(V,1);
dim_time = size(U,1);

if Norm == 2
    norms = sqrt(sum(V.^2,1));
    norms = max(norms,np_eps);
else
    norms = max(V,[],1);
    norms = max(norms,np_eps);
end

if NormV
    U = U.*repmat(norms,dim_time,1);
    V = V./repmat(norms,dim_space,1);
else
    U = U./repmat(norms,dim_time,1);
    V = V.*repmat(norms,dim_space,1);
end
end
function [Q, R] = globalqr(X, param)
% [Q, R] = GLOBALQR(X, param) computes a QR factorization based on the
% global inner product.
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
% Default
if nargin == 1
    param.global_scale = 1;
end
if ~isfield(param, 'global_scale')
    param.global_scale = 1;
end

% Set scaling factor
s = size(X,2);
if param.global_scale
    scl = sqrt(s);
else
    scl = 1;
end

R = norm(X,'fro') / scl;
Q = X / R;
R = R * eye(s);

end
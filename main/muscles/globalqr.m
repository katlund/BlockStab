function [Q, R] = globalqr(X)
% [Q, R] = GLOBALQR(X) computes a QR factorization based on the
% global inner product.
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
% Set scaling factor
s = size(X,2);

R = norm(X,'fro') / s;
Q = X / R;
R = R * eye(s);

end
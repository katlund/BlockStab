function [Q, R] = globalqr(X)
% [Q, R] = GLOBALQR(X) computes a QR factorization based on the
% global inner product.
%
% Part of [BlockStab](https://github.com/katlund) package.  Check README
% for how to properly cite and reuse this file.

%%
% Set scaling factor
s = size(X,2);

R = norm(X,'fro') / s;
Q = X / R;
R = R * eye(s);

end
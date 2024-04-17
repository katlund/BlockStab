function [Q, R] = globalqr(X)
% [Q, R] = GLOBALQR(X) computes a QR factorization based on the
% global inner product.

%%
% Set scaling factor
s = size(X,2);

R = norm(X,'fro') / s;
Q = X / R;
R = R * eye(s);

end
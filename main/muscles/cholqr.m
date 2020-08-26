function [Q, R] = cholqr(X)
% [Q, R] = CHOLQR(X) computes Cholesky QR factorization of the m x s matrix
% X as described in [Stathopoulos & Wu 2002].

%%
[m, s] = size(X);
A = X' * X;
[~, flag] = chol(A);
if ~flag
    R = chol(A);
    Q = X / R;
else
    Q = NaN(m,s);
    R = NaN(s);
end
end
function [Q, R] = cholqr(X)
% [Q, R] = CHOLQR(X) computes Cholesky QR factorization of the m x s matrix
% X as described in [Stathopoulos & Wu 2002].
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
A = X' * X;
R = chol_nan(A);
Q = X / R;
end
function [Q, R] = cholqr_vpa(X)
% [Q, R] = CHOLQR_VPA(X) computes Cholesky QR factorization of the m x s
% matrix X with VPA.

%%
A = X' * X;
R = chol_free_vpa(A);
Q = X / R;
end
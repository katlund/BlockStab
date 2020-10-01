function [Q, R] = cholqr_free(X)
% [Q, R] = CHOLQR_FREE(X) computes the Cholesky QR factorization of the m x
% s matrix X as described in [Stathopoulos & Wu 2002] with CHOL_FREE
% instead of Matlab's built-in CHOL.

%%
A = X' * X;
R = chol_free(A);
Q = X / R;
end
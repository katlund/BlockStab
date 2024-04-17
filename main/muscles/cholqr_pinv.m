function [Q, R] = cholqr_pinv(X)
% [Q, R] = CHOLQR_PINV(X) computes the Cholesky QR factorization of the m x
% s matrix X with CHOL_FREE instead of Matlab's built-in CHOL, and using
% the Moore-Penrose pseudoinverse to invert R, instead of direct inversion.

%%
A = X' * X;
R = chol_free(A);
Q = (X / A) * R';   % Moore-Penrose pseudoinverse
end
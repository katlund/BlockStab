function [Q, R] = cholqr_mp(X)
% [Q, R] = CHOLQR_mp(X) computes Cholesky QR factorization of the m x s
% matrix X with MP (Advanpix).

%%
A = X' * X;
R = chol_free_mp(A);
Q = X / R;
end
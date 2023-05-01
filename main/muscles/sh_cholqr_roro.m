function [Q, R] = sh_cholqr_roro(X)
% [Q, R] = SH_CHOLQR_RORO(X, verbose) computes a shifted Cholesky QR
% factorization with reorthonormalization of the m x s matrix X.  This
% algorithm is equivalent to Algorithm 4.2 (shiftedCholeskyQR3) of [Fukaya,
% et al. 2020].
%
% See INTRAORTHO for more details about the parameters.
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
[m, s] = size(X);

% Shifted CholQR (round 1)
sh = 11 * (m*s + s*(s+1)) * eps * norm(X,2)^2;
A = X' * X + sh * eye(s);
[~, flag] = chol(A);
if flag == 0
    R1 = chol(A);
    Q = X / R1;
else
    Q = NaN(m,s);
    R = NaN(s);
    return;
end

% Rounds 2 and 3
[Q, R2] = IntraOrtho(Q, 'cholqr_ro');

% Combine
R = R2 * R1;
end
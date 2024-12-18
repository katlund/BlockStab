function [Q, R] = sh_cholqr_roro(X, param)
% [Q, R] = SH_CHOLQR_RORO(X, param) computes a shifted Cholesky QR
% factorization with reorthonormalization of the m x s matrix X.  This
% algorithm is equivalent to Algorithm 4.2 (shiftedCholeskyQR3) of [Fukaya,
% et al. 2020].
%
% See INTRAORTHO and CHOL_SWITCH for more details about the parameters.
%
% Part of [BlockStab](https://github.com/katlund/BlockStab) package.  Check README
% for how to properly cite and reuse this file.
%
% Part of [BlockStab](https://github.com/katlund/BlockStab) package.  Check README
% for how to properly cite and reuse this file.


%%
[m, s] = size(X);

% Shifted CholQR (round 1)
sh = 11 * (m*s + s*(s+1)) * eps * norm(X,2)^2;
A = X' * X + sh * eye(s);
R1 = chol_switch(A, param);
Q = X / R1;

% Rounds 2 and 3
[Q, R2] = cholqr_ro(Q, param);

% Combine
R = R2 * R1;
end
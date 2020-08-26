function [Q, R] = sh_cholqr_roro(X)
% [Q, R] = SH_CHOLQR_RORO(X) computes a shifted Cholesky QR factorization
% with reorthonormalization of the m x s matrix X.  This algorithm is
% equivalent to Algorithm 4.2 (shiftedCholeskyQR3) of [Fukaya, et. al.
% 2018].

%%
[m, s] = size(X);

% Shifted CholQR (round 1)
sh = 11*(m*s + s*(s+1))*eps*norm(X,2)^2;
A = X'*X + sh*eye(s);
[~, p] = chol(A);
if ~p
    R1 = chol(A);
    Q = X/R1;
else
    Q = NaN(m,s);
    R = NaN(s);
    return;
end

% Rounds 2 and 3
[Q, R2] = IntraOrtho(Q, 'cholqr_ro');

% Combine
R = R2*R1;

% fprintf('%s:\n',mfilename);
% fprintf('||I - Q''*Q|| = %0.5e\n', norm(eye(s) - Q'*Q));
% fprintf('||Q*R - X|| = %0.5e\n', norm(Q*R - X));
end
function [Q, R] = cholqr(X)
% [Q, R] = CHOLQR(X) computes Cholesky QR factorization of the m x s matrix
% X as described in [Stathopoulos & Wu 2002].

%%
[m, s] = size(X);
A = X'*X;
[~, p] = chol(A);
if ~p
    R = chol(A);
    Q = X/R;
else
    Q = NaN(m,s);
    R = NaN(s);
end

% fprintf('%s:\n',mfilename);
% fprintf('||I - Q''*Q|| = %0.5e\n', norm(eye(s) - Q'*Q));
% fprintf('||Q*R - X|| = %0.5e\n', norm(Q*R - X));
end
function [Q, R, T] = mgs_svl_trans(X)
% [Q, R, T] = MGS_SVL_TRANS(X) performs an unstable Modified Gram-Schmidt
% with Schreiber-Van-Loan reformulation on the m x s matrix X.

%%
% Pre-allocate memory for Q, R, and T
[m, s] = size(X);
R = zeros(s,s);
Q = zeros(m,s);
T = eye(s,s);

R(1,1) = norm(X(:,1));
Q(:,1) = X(:,1) / R(1,1);
for k = 1:s-1
    w = X(:,k+1);
    R(1:k,k+1) = T(1:k,1:k) * (Q(:,1:k)' * w);
    
    w = w - Q(:,1:k) * R(1:k,k+1);
    R(k+1,k+1) = norm(w);
    Q(:,k+1) = w / R(k+1,k+1);
    
    T(1:k,k+1) = -T(1:k,1:k) * (Q(:,1:k)' * Q(:,k+1));
end

% fprintf('%s:\n',mfilename);
% fprintf('||I - Q''*Q|| = %0.5e\n', norm(eye(s) - Q'*Q));
% fprintf('||Q*R - X|| = %0.5e\n', norm(Q*R - X));
end
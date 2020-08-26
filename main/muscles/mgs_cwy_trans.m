function [Q, R, T] = mgs_cwy_trans(X)
% [Q, R, T] = MGS_CWY_TRANS(X) performs an unstable Modified Gram-Schmidt
% with Compact WY form on the m x s matrix X.

%%
% Pre-allocate memory for Q, R, and T
[m, s] = size(X);
Q = zeros(m,s);
R = zeros(s,s);
T = eye(s,s);

u = X(:,1);
for k = 1:s-1
    w = X(:,k+1);
    % Only synchronization point-- block inner product allows us to
    % represent what is normally 3 inner products all at once!
    if k == 1
        tmp = u' * [u w];
        R(k,k) = sqrt(tmp(k,1));
    elseif k > 1
        tmp = [Q(:,1:k-1) u]' * [u w];
        R(k,k) = sqrt(tmp(k,1));
        T(1:k-1,k) = -T(1:k-1,1:k-1) * (tmp(1:k-1,1) / R(k,k));
    end
    
    tmp(k,2) = tmp(k,2) / R(k,k);
    R(1:k,k+1) = T(1:k,1:k) * tmp(:,2);
    Q(:,k) = u / R(k,k);
    u = w - Q(:,1:k) * R(1:k,k+1);
end
R(s,s) = norm(u);
Q(:,s) = u / R(s,s);

% fprintf('%s:\n',mfilename);
% fprintf('||I - Q''*Q|| = %0.5e\n', norm(eye(s) - Q'*Q));
% fprintf('||Q*R - X|| = %0.5e\n', norm(Q*R - X));
end

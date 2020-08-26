function [Q, R, T] = mgs_cwy(X)
% [Q, R, T] = MGS_CWY(X) performs Modified Gram-Schmidt with Compact WY
% form on the m x s matrix X. Same as Algorithm 6 from [Swirydowicz et. al.
% 2020].

%%
% Pre-allocate memory for Q, R, and T
[m, s] = size(X);
Q = zeros(m,s);
R = zeros(s,s);
T = eye(s,s);

q_tmp = X(:,1);
for k = 1:s-1
    % Pull out vector (keeps MATLAB from copying full X repeatedly)
    xk = X(:,k+1);
    
    % Compute temporary quantities -- the only sync point!
    if k == 1
        tmp = q_tmp' * [q_tmp xk];
        r_diag = sqrt(tmp(k,1));
    else
        tmp = [Q(:,1:k-1) q_tmp]' * [q_tmp xk];
        r_diag = sqrt(tmp(k,1));
        T(1:k-1,k) = -T(1:k-1,1:k-1) * (tmp(1:k-1,1) / r_diag);
    end
    
    R(k,k) = r_diag;
    tmp(k,2) = tmp(k,2) / r_diag;
    R(1:k,k+1) = T(1:k,1:k)' * tmp(:,2);
    
    Q(:,k) = q_tmp / r_diag;
    q_tmp = xk - Q(:,1:k) * R(1:k,k+1);
    
    fprintf('%d: LOO: %2.4e |', k, norm( eye(k) - Q(:, 1:k)' * Q(:, 1:k) ) );
    fprintf('  Res: %2.4e\n', norm( X(:,1:k) - Q(:,1:k) * R(1:k,1:k) ) / norm(X(:,1:k)) );
end
R(s,s) = norm(q_tmp);
Q(:,s) = q_tmp / R(s,s);

fprintf('%d: LOO: %2.4e |', k+1, norm( eye(s) - Q' * Q ) );
fprintf('  Res: %2.4e\n', norm( X - Q * R ) / norm(X) );
end

function [Q, R] = cgs_iro_ls(X)
% [Q, R] = CGS_IRO_LS(X) performs Classical Gram-Schmidt with
% Low-Synchronzation and Reorthogonalization on the m x s matrix X.
% CGS_IRO_LS is Algorithm 3 from [Swirydowicz, et. al., 2020].

% Note: as much as possible, we avoid calling entries of R or vectors of X
% or Q to keep MATLAB from copying these structures when doing arithmetic.
% While this may also speed up the code, the main reason for writing the
% routine this way is that it makes stability analysis later on more
% transparent, so that we can distinguish between temporary quantities and
% "finished products."

%%
% Pre-allocate memory for Q and R
[m,s] = size(X);
Q = zeros(m,s);
R = zeros(s,s);

q_tmp = X(:,1);
for k = 2:s
    % Pull out vector (keeps MATLAB from copying full X repeatedly)
    xk = X(:,k);
    
    % Compute temporary quantities -- the only sync point!
    if k == 2
        r_tmp = q_tmp' * [q_tmp xk];
    else
        tmp = [Q(:,1:k-2) q_tmp]' * [q_tmp xk];
        w = tmp(1:k-2,1);
        z = tmp(1:k-2,2);
        r_tmp = tmp(end,:) - w'*[w z];
    end
    
    % Pythagorean trick for R diagonals; assign finished entry
    r_diag = sqrt(r_tmp(1));
    
    % Assign finished entries of R
    R(k-1,k-1) = r_diag;
    R(k-1,k) = r_tmp(2) / r_diag;
    
    if k == 2
        % Finish normalizing Q(:,k-1)
        Q(:,k-1) = q_tmp / r_diag;
    else
        % Assign finished entries of R
        R(1:k-2,k-1) = R(1:k-2,k-1) + w;
        R(1:k-2,k) = z;
        
        % Finish normalizing Q(:,k-1)
        Q(:,k-1) = (q_tmp - Q(:,1:k-2) * w) / r_diag;
    end
    
    % Set up temporary vector for next iteration
    q_tmp = xk - Q(:,1:k-1) * R(1:k-1,k);
    
    fprintf('%d: LOO: %2.4e |', k-1, norm( eye(k-1) - Q(:, 1:k-1)' * Q(:, 1:k-1) ) );
    fprintf('  Res: %2.4e\n', norm( X(:,1:k-1) - Q(:,1:k-1) * R(1:k-1,1:k-1) ) / norm(X(:,1:k-1)) );
end

% Finish normalizing last basis vector and assign last diagonal entry of R.
%  We can do this in one sync, no norm needed.
tmp = [Q(:,1:s-1) q_tmp]' * q_tmp;
w = tmp(1:s-1,1);
r_tmp = tmp(end,1) - w'*w;
r_diag = sqrt(r_tmp);
R(s,s) = r_diag;
R(1:s-1,s) = R(1:s-1,s) + w;
Q(:,s) = (q_tmp - Q(:,1:s-1) * w) / r_diag;

fprintf('%d: LOO: %2.4e |', s, norm( eye(s) - Q' * Q ) );
fprintf('  Res: %2.4e\n', norm( X - Q * R ) / norm(X) );
end
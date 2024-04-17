function [Q, R] = cgs_iro_ls(X, verbose)
% [Q, R] = CGS_IRO_LS(X, verbose) performs Classical Gram-Schmidt with
% Low-Synchronzation and Reorthogonalization on the m x s matrix X.
% CGS_IRO_LS is Algorithm 3 from [Swirydowicz, et. al., 2020], with some
% cosmetic modifications.
%
% See INTRAORTHO for more details about the parameters.
%
% Part of [BlockStab](https://github.com/katlund) package.  Check README
% for how to properly cite and reuse this file.

%%
% Default: debugging off
if nargin < 2
    verbose = 0;
end

% Pre-allocate memory for Q and R
[m,s] = size(X);
Q = zeros(m,s);
R = zeros(s,s);

q_tmp = X(:,1);

if verbose
    fprintf('         LOO      |    RelRes\n');
    fprintf('-----------------------------------\n');
end

for k = 2:s
    % Pull out vector (keeps MATLAB from copying full X repeatedly)
    u = X(:,k);
    
    % Compute temporary quantities -- the only sync point!
    if k == 2
        r_tmp = q_tmp' * [q_tmp u];
    else
        tmp = [Q(:,1:k-2) q_tmp]' * [q_tmp u];
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
    q_tmp = u - Q(:,1:k-1) * R(1:k-1,k);
    
    if verbose
        fprintf('%3.0d:', k-1);
        fprintf('  %2.4e  |',...
            norm( eye(k-1) - Q(:,1:k-1)' * Q(:,1:k-1) ) );
        fprintf('  %2.4e\n',...
            norm( X(:,1:k-1) - Q(:,1:k-1) * R(1:k-1,1:k-1) ) / norm(X(:,1:k-1)) );
    end
end

% Finish normalizing last basis vector and assign last diagonal entry of R.
% We can do this in one sync, no norm needed.
tmp = [Q(:,1:s-1) q_tmp]' * q_tmp;
w = tmp(1:s-1,1);
r_tmp = tmp(end,1) - w'*w;
r_diag = sqrt(r_tmp);
R(s,s) = r_diag;
R(1:s-1,s) = R(1:s-1,s) + w;
Q(:,s) = (q_tmp - Q(:,1:s-1) * w) / r_diag;

if verbose
    fprintf('%3.0d:', k);
    fprintf('  %2.4e  |', norm( eye(s) - Q' * Q ) );
    fprintf('  %2.4e\n', norm( X - Q * R ) / norm(X) );
end
end
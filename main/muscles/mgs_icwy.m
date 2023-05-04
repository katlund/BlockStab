function [Q, R, T] = mgs_icwy(X, verbose)
% [Q, R, T] = MGS_ICWY(X, verbose) performs Modified Gram-Schmidt with
% Inverse Compact WY form on the m x s matrix X. Same as Algorithm 5 from
% [Swirydowicz et. al. 2020].
%
% See INTRAORTHO for more details about the parameters.
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
% Default: debugging off
if nargin < 2
    verbose = 0;
end

% Pre-allocate memory for Q, R, and T
[m, s] = size(X);
Q = zeros(m,s);
R = zeros(s,s);
T = eye(s,s);

q_tmp = X(:,1);

if verbose
    fprintf('         LOO      |    RelRes\n');
    fprintf('-----------------------------------\n');
end

for k = 1:s-1
    % Pull out vector (keeps MATLAB from copying full X repeatedly)
    u = X(:,k+1);
    
    % Compute temporary quantities -- the only sync point!
    if k == 1
        tmp = q_tmp' * [q_tmp u];
        r_diag = sqrt(tmp(k,1));
    else
        tmp = [Q(:,1:k-1) q_tmp]' * [q_tmp u];
        r_diag = sqrt(tmp(k,1));
        T(1:k-1, k) = tmp(1:k-1,1) / r_diag;
    end
    
    R(k,k) = r_diag;
    tmp(k,2) = tmp(k,2) / r_diag;
    R(1:k,k+1) = T(1:k,1:k)' \ tmp(:,2);
    
    Q(:,k) = q_tmp / r_diag;
    q_tmp = u - Q(:,1:k) * R(1:k,k+1);
    
    if verbose
        fprintf('%3.0d:', k);
        fprintf('  %2.4e  |',...
            norm( eye(k) - Q(:,1:k)' * Q(:,1:k) ) );
        fprintf('  %2.4e\n',...
            norm( X(:,1:k) - Q(:,1:k) * R(1:k,1:k) ) / norm(X(:,1:k)) );
    end
end
R(s,s) = norm(q_tmp);
Q(:,s) = q_tmp / R(s,s);

if verbose
    fprintf('%3.0d:', k+1);
    fprintf('  %2.4e  |', norm( eye(s) - Q' * Q ) );
    fprintf('  %2.4e\n', norm( X - Q * R ) / norm(X) );
end

% Barlow check
% for k = 1:s
%    fprintf('Gamma: %2.4e\n', norm( (eye(k) - inv(T(1:k,1:k)) ) * R(1:k,1:k) )); 
% end
end

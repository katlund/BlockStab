function [Q, R] = cgs_iro(X, verbose)
% [Q, R] = CGS_IRO(X, verbose) performs Classical Gram-Schmidt with Inner
% ReOrthogonalization on the m x s matrix X.
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

% Pre-allocate memory for Q and R
[m, s] = size(X);
R = zeros(s);
Q = zeros(m,s);

R(1,1) = norm(X(:,1));
Q(:,1) = X(:,1) / R(1,1);

if verbose
    fprintf('         LOO      |    RelRes\n');
    fprintf('-----------------------------------\n');
    fprintf('%3.0d:', 1);
    fprintf('  %2.4e  |',...
        norm( 1 - Q(:,1)' * Q(:,1) ) );
    fprintf('  %2.4e\n',...
        norm( X(:,1) - Q(:,1) * R(1,1) ) / norm(X(:,1)) );
end

for k = 1:s-1
    w = X(:,k+1);
    
    % First CGS step
    R1 = Q(:,1:k)' * w;
    w = w - Q(:,1:k) * R1;
    r1 = norm(w);
    w = w / r1;
    
    % Second CGS step
    R(1:k,k+1) = Q(:,1:k)' * w;
    w = w - Q(:,1:k) * R(1:k,k+1);
    R(k+1,k+1) = norm(w);
    Q(:,k+1) = w / R(k+1,k+1);
    
    % Combine both steps
    R(1:k,k+1) = R1 + R(1:k,k+1) * r1;
    R(k+1,k+1) = R(k+1,k+1) * r1;
    
    if verbose
        fprintf('%3.0d:', k+1);
        fprintf('  %2.4e  |',...
            norm( eye(k+1) - Q(:,1:k+1)' * Q(:,1:k+1) ) );
        fprintf('  %2.4e\n',...
            norm( X(:,1:k+1) - Q(:,1:k+1) * R(1:k+1,1:k+1) ) / norm(X(:,1:k+1)) );
    end
end
end
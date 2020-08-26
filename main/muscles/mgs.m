function [Q, R] = mgs(X, verbose)
% [Q, R] = MGS(X, verbose) performs Modified Gram-Schmidt on the m x s
% matrix X.
%
% See INTRAORTHO for more details about the parameters.

%%
% Default: debugging off
if nargin < 2
    verbose = 0;
end

% Pre-allocate memory for Q and R
[m, s] = size(X);
R = zeros(s,s);
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
    for j = 1:k
        R(j,k+1) = Q(:,j)' * w;
        w = w - Q(:,j) * R(j,k+1);
    end
    R(k+1,k+1) = norm(w);
    Q(:,k+1) = w / R(k+1,k+1);
    
    if verbose
        fprintf('%3.0d:', k+1);
        fprintf('  %2.4e  |',...
            norm( eye(k+1) - Q(:,1:k+1)' * Q(:,1:k+1) ) );
        fprintf('  %2.4e\n',...
            norm( X(:,1:k+1) - Q(:,1:k+1) * R(1:k+1,1:k+1) ) / norm(X(:,1:k+1)) );
    end
end
end
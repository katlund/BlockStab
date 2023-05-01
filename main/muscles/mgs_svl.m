function [Q, R, T] = mgs_svl(X, verbose)
% [Q, R, T] = MGS_SVL(X, verbose) performs Modified Gram-Schmidt with
% Schreiber-Van-Loan reformulation on the m x s matrix X. MGS_SVL is
% equivalent to MGS2 as described in [Barlow 2019].
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
    R(1:k,k+1) = T(1:k,1:k)' * (Q(:,1:k)' * w);
    w = w - Q(:,1:k) * R(1:k,k+1);
    R(k+1,k+1) = norm(w);
    Q(:,k+1) = w / R(k+1,k+1);
    T(1:k,k+1) = -T(1:k,1:k) * (Q(:,1:k)' * Q(:,k+1));
    
    if verbose
        fprintf('%3.0d:', k+1);
        fprintf('  %2.4e  |',...
            norm( eye(k+1) - Q(:,1:k+1)' * Q(:,1:k+1) ) );
        fprintf('  %2.4e\n',...
            norm( X(:,1:k+1) - Q(:,1:k+1) * R(1:k+1,1:k+1) ) / norm(X(:,1:k+1)) );
    end
end

% Barlow check
% for k = 1:s
%    fprintf('Gamma: %2.4e\n', norm( (eye(k) - T(1:k,1:k) ) * R(1:k,1:k) ));
% end
end
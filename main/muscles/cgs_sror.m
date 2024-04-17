function [Q, R] = cgs_sror(X, rpltol, verbose)
% [Q, R] = CGS_SROR(X, rpltol, verbose) performs Classical Gram-Schmidt
% with Selective ReOrthogonalization and Replacement on the m x s matrix X
% as described in [Stewart 2008]. The core part of the routine is contained
% in CGS_STEP_SROR.
%
% See INTRAORTHO for more details about the parameters.
%
% Part of [BlockStab](https://github.com/katlund) package.  Check README
% for how to properly cite and reuse this file.

%%
% Default for rpltol
if isempty(rpltol)
    rpltol = 1;
end

% Default: debugging off
if nargin < 3
    verbose = 0;
end

% Pre-allocate memory for Q and R
[m, s] = size(X);
Q = zeros(m, s);
R = zeros(s);

[Q(:,1), ~, R(1,1)] = cgs_step_sror(zeros(m,0), X(:,1), rpltol);

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
    [Q(:,k+1), R(1:k,k+1), R(k+1,k+1)] =...
        cgs_step_sror(Q(:,1:k), X(:,k+1), rpltol);
    
    if verbose
        fprintf('%3.0d:', k+1);
        fprintf('  %2.4e  |',...
            norm( eye(k+1) - Q(:,1:k+1)' * Q(:,1:k+1) ) );
        fprintf('  %2.4e\n',...
            norm( X(:,1:k+1) - Q(:,1:k+1) * R(1:k+1,1:k+1) ) / norm(X(:,1:k+1)) );
    end
end
end
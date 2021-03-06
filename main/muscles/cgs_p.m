function [Q, R] = cgs_p(X, verbose)
% [Q, R] = CGS_P(X, verbose) performs Classical Gram-Schmidt with
% Pythagorean modification (from [Smoktunowicz et al., 2006]) on the m x s
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
R = zeros(s);
Q = zeros(m,s);

w = X(:,1);
r_diag = norm(w);
Q(:,1) = w / r_diag;
R(1,1) = r_diag;

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
    xk = X(:,k+1);
    r_vec = Q(:,1:k)' * xk;
    w = xk - Q(:,1:k) * r_vec;
    R(1:k,k+1) = r_vec;
    
    psi = norm(xk);
    phi = norm(r_vec);
    r_diag = sqrt(psi - phi) * sqrt(psi + phi);
    Q(:,k+1) = w / r_diag;
    R(k+1,k+1) = r_diag;
    
    if verbose
        fprintf('%3.0d:', k+1);
        fprintf('  %2.4e  |',...
            norm( eye(k+1) - Q(:,1:k+1)' * Q(:,1:k+1) ) );
        fprintf('  %2.4e\n',...
            norm( X(:,1:k+1) - Q(:,1:k+1) * R(1:k+1,1:k+1) ) / norm(X(:,1:k+1)) );
    end
end
end
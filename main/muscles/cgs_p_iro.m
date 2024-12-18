function [Q, R] = cgs_p_iro(X, verbose)
% [Q, R] = CGS_P_IRO(X, verbose) performs Classical Gram-Schmidt with
% Pythagorean correction and Inner ReOrthogonalization on the m x s matrix
% X.  See CGS_P for the base version.
%
% See INTRAORTHO for more details about the parameters.
%
% Part of [BlockStab](https://github.com/katlund/BlockStab) package.  Check README
% for how to properly cite and reuse this file.

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
R(1,1) = norm(w);
Q(:,1) = w / R(1,1);

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
    % Set up next vector
    w = X(:,k+1);

    % First CGS-P step
    tmp = [Q(:,1:k) w]' * w;
    R1 = tmp(1:k);
    w = w - Q(:,1:k) * R1;
    
    psi = sqrt(tmp(k+1));
    phi = norm(R1);
    r1 = sqrt(psi - phi) * sqrt(psi + phi);
    w = w / r1;
    
    % Second CGS-P step
    tmp = [Q(:,1:k) w]' * w;
    R(1:k,k+1) = tmp(1:k);
    w = w - Q(:,1:k) * R(1:k,k+1);
    
    psi = sqrt(tmp(k+1));
    phi = norm(R(1:k,k+1));
    R(k+1,k+1) = sqrt(psi - phi) * sqrt(psi + phi);
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
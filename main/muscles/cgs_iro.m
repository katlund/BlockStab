function [Q, R] = cgs_iro(X)
% [Q, R] = CGS_IRO(X) performs Classical Gram-Schmidt with Inner
% ReOrthogonalization on the m x s matrix X.

%%
% Pre-allocate memory for Q and R
[m, s] = size(X);
R = zeros(s);
Q = zeros(m,s);

% Classical Gram-Schmidt with in situ reorthognalization
R(1,1) = norm(X(:,1));
Q(:,1) = X(:,1)/R(1,1);
for k = 1:s-1
    w = X(:,k+1);
    
    % First CGS step
    R1 = Q(:,1:k)'*w;
    w = w - Q(:,1:k)*R1;
    r1 = norm(w);
    w = w/r1;
    
    % Second CGS step
    R(1:k,k+1) = Q(:,1:k)'*w;
    w = w - Q(:,1:k)*R(1:k,k+1);
    R(k+1,k+1) = norm(w);
    Q(:,k+1) = w/R(k+1,k+1);
    
    % Combine both steps
    R(1:k,k+1) = R1 + R(1:k,k+1) * r1;
    R(k+1,k+1) = R(k+1,k+1) * r1;
end

% fprintf('%s:\n',mfilename);
% fprintf('||I - Q''*Q|| = %0.5e\n', norm(eye(s) - Q'*Q));
% fprintf('||Q*R - X|| = %0.5e\n', norm(Q*R - X));
end
function [Q, R] = cgs(X)
% [Q, R] = CGS(X) performs Classical Gram-Schmidt on the m x s matrix X

%%
% Pre-allocate memory for Q and R
[m, s] = size(X);
R = zeros(s);
Q = zeros(m,s);

% Classical Gram-Schmidt
w = X(:,1);
R(1,1) = norm(w);
Q(:,1) = w/R(1,1);
for k = 1:s-1
    w = X(:,k+1);
    R(1:k,k+1) = Q(:,1:k)'*w;
    w = w - Q(:,1:k)*R(1:k,k+1);
    R(k+1,k+1) = norm(w);
    Q(:,k+1) = w/R(k+1,k+1);
end

% fprintf('%s:\n',mfilename);
% fprintf('||I - Q''*Q|| = %0.5e\n', norm(eye(s) - Q'*Q));
% fprintf('||Q*R - X|| = %0.5e\n', norm(Q*R - X));
end
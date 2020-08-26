function [Q, R] = mgs(X)
% [Q, R] = MGS(X) performs Modified Gram-Schmidt on the m x s matrix X.

%%
% Pre-allocate memory for Q and R
[m, s] = size(X);
R = zeros(s,s);
Q = zeros(m,s);

% Modified Gram-Schmidt
R(1,1) = norm(X(:,1));
Q(:,1) = X(:,1)/R(1,1);
for k = 1:s-1
    w = X(:,k+1);
    for j = 1:k
        R(j,k+1) = Q(:,j)'*w;
        w = w - Q(:,j)*R(j,k+1);
    end
    R(k+1,k+1) = norm(w);
    Q(:,k+1) = w/R(k+1,k+1);
end

% fprintf('%s:\n',mfilename);
% fprintf('||I - Q''*Q|| = %0.5e\n', norm(eye(s) - Q'*Q));
% fprintf('||Q*R - X|| = %0.5e\n', norm(Q*R - X));
end
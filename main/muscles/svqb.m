function [Q, R] = svqb(X)
% [Q, R] = SVQB(X) computes the QR factorization of the m x s matrix X as
% described in [Stathopoulos & Wu 2002].

%%
[m,s] = size(X);

A = X'*X;
B = diag(1./sqrt(diag(A)));
C = B*A*B;

if sum(isnan(C),'all') > 0
    Q = NaN(m,s);
    R = NaN(s);
else
    [V, D] = eig(C, 'vector');
    tau = eps*max(D);
    for i = 1:s
        if D(i) < tau
            D(i) = tau;
        end
    end
    R = (diag(sqrt(D))/V)*diag(sqrt(diag(A)));
    Q = X*(B*(V*diag(1./sqrt(D))));
end

% fprintf('%s:\n',mfilename);
% fprintf('||I - Q''*Q|| = %0.5e\n', norm(eye(s) - Q'*Q));
% fprintf('||Q*R - X|| = %0.5e\n', norm(Q*R - X));
end
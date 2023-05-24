function [Q, R] = iter_cholqr(X)
% [Q, R] = ITER_CHOLQR(X) computes the Iterated Cholesky QR factorization
% of the m x s matrix X.  This algorithm is equivalent to Algorithm 4.1 of
% [Fukaya, et. al. 2018].
%
% The number of iterations is limited to 10.
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
[m, s] = size(X);
Q = X;
I = eye(s);
R = I;
normX2 = norm(X,2)^2;

iter = 0;
loss_ortho = norm(Q' * Q - I, 'fro');
if isnan(loss_ortho)
    Q = NaN(m,s);
    R = NaN(s,s);
    fprintf('%s failed to converge.\n', mfilename);
    return
end
while (abs(loss_ortho - sqrt(s)*eps) > eps) && (iter < 10)
    iter = iter + 1;
    A = Q' * Q;
    [R2, flag] = chol(A);
    if flag ~= 0 % need to shift
        sh = 11 * (m*s + s * (s+1)) *eps * normX2;
        R2 = chol(A + sh*I);
    end
    Q = Q/R2;
    R = R2*R;
    loss_ortho = norm(Q' * Q - I, 'fro'); % sync point!
    if loss_ortho >= 1
        % then there's probably a zero column and the difference will never
        % be less than 1
        Q = NaN(m,s);
        R = NaN(s,s);
        fprintf('%s failed to converge.\n', mfilename);
        return
    end
end
if (abs(loss_ortho - sqrt(s)*eps) <= eps)
    fprintf('\t%s converged in %d iterations.\n', mfilename, iter);
else
    fprintf('%s failed to converge in 10 iterations. Solution may not be reliable.\n', mfilename);
end
end
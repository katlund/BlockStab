function R = chol_free_vpa(A)
% R = CHOL_FREE_VPA(A) computes the Cholesky factor of A according to
% Algorithm 10.2 from [Higham 2002], using simulated quadruple precision
% via MATLAB vpa. Unlike Matlab's built-in CHOL, there is no fail-safe for
% numerical stability, so that we may more freely observe the algorithm's
% behavior for matrices that are not numerically positive definite.

%%
A = vpa(A,32);
s = size(A,1);
R = vpa(zeros(s,s),32);
for j = 1:s
    for i = 1:j-1
        rij = A(i,j) - R(1:i-1,i)' * R(1:i-1,j);
        R(i,j) = rij / R(i,i);
    end
    R(j,j) = sqrt( A(j,j) - R(1:j-1,j)' * R(1:j-1,j) );
end
end
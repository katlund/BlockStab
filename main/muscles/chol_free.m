function R = chol_free(A)
% R = CHOL_FREE(A) computes the Cholesky factor of A according to Algorithm
% 10.2 from [Higham 2002].  Unlike Matlab's built-in CHOL, there is no
% fail-safe for numerical stability, so that we may more freely observe the
% algorithm's behavior for matrices that are not numerically positive
% definite.

%%
s = size(A,1);
R = zeros(s,s);
for j = 1:s
    for i = 1:j-1
        rij = A(i,j) - R(1:i-1,i)' * R(1:i-1,j);
        R(i,j) = rij / R(i,i);
    end
    R(j,j) = sqrt( A(j,j) - R(1:j-1,j)' * R(1:j-1,j) );
end
end
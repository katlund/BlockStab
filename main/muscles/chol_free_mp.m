function R = chol_free_mp(A)
% R = CHOL_FREE_MP(A) computes the Cholesky factor of A according to Algorithm
% 10.2 from [Higham 2002], using simulated quadruple precision via the Advanpix 
% Multiprecision Computing Toolbox.  Unlike Matlab's built-in CHOL, there is no
% fail-safe for numerical stability, so that we may more freely observe the
% algorithm's behavior for matrices that are not numerically positive
% definite.

%%
A=mp(A,34);
s = size(A,1);
R = mp(zeros(s,s),34);
for j = 1:s
    for i = 1:j-1
        rij = A(i,j) - R(1:i-1,i)' * R(1:i-1,j);
        R(i,j) = rij / R(i,i);
    end
    R(j,j) = sqrt( A(j,j) - R(1:j-1,j)' * R(1:j-1,j) );
end
end
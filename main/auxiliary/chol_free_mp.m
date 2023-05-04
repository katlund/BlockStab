function R = chol_free_mp(A, param)
% R = CHOL_FREE_MP(A, param) computes the Cholesky factor of A according to
% Algorithm 10.2 from [Higham 2002].  Unlike Matlab's built-in CHOL, there
% is no fail-safe for numerical stability, so that we may more freely
% observe the algorithm's behavior for matrices that are not numerically
% positive definite.
%
% This mixed precision version furthermore selectively uses simulated quad
% precision via Advanpix or the Symbolic Math Toolbox as specified by
% param.mp, with returned quantities cast to double. See MP_SWITCH for more
% details about the param struct.
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
% Set up quad-precision subroutine
switch param.mp_package
    case 'advanpix'
        param.mp_digits = 34;
    case 'symbolic math'
        param.mp_digits = 32;
end
qp = @(x) mp_switch(x, param);

A = qp(A);
s = size(A,1);
R = qp(zeros(s,s));

% CHOL_FREE
for j = 1:s
    for i = 1:j-1
        rij = A(i,j) - R(1:i-1,i)' * R(1:i-1,j);
        R(i,j) = rij / R(i,i);
    end
    R(j,j) = sqrt( A(j,j) - R(1:j-1,j)' * R(1:j-1,j) );
end
end
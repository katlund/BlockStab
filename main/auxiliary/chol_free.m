function R = chol_free(A, param)
% R = CHOL_FREE(A, param) computes the Cholesky factor of A according to
% Algorithm 10.2 from [Higham 2002].  Unlike Matlab's built-in CHOL, there
% is no fail-safe for numerical stability, so that we may more freely
% observe the algorithm's behavior for matrices that are not numerically
% positive definite.
%
% When no param struct is provided, or when it is provided as an empty set,
% CHOL_FREE is run in standard double precision.
%
% To specify a multiprecision implementation, param with the following
% fields must be provided:
% - .mp_package: 'advanpix', 'symbolic math', or 'none'
% - .mp_digits: the desired number of digits, according to the chosen
%   package specifications.  The default for 'advanpix' is 34 (quad) and
%   for 'symbolic math' is 32 (quad).
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
% Defaults
if nargin == 1
    qp = @(x) x;
elseif nargin == 2
    if isempty(param)
        qp = @(x) x;
    else
        qp = @(x) mp_switch(x, param);
    end
end

% Allocate space for R
A = qp(A);
s = size(A,1);
R = qp(zeros(s,s));

for j = 1:s
    for i = 1:j-1
        R(i,j) = ( A(i,j) - R(1:i-1,i)' * R(1:i-1,j) ) / R(i,i);
    end
    R(j,j) = sqrt( A(j,j) - R(1:j-1,j)' * R(1:j-1,j) );
end
end
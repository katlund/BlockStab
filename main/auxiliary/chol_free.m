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
% - .mp_pair: a cell pair of precisions, with the second being the higher
% precision (and what this routine will use to compute Cholesky)
%
% Part of [BlockStab](https://github.com/katlund) package.  Check README
% for how to properly cite and reuse this file.

%%
% Defaults
if nargin == 1
    p2 = @(x) x;
elseif nargin == 2
    if isempty(param)
        p2 = @(x) x;
    elseif ~isfield(param, 'mp_package')
        p2 = @(x) x;
    else
        p2 = @(x) mp_switch(x, param.mp_package, param.mp_pair{2});
    end
end

% Allocate space for R
A = p2(A);
s = size(A,1);
R = p2(zeros(s,s));

for j = 1:s
    for i = 1:j-1
        R(i,j) = ( A(i,j) - R(1:i-1,i)' * R(1:i-1,j) ) / R(i,i);
    end
    R(j,j) = sqrt( A(j,j) - R(1:j-1,j)' * R(1:j-1,j) );
end
end
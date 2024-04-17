function [Q, R] = cholqr_ro(X, param)
% [Q, R] = CHOLQR_RO(X, param) computes Cholesky QR factorization with
% (outer) ReOrthogonalization of the m x s matrix X.  The default Cholesky
% subroutine is CHOL_NAN; the other option is CHOL_FREE, which can be
% further forced to operate in multiprecision with the following fields
% specified in the param struct:
% - .chol = 'chol_nan' or 'chol_free'
% - .mp_package = 'advanpix' or 'symbolic toolbox'
% - .mp_pair: a cell pair of precisions, with the second being the higher
%    precision
%
% Part of [BlockStab](https://github.com/katlund) package.  Check README
% for how to properly cite and reuse this file.

%%
% Default
if nargin == 1
    param.chol = 'chol_nan';
end

[Q1, R1] = cholqr(X, param);
[Q, R] = cholqr(Q1, param);
R = R * R1;
end
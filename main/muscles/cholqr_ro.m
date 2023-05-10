function [Q, R] = cholqr_ro(X, param)
% [Q, R] = CHOLQR_RO(X, param) computes Cholesky QR factorization with
% (outer) ReOrthogonalization of the m x s matrix X.  The default Cholesky
% subroutine is CHOL_NAN; the other option is CHOL_FREE, which can be
% further forced to operate in mixed precision with the following fields
% specified in the param struct:
% - .chol = 'chol_nan' or 'chol_free'
% - .mp_package = 'advanpix' or 'symbolic toolbox'
% - .mp_digits = desired number of digits (toolbox-dependent)
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
% Default
if nargin == 1
    param.chol = 'chol_nan';
end

[Q1, R1] = cholqr(X, param);
[Q, R] = cholqr(Q1, param);
R = R*R1;

end
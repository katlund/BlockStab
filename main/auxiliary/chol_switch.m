function [R, nan_flag] = chol_switch(A, param)
% [R, nan_flag] = CHOL_SWITCH(A, param) switches which Cholesky subroutine
% is called based on param.chol.  The default is CHOL_NAN, which returns
% NaN whenever CHOL throws an error.  Other options include the following:
% - CHOL, the built-in MATLAB routine, which will throw errors when A is
% not numerically positive definite
% - CHOL_FREE, which is not optimized for performance but will not throw
% errors or return NaN when A is not numerically positive definite; a param
% struct with the following fields must be specified to handle arbitrary
% precision:
% - .mp_package: 'advanpix', 'symbolic math', or 'none'
% - .mp_pair: a cell apir of precisions
% When these fields are left undefined, standard double precision is used.
%
% Note: both CHOL and CHOL_NAN are overloaded to operate in the precision
% of the provided argument, i.e., if A is already stored in quadruple
% precision, then the Cholesky factor will be returned in that same
% precision.  CHOL_FREE, however, will only operate in double unless param
% is provided.

%%
% Defaults
if nargin == 1
    param.chol = 'chol_nan';
elseif nargin ==2
    if ~isfield(param, 'chol')
        param.chol = 'chol_nan';
    else
        if isempty(param.chol)
            param.chol = 'chol_nan';
        end
    end
end
if isempty(param)
    param.chol = 'chol_nan';
end

% Switch
switch lower(param.chol)
    case 'chol'
        R = chol(A);

    case 'chol_nan'
        [R, nan_flag] = chol_nan(A);

    case 'chol_free'
        R = chol_free(A, param);
        nan_flag = NaN;
end
end
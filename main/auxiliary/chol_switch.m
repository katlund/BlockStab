function [R, nan_flag] = chol_switch(A, param)
% [R, nan_flag] = CHOL_SWITCH(A, param) switches which Cholesky subroutine
% is called based on param.chol.  The default is CHOL_NAN.

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

% Switch
switch lower(param.chol)
    case 'chol_nan'
        [R, nan_flag] = chol_nan(A);

    case 'chol_free'
        R = chol_free(A);
        nan_flag = NaN;

    case 'chol_free_mp'
        R = chol_free_mp(A);
        nan_flag = NaN;
end
end
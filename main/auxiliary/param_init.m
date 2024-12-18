function param = param_init(param)
% Initialize missing parameters for BGS and INTRAORTHO.
%
% Part of [BlockStab](https://github.com/katlund/BlockStab) package.  Check README
% for how to properly cite and reuse this file.

%%
if nargin == 0
    param.chol = 'chol_nan';
    param.rpltol = 1;
    param.verbose = 0;
elseif nargin == 1
    if isempty(param)
        param.chol = 'chol_nan';
        param.rpltol = 1;
        param.verbose = 0;
    else
        if ~isfield(param, 'chol')
            param.chol = 'chol_nan';
        end
        if ~isfield(param, 'rpltol')
            param.rpltol = 1;
        end
        if ~isfield(param, 'verbose')
            param.verbose = 0;
        end
    end
end

end
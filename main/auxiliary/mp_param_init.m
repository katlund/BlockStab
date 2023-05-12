function param = mp_param_init(param)
% Initialize missing parameters for mixed precision skeletons

%%
if nargin == 0
    param.chol = 'chol_free';
    param.mp_package = 'advanpix';
    param.mp_digits = 34;
    param.verbose = 0;
elseif nargin == 1
    param.chol = 'chol_free';
    if ~isfield(param, 'mp_package')
        param.mp_package = 'advanpix';
    end
    if ~isfield(param, 'mp_digits')
        switch param.mp_package
            case 'advanpix'
                param.mp_digits = 34;
            case {'symbolic math', 'vpa', 'symbolic toolbox'}
                param.mp_digits = 32;
        end
    end
    if ~isfield(param, 'verbose')
        param.verbose = 0;
    end
end

end
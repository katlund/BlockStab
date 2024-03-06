function param = mp_param_init(param)
% Initialize missing parameters for multiprecision skeletons

%%
if nargin == 0
    param.chol = 'chol_free';
    param.mp_package = 'none';
    param.mp_pair = {'single', 'double'};
    param.verbose = 0;
elseif nargin == 1
    param.chol = 'chol_free';
    if ~isfield(param, 'mp_package')
        if ~isfield(param, 'mp_pair')
            % Assume user doesn't have either advanpix or symbolic
            param.mp_package = 'none';
            param.mp_pair = {'single', 'double'};
        else
            if strcmp(param.mp_pair{1}, 'quad') || strcmp(param.mp_pair{2}, 'quad')
                param.mp_package = 'advanpix';
            else
                param.mp_package = 'none';
            end
        end
    else
        if ~isfield(param, 'mp_pair')
            switch param.mp_package
                case {'advanpix', 'symbolic math', 'vpa', 'symbolic toolbox'}
                    param.mp_pair = {'double', 'quad'};
                otherwise
                    param.mp_pair = {'single', 'double'};
            end
        end
    end
    if ~isfield(param, 'verbose')
        param.verbose = 0;
    end
end

end
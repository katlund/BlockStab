function param = mp_param_init(param)
% Initialize missing parameters for multiprecision skeletons.
%
% Part of [BlockStab](https://github.com/katlund/BlockStab) package.  Check README
% for how to properly cite and reuse this file.

%%
% Check if Advanpix and/or Symbolic Math Toolbox are is installed
advanpix_installed = exist('mp', 'file') == 2;
symbolic_installed = exist('vpa', 'file') == 2;

if nargin == 0
    param.chol = 'chol_free';
    param.mp_package = 'none';
    param.mp_pair = {'single', 'double'};
    param.verbose = 0;
elseif nargin == 1
    if ~isfield(param, 'chol')
        param.chol = 'chol_free';
    end
    if ~isfield(param, 'mp_package')
        if ~isfield(param, 'mp_pair')
            % Give preference to Advanpix
            if advanpix_installed
                param.mp_package = 'advanpix';
                param.mp_pair = {'double', 'quad'};
            else
                if symbolic_installed
                    param.mp_package = 'vpa';
                    param.mp_pair = {'double', 'quad'};
                else
                    warning('No multiprecision toolbox installed.  Setting .mp_pair to {''single'', ''double''} for MP algorithms.')
                    param.mp_package = 'none';
                    param.mp_pair = {'single', 'double'};
                end
            end
        else
            % Switch based on what's installed
            if strcmp(param.mp_pair{1}, 'quad') || strcmp(param.mp_pair{2}, 'quad')
                % Given preference to Advanpix
                if advanpix_installed
                    param.mp_package = 'advanpix';
                else
                    if symbolic_installed
                        param.mp_package = 'vpa';
                    else
                        warning('No multiprecision toolbox installed.  Setting .mp_pair to {''single'', ''double''} for MP algorithms.')
                        param.mp_package = 'none';
                        param.mp_pair = {'single', 'double'};
                    end
                end
            else
                param.mp_package = 'none';
            end
        end
    else
        % Check that mp_package matches what's installed
        switch param.mp_package
            case 'advanpix'
                if ~advanpix_installed
                    if symbolic_installed
                        warning('Advanpix not installed. Switching to Symbolic Math Toolbox for MP algorithms.')
                        param.mp_package = 'vpa';
                    else
                        warning('No multiprecision toolbox installed.')
                        param.mp_package = 'none';
                    end
                end
            case {'symbolic math', 'vpa', 'symbolic toolbox'}
                if ~symbolic_installed
                    if advanpix_installed
                        warning('Symbolic Math Toolbox not installed. Switching to Advanpix for MP algorithms.')
                        param.mp_package = 'advanpix';
                    else
                        warning('No multiprecision toolbox installed.')
                        param.mp_package = 'none';
                    end
                end
        end

        % Then deal with missing mp_pair
        if ~isfield(param, 'mp_pair')
            switch param.mp_package
                case {'advanpix', 'symbolic math', 'vpa', 'symbolic toolbox'}
                    warning('.mp_pair not specified.  Set to {''double'', ''quad''} for MP algorithms.')
                    param.mp_pair = {'double', 'quad'};

                otherwise
                    warning('No viable .mp_package specified.  Setting .mp_pair to {''single'', ''double''} for MP algorithms.')
                    param.mp_pair = {'single', 'double'};
            end
        else
            % Switch to default when quad requested but neither MP package installed
            if (strcmp(param.mp_pair{1}, 'quad') || strcmp(param.mp_pair{2}, 'quad')) ...
                    && ~(advanpix_installed || symbolic_installed)
                warning('No multiprecision toolbox installed.  Setting .mp_pair to {''single'', ''double''} for MP algorithms.')
                param.mp_package = 'none';
                param.mp_pair = {'single', 'double'};
            end
        end
    end
    if ~isfield(param, 'verbose')
        param.verbose = 0;
    end
end

end
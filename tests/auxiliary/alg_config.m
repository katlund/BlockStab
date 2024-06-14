function [skel, musc, param] = alg_config(config_file)
% [skel, musc, param] = ALG_CONFIG(config_file) configures cells skel,
% musc, and param for RUNKAPPAPLOT.
%
% Part of [BlockStab](https://github.com/katlund) package.  Check README
% for how to properly cite and reuse this file.

%%
% Convert .json to struct
alg_struct = jsondecode(fileread(config_file));

% Extract top-level fieldnames
alg_opts = fieldnames(alg_struct);
n_alg_opts = length(alg_opts);

% Extract all muscles and skeletons
musc_list = {dir('../main/muscles/').name};
musc_list{end+1} = 'houseqr.m';
try
    musc_list(1:2) = [];
catch
    error(['alg_config tries to list all muscle and skeleton routines ' ...
        'using ''dir''. Please ensure you are calling alg_config from the ' ...
        '''tests'' folder (i.e., from ''RunKappaPlot'' in the ''tests'' ' ...
        'folder) otherwise this operation will not work.'])
end
skel_list = {dir('../main/skeletons/').name};
skel_list(1:2) = [];

% Construct algorithm configurations
skel = [];
musc = [];
param = [];
for i = 1:n_alg_opts
    %% Set up muscle as a main solver
    if any(strcmp(musc_list, [alg_opts{i} '.m']))
        if isempty(alg_struct.(alg_opts{i}))
            % Then set defaults for this muscle
            musc{end+1} = alg_opts{i}; %#ok<*AGROW>
            param{end+1} = [];
            skel{end+1} = [];

        else
            % Extract param_opts
            param_opts = alg_struct.(alg_opts{i}).param;
            n_param_opts = length(param_opts);

            % Prescribe everything
            if n_param_opts > 0
                for j = 1:n_param_opts
                    musc{end+1} = alg_opts{i};
                    param{end+1} = param_opts(j);
                    skel{end+1} = [];
                end
            else
                musc{end+1} = alg_opts{i};
                param{end+1} = [];
                skel{end+1} = [];
            end
        end

    %% Set up skeleton-muscle as a main solver
    elseif any(strcmp(skel_list, [alg_opts{i} '.m']))

        if isempty(alg_struct.(alg_opts{i}))
            % Assign skeleton and set defaults for everything else
            skel{end+1} = alg_opts{i};
            param{end+1} = [];
            musc{end+1} = [];

        else
            % Identify all parameter and muscle options
            musc_opts = fieldnames(alg_struct.(alg_opts{i}));
            ind = strcmp(musc_opts, 'param');

            if any(ind)
                % Remove param from musc_opts
                musc_opts(ind) = [];

                % Extract param_opts
                param_opts = alg_struct.(alg_opts{i}).param;
                n_param_opts = length(param_opts);
            else
                param_opts = [];
                n_param_opts = 1;
            end

            % Cycle through all possible combinations
            for j = 1:n_param_opts
                % Prescribe muscles
                switch alg_opts{i}
                    case 'bcgs_sror'
                        % Prescribe skeleton
                        skel{end+1} = alg_opts{i};
    
                        % Prescribe parameters
                        if isempty(param_opts)
                            param{end+1} = [];
                        else
                            param{end+1} = param_opts(j);
                        end
                        
                        % Muscle is fixed as cgs_sror
                        musc{end+1} = 'cgs_sror';
                        
                    case {'bcgs_iro_ls', 'bcgs_iro_ls_mp'}
                        % Prescribe skeleton
                        skel{end+1} = alg_opts{i};
    
                        % Prescribe parameters
                        if isempty(param_opts)
                            param{end+1} = [];
                        else
                            param{end+1} = param_opts(j);
                        end

                        % All muscles ignored
                        musc{end+1} = [];
                        
                    otherwise
                        n_musc_opts = length(musc_opts);
                        
                        if n_musc_opts > 0
                            for k = 1:n_musc_opts
                                % Vet remaining muscle options.  Keep only
                                % ones that match a built-in or that have
                                % an 'io_a' field.
                                if ~any(strcmp(musc_list, [musc_opts{k} '.m']))
                                    if isstruct(alg_struct.(alg_opts{i}).(musc_opts{k}))
                                        % Then configure a multi-IO
                                        skel{end+1} = alg_opts{i};
                                        musc{end+1} = alg_struct.(alg_opts{i}).(musc_opts{k});
                                        musc{end}.id = musc_opts{k};

                                        % Prescribe skeleton-level parameters
                                        if isempty(param_opts)
                                            param{end+1} = [];
                                        else
                                            param{end+1} = param_opts(j);
                                        end
                                        
                                    else
                                        skel{end+1} = [];
                                        musc{end+1} = [];
                                        param{end+1} = [];
                                        fprintf(['%s: ''%s'' is neither a built-in muscle nor a multi-IO\n' ...
                                            '\tand will therefore be excluded from the configuration.\n'], ...
                                            mfilename, musc_opts{k});
                                    end
                                else
                                    % Treat like a built-in muscle
                                    if isempty(alg_struct.(alg_opts{i}).(musc_opts{k}))
                                        % No muscle-level parameters
                                        skel{end+1} = alg_opts{i};
                                        musc{end+1} = musc_opts{k};
                    
                                        % Prescribe skeleton-level parameters
                                        if isempty(param_opts)
                                            param{end+1} = [];
                                        else
                                            param{end+1} = param_opts(j);
                                        end
                                    else
                                        % Extract muscle-level parameters
                                        musc_param_opts = ...
                                            alg_struct.(alg_opts{i}).(musc_opts{k}).param;
                                        n_musc_param_opts = length(musc_param_opts);
                            
                                        for ll = 1:n_musc_param_opts
                                            % Prescribe skeleton & muscle
                                            skel{end+1} = alg_opts{i};
                                            musc{end+1} = musc_opts{k};
                        
                                            % Prescribe appended parameters
                                            if isempty(param_opts)
                                                if isempty(musc_param_opts(ll))
                                                    param{end+1} = [];
                                                else
                                                    param{end+1} = musc_param_opts(ll);
                                                end
                                            else
                                                if isempty(musc_param_opts(ll))
                                                    param{end+1} = param_opts(j);
                                                else
                                                    param{end+1} = catstruct(...
                                                    param_opts(j),...
                                                    musc_param_opts(ll));
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        else
                            % Prescribe skeleton and default muscle
                            skel{end+1} = alg_opts{i};
                            musc{end+1} = [];
        
                            % Prescribe skeleton-level parameters
                            if isempty(param_opts)
                                param{end+1} = [];
                            else
                                param{end+1} = param_opts(j);
                            end
                        end
                end
            end
        end
    else
        warning('%s is not a known muscle or skeleton.', alg_opts{i})
    end

    % Initialize param so that any corrections are included in the report
    param{end} = param_init(param{end});
    param{end} = mp_param_init(param{end});
end

% Remove entries when skel, musc, and param are simultaneously empty
ind = cellfun(@isempty, skel) & cellfun(@isempty, musc) & cellfun(@isempty, param);
skel = skel(~ind);
musc = musc(~ind);
param = param(~ind);

%% Verify configurations for debugging
if nargin == 0
    for i = 1:length(skel)
        if isstruct(musc{i})
            fprintf('(%d) skel: %s, musc: %s\n', i, skel{i}, musc{i}.id);
        else
            fprintf('(%d) skel: %s, musc: %s\n', i, skel{i}, musc{i});
        end
        if ~isempty(param{i})
            disp(param{i})
        else
            fprintf('\n')
        end
    end
end
end
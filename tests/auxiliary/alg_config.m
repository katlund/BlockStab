function [skel, musc, param] = alg_config(config_file)
% [skel, musc, param] = ALG_CONFIG(config_file) configures cells skel,
% musc, and param for RUNKAPPAPLOT.

%%
% Default for debugging
if nargin == 0
    config_file = 'demo.json';
    disp(config_file)
end

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
            for j = 1:n_param_opts
                musc{end+1} = alg_opts{i};
                param{end+1} = param_opts(j);
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
                n_param_opts = 0;
            end

            if n_param_opts > 0
                % Cycle through all possible combinations
                for j = 1:n_param_opts
                    % Prescribe muscles
                    switch alg_opts{i}
                        case 'bcgs_sror'
                            % Prescribe skeleton
                            skel{end+1} = alg_opts{i};
        
                            % Prescribe parameters
                            param{end+1} = param_opts(j);
                            
                            % Muscle is fixed as cgs_sror
                            musc{end+1} = 'cgs_sror';
                            
                        case {'bcgs_iro_ls', 'bcgs_iro_ls_mp'}
                            % Prescribe skeleton
                            skel{end+1} = alg_opts{i};
        
                            % Prescribe parameters
                            param{end+1} = param_opts(j);
    
                            % All muscles ignored
                            musc{end+1} = [];
                            
                        otherwise
                            n_musc_opts = length(musc_opts);
                            if n_musc_opts > 0
                                for k = 1:n_musc_opts
                                    if isempty(alg_struct.(alg_opts{i}).(musc_opts{k}))
                                        % No muscle-level parameters
                                        skel{end+1} = alg_opts{i};
                                        musc{end+1} = musc_opts{k};
                    
                                        % Prescribe skeleton-level parameters
                                        param{end+1} = param_opts(j);
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
                                            param{end+1} = catstruct(...
                                                param_opts(j),...
                                                musc_param_opts(ll));
                                        end
                                    end
                                end
                            else
                                % Prescribe skeleton and default muscle
                                skel{end+1} = alg_opts{i};
                                musc{end+1} = [];
            
                                % Prescribe skeleton-level parameters
                                param{end+1} = param_opts(j);
                            end
                    end
                end
            else
                % Assign default skeleton parameters
                    switch alg_opts{i}
                        case 'bcgs_sror'
                            % Prescribe skeleton
                            skel{end+1} = alg_opts{i};
        
                            % Prescribe parameters
                            param{end+1} = [];
                            
                            % Muscle is fixed as cgs_sror
                            musc{end+1} = 'cgs_sror';
                            
                        case {'bcgs_iro_ls', 'bcgs_iro_ls_mp'}
                            % Prescribe skeleton
                            skel{end+1} = alg_opts{i};
        
                            % Prescribe parameters
                            param{end+1} = [];
    
                            % All muscles ignored
                            musc{end+1} = [];
                            
                        otherwise
                            n_musc_opts = length(musc_opts);
                            if n_musc_opts > 0
                                for k = 1:n_musc_opts
                                    if isempty(alg_struct.(alg_opts{i}).(musc_opts{k}))
                                        % No muscle-level parameters
                                        skel{end+1} = alg_opts{i};
                                        musc{end+1} = musc_opts{k};
                    
                                        % Prescribe skeleton-level parameters
                                        param{end+1} = [];
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
                                            param{end+1} = musc_param_opts(ll);
                                        end
                                    end
                                end
                            else
                                % Prescribe skeleton and default muscle
                                skel{end+1} = alg_opts{i};
                                musc{end+1} = [];
            
                                % Prescribe skeleton-level parameters
                                param{end+1} = [];
                            end
                    end
            end
        end
    else
        warning('%s is not a known muscle or skeleton.', alg_opts{i})
    end
end

%% Verify configurations for debugging
if nargin == 0
    for i = 1:length(skel)
        fprintf('(%d) skel: %s, musc: %s\n', i, skel{i}, musc{i});
        if ~isempty(param{i})
            disp(param{i})
        else
            fprintf('\n')
        end
    end
end
end
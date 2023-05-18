function [skel, musc, param] = alg_config(config_file)
% [skel, musc, param] = ALG_CONFIG(json_file)

%%
% Default for debugging
if nargin == 0
    config_file = 'sample_config.json';
end

% Convert .json to struct
alg_struct = jsondecode(fileread(config_file));

% Extract top-level fieldnames
alg_opts = fieldnames(alg_struct);
n_alg_opts = length(alg_opts);

% Extract all muscles and skeletons
musc_list = {dir('..\main\muscles\').name};
musc_list(1:2) = [];
skel_list = {dir('..\main\skeletons\').name};
skel_list(1:2) = [];

% Construct algorithm configurations
skel = [];
musc = [];
param = [];
for i = 1:n_alg_opts
    % Set up muscle as a main solver
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

    % Set up skeleton-muscle as a main solver
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
                n_param_opts = 1;
            end
            n_musc_opts = length(musc_opts);

            % Cycle through all possible combinations
            for j = 1:n_param_opts
                for k = 1:n_musc_opts
                    % Prescribe muscles
                    switch alg_opts{i}
                        case 'bcgs_sror'
                            % Prescribe skeleton
                            skel{end+1} = alg_opts{i};
        
                            % Prescribe parameters
                            if n_param_opts == 1
                                param{end+1} = [];
                            else
                                param{end+1} = param_opts(j);
                            end

                            % Muscle is fixed as cgs_sror
                            musc{end+1} = 'cgs_sror';
                            break
                            
                        case {'bcgs_iro_ls', 'bcgs_iro_ls_mp'}
                            % Prescribe skeleton
                            skel{end+1} = alg_opts{i};
        
                            % Prescribe parameters
                            if n_param_opts == 1
                                param{end+1} = [];
                            else
                                param{end+1} = param_opts(j);
                            end

                            % All muscles ignored
                            musc{end+1} = [];
                            break
                            
                        otherwise
                            if isempty(alg_struct.(alg_opts{i}).(musc_opts{k}))
                                % No muscle-level parameters
                                skel{end+1} = alg_opts{i};
                                musc{end+1} = musc_opts{k};
            
                                % Prescribe skeleton-level parameters
                                if n_param_opts == 1
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
                                    if n_param_opts == 1
                                        param{end+1} =...
                                            musc_param_opts(ll);
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
        end
    else
        warning('%s is not a known muscle or skeleton.', alg_opts{i})
    end
end

% Verify configurations for debugging
if nargin == 0
    for i = 1:length(skel)
        fprintf('skel: %s, musc: %s\n', skel{i}, musc{i});
        disp(param{i})
    end
end
end
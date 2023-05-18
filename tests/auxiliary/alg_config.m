function alg_list = alg_config(config_file)
% alg_list = ALG_CONFIG(json_file)

%%
% Param list; any other provided parameters are ignored
param_list = {'chol', 'global_scale', 'mp_package', ...
    'mp_digits', 'rpltol', 'verbose'};

% Convert .json to struct
alg_struct = jsondecode(fileread(config_file));

% Extract top-level fieldnames
fn1 = fieldnames(alg_struct);
n_fn1 = length(fn1);

% Extract all muscles and skeletons
musc_list = {dir('..\main\muscles\').name};
musc_list(1:2) = [];
skel_list = {dir('..\main\skeletons\').name};
skel_list(1:2) = [];

% Construct algorithm configurations
skel_config = [];
musc_config = [];
param_config = [];
for i = 1:n_fn1
    % Set up muscle
    if any(strcmp(musc_list, [fn1{i} '.m']))


    % Set up skeleton-muscle
    elseif any(strcmp(skel_list, [fn1{i} '.m']))
        % Extract parameters and muscles for this specific skeleton
        fn2 = fieldnames(alg_struct.(fn1{i}));

        switch fn1{i}
            case 'bcgs_sror'
                % Only cgs_sror, cgs_ror allowed
                % No chol options, no mp options, no global; only rpltol is
                % parameter


            case {'bcgs_iro_ls', 'bcgs_iro_ls_mp'}
                % Muscle is fixed
                

            otherwise
                % Split into muscles and parameters, which apply for all
                % muscles at this level


        end

    else
        warning('%s is not a known muscle or skeleton.', fn1{i})
    end

end

% Return algorithm list
alg_list = struct( ...
    'skel', skel_config, ...
    'musc', musc_config, ...
    'param', param_config);

end
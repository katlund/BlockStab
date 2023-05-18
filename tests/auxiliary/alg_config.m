function alg_list = alg_config(config_file)
% alg_list = ALG_CONFIG(json_file)

%%
% Param list; any other provided parameters are ignored
param_list = {'chol', 'mp_package', 'mp_digits', 'rpltol'};

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
    % Set up muscle as main solver
    if any(strcmp(musc_list, [fn1{i} '.m']))
        if isempty(alg_struct.(fn1{i}))
            % Then set defaults for this muscle
            skel_config{end+1} = [];
            musc_config{end+1} = fn1{i};
            param_config{end+1} = [];

        else
            % Extract parameter arrays
            fn2 = fieldnames(alg_struct.(fn1{i}));
            for k = 1:length(fn2)
                skel_config{end+1} = [];
                musc_config{end+1} = fn1{i};
                param_config{end+1} = alg_struct.(fn1{i}).(fn2{k});
            end

    % Set up skeleton-muscle as main solver
    elseif any(strcmp(skel_list, [fn1{i} '.m']))

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
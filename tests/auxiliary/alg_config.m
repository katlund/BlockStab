function alg_list = alg_config(config_file)
% alg_list = ALG_CONFIG(json_file)

%%
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

        if isempty(alg_struct.(fn1{i}))
            % Then set muscle and param defaults
            skel_config{end+1} = fn1{i};
            musc_config{end+1} = [];
            param_config{end+1} = [];

        else
            % Extract muscles for this specific skeleton
            fn2 = fieldnames( alg_struct.(fn1{i}) );
            n_fn2 = length(fn2);

            for j = 1:n_fn2
                if isempty( alg_struct.(fn1{i}).(fn2{j}) )
                    % Then set param defaults
                    skel_config{end+1} = fn1{i};
                    musc_config{end+1} = fn2{j};
                    param_config{end+1} = [];

                else
                    % Extract param options for this specific muscle
                    fn3 = fieldnames( alg_struct.(fn1{i}).(fn2{j}) );
                    n_fn3 = length(fn3);
                    for k = 1:n_fn3
                        skel_config{end+1} = fn1{i};
                        musc_config{end+1} = fn2{j};
                        param = [];
                        param_config{end+1} = param;
                    end
                end

            end
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
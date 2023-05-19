function run_data = RunKappaPlot(options, config_file)
% run_data = RUNKAPPAPLOT(options, config_file) is a test manager for all
% KappaPlot tests.
%
% INPUTS:
% - options struct with the following fields:
%   - .mat_type: 'default', 'glued', 'laeuchli', and 'monomial';
%      trial matrices are stored in the cell XX
%      default: 'default'
%   - .scale: vector determining how condition numbers of specified
%      mat_type vary
%      default: built-in for each mat_type
%   - .num_rows: (m) scalar denoting the number of rows in each trial
%      matrix
%      default: 100
%   - .num_partitions: (p) scalar denoting the number of block partitions
%      in each trial matrix
%      default: 10
%   - .block_size: (s) scalar denoting size of block partitions in each
%      trial matrix
%      default: 2
%   - .save_eps: boolean for whether to save figures as .eps files
%      default: false
%   - .save_fig: boolean for whether to save figures as .fig files
%      default: false
%   - .tex_report: boolean for whether to generate a TeX report
%      default: false
%
% - config_file: string specifying a .json file encoding algorithm
%   configurations; tips to keep in mind:
%   - all fields must be lowercase
%   - the top-level field specifies either a muscle or a skeleton
%       - when a muscle alone is specified, its subfields must correspond
%       to the possible parameter fields outlined in INTRAORTHO
%           - the problem will be solved column-wise regardless of
%           block-partitioning structure
%          - each field can be made multivalued by assigning a cell array
%          of desired parameters; a Cartesian product of all parmeters will
%          be run for the muscle
%       - when a skeleton is specified, its subfields must correspond
%       to the possible parameter fields outlined in BGS; additionally
%       subfields corresponding to desired muscles must be included (and
%       their subfields prescribed as above)
%           - the problem will be solved block-wise according to
%           options.block_size
%           - each parameter subfield can be made multivalued by assigning
%           a cell array; a Cartesian product of all parmeters will be run
%   default: see demo_config.json
%
% Note on problem dimensions: trial matrices are constructed with dimension
% m x ps.  A non-block problem can be constructed from either p = 1 or s =
% 1.
%
% OUTPUT: run_data struct with the following fields:
% - condXX: vector of condition numbers of each XX{i}
% - datetime: date and time when RUNKAPPAPLOT was called
% - dir_str: directy where plots, run_data, and TeX report are saved
% - lgd: cell of legend entries corresponding to each algorithm
%   configuration
% - loss_ortho: matrix of loss of orthogonality for each XX{i} and
%   algorithm configuration pair
% - musc: cell of muscle strings generated by ALG_CONFIG; empty entries
%   indicate default
% - normXX
% - options: vector of norm(XX{i})
% - param: cell of param structs generated by ALG_CONFIG; empty entries
%   indicate default
% - rel_chol_res: matrix of relative Cholesky residuals for each XX{i} and
%   algorithm configuration pair
% - rel_res: matrix of relative residuals for each XX{i} and
%   algorithm configuration pair
% - skel: cell of skeleton strings generated by ALG_CONFIG; empty entries
%   indicate a pure muscle algorithm
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
% Set time and date
dtnow = datetime('now');

% Defaults
if nargin == 0
    options = options_init;
    config_file = 'demo_config.json';
elseif nargin == 1
    options = options_init(options);
    config_file = 'demo_config.json';
elseif nargin == 2
    options = options_init(options);
end

% Extract dimensions
m = options.num_rows;
p = options.num_partitions;
s = options.block_size;
n = p*s;
I = eye(n);

%% Build algorithm configurations
[skel, musc, param] = alg_config(config_file);

%% Set up matrices
n_mat = length(options.scale);
XX = cell(n_mat,1);
switch lower(options.mat_type)

    case 'default'
        rng(1); U = orth(randn(m,n));
        rng(2); V = orth(randn(n,n));
        for i = 1:n_mat
            Sigma = diag(logspace(0, options.scale(i), n)');
            XX{i} = U * Sigma * V';
        end
        
    case 'glued'
        % Fix glued dimensions
        factors = factor(n);
        mid_ind = round(length(n)/2);
        r = prod(factors(1:mid_ind));
        t = prod(factors(mid_ind+1:end));
        for i = 1:n_mat
            % Create glued matrix; if s matches glued block size, results
            % are skewed
            XX{i} = create_glued_matrix(m, r, t,...
                .5 * options.scale(i), options.scale(i));
        end
        
    case 'laeuchli'
        for i = 1:n_mat
            XX{i} = laeuchli(m, n, options.scale(i));
        end
        
    case 'monomial'
        A = spdiags(linspace(.1,1,m)',0,m,m);
        for i = 1:n_mat
            mat_s = options.scale(i); mat_p = n/mat_s;
            rng(4); Y = rand(m,mat_p*mat_s);
            pp = 1:mat_p;
            Y(:,pp) = Y(:,pp) / norm(Y(:,pp));
            for k = 2:mat_s
                pp = pp + mat_p;
                Y(:,pp) = A * Y(:,pp - mat_p);
            end
            
            % Reshape
            Z = zeros(m,n);
            ind = 1:mat_s:n;
            kk = 1:mat_p;
            for k = 1:mat_s
                Z(:,ind) = Y(:,kk);
                ind = ind + 1;
                kk = kk + mat_p;
            end
            XX{i} = Z;
        end
end

%% Loop through alg_list and compute loo, rel_res, rel_chol_res
n_alg = length(skel);
loss_ortho = zeros(n_mat, n_alg);
rel_res = zeros(n_mat, n_alg);
rel_chol_res = zeros(n_mat, n_alg);
condXX = zeros(n_mat, 1);
normXX = zeros(n_mat, 1);

for i = 1:n_mat
    % Display matrix ID
    fprintf('Matrix #%d\n', i);

    % Compute cond(XX)
    condXX(i) = cond(XX{i});
    
    % Compute norm(XX)
    normXX(i) = norm(XX{i}, 2);

    for j = 1:n_alg
        if isempty(skel{j})
            % Call IntraOrtho muscle
            [QQ, RR] = IntraOrtho(XX{i}, musc{j}, param{j});

            % Display algorithm configuration
            fprintf('{%d} musc: %s\n', j, musc{j});
            if ~isempty(param{j})
                disp(param{j})
            else
                fprintf('\n')
            end
            
        else
            % Call BGS skeleton-muscle configuration
            [QQ, RR] = BGS(XX{i}, s, skel{j}, musc{j}, param{j});

            % Display algorithm configuration
            fprintf('{%d} skel: %s, ', j, skel{j})
            if ~isempty(musc{j})
                fprintf('musc: default\n');
            else
                fprintf('musc: %s\n', musc{j});
            end
                
            if ~isempty(param{j})
                disp(param{j})
            else
                fprintf('\n')
            end
        end
    
        % Compute loss of orthonormality
        loss_ortho(i,j) = norm(I - QQ' * QQ, 2);
    
        % Compute relative residual
        rel_res(i,j) = norm(XX{i} - QQ * RR, 2) / normXX(j);
        
        % Compute relative residual for Cholesky relation
        rel_chol_res(i,j) = norm(XX{i}' * XX{i} - RR' * RR, 2) / normXX(j)^2;
        
        % Clear computed variables before next run
        clear QQ RR
    end
    fprintf('\n')
end

%% Save data
% Create directory
dir_str = sprintf('results/%s_m%d_p%d_s%d', options.mat_type, m, p, s);
mkdir(dir_str)

% Build legend
lgd = alg_string(skel, musc);

% Build struct
run_data = struct( ...
    'condXX', {condXX}, ...
    'datetime', dtnow, ...
    'dir_str', dir_str, ...
    'lgd', {lgd}, ...
    'loss_ortho', {loss_ortho}, ...
    'musc', {musc}, ...
    'normXX', {normXX}, ...
    'options', options, ...
    'param', {param}, ...
    'rel_chol_res', {rel_chol_res}, ...
    'rel_res', {rel_res}, ...
    'skel', {skel});

% Save run_data
save_str = sprintf('%s/run_data', dir_str);
save(save_str, 'run_data');

%% Generate plots
alg_cmap = lines(n_alg);
symb = {'s-', 'o-', '*-', '^-', 'p-', '.-', 'h-', 'd-'};
n_symb = length(symb);

% Initialize figures and axes
fg = cell(1,3); ax = cell(1,3);
for k = 1:3
    fg{k} = figure;
    ax{k} = gca;
    hold on;
end

% Plot data
for j = 1:n_alg
    k = mod(j,n_symb) + 1;
    plot(ax{1}, condXX, loss_ortho(:,j),...
        symb{k}, 'Color', alg_cmap(j,:), 'MarkerSize', 10, 'LineWidth', 1);
    plot(ax{2}, condXX, rel_res(:,j),... 
        symb{k}, 'Color', alg_cmap(j,:), 'MarkerSize', 10, 'LineWidth', 1);
    plot(ax{3}, condXX, rel_chol_res(:,j),...
        symb{k}, 'Color', alg_cmap(j,:), 'MarkerSize', 10, 'LineWidth', 1);
end
% Plot comparison lines
plot(ax{1}, condXX, eps*condXX, 'k--', condXX, eps*(condXX.^2), 'k-')

% Make plots pretty and save them
plot_str = {'loss_ortho', 'rel_res', 'rel_chol_res'};
for k = 1:3
    % Aesthetics
    set(ax{k}, 'Yscale', 'log', 'Xscale', 'log',...
        'XGrid', 'on', 'YGrid', 'on',...
        'XMinorGrid', 'off', 'YMinorGrid', 'off',...
        'FontSize', 14);
    
    % X-axis label
    xlabel(ax{k}, '$\kappa(\mathcal{X})$', ...
        'Interpreter', 'Latex', ...
        'FontSize', 16)
    
    % Legends and titles
    if k == 1
        lgd_loo = lgd;
        lgd_loo{end+1} = '$O(\varepsilon) \kappa(\mathcal{X})$'; %#ok<*AGROW> 
        lgd_loo{end+1} = '$O(\varepsilon) \kappa^2(\mathcal{X})$';
        legend(ax{k}, lgd_loo, 'Location', 'NorthWest', ...
            'Interpreter', 'Latex', ...
            'FontSize', 10, ...
            'EdgeColor','none', ...
            'Color','none');
        title(ax{k}, ...
            ['Loss of Orthogonality: ' ...
            '$\Vert I - \bar\mathcal{Q}^T \bar\mathcal{Q}\Vert$'], ...
            'Interpreter', 'Latex', ...
            'FontSize', 16);
        movegui(fg{1},'northwest')
    else
        legend(ax{k}, lgd, 'Location', 'NorthWest', ...
            'Interpreter', 'Latex', ...
            'FontSize', 10, ...
            'EdgeColor', 'none', ...
            'Color', 'none');
        if k == 2
            title(ax{k}, ...
                ['Relative Residual: ' ...
                '$\Vert \mathcal{X} - ' ...
                '\bar\mathcal{Q}\bar\mathcal{R}\Vert/\Vert X\Vert$'], ...
                'Interpreter', 'Latex', ...
                'FontSize', 16);
            movegui(fg{2},'northeast')
        elseif k == 3
            title(ax{k}, ...
                ['Relative Cholesky Residual: ' ...
                '$\Vert \mathcal{X}^T \mathcal{X} - ' ...
                '\bar\mathcal{R}^T \bar\mathcal{R}\Vert/ ' ...
                '\Vert \mathcal{X}\Vert^2$'], ...
                'Interpreter', 'Latex', ...
                'FontSize', 16);
            movegui(fg{3},'south')
        end
    end

    % Save figures
    save_str = sprintf('%s/%s', dir_str, plot_str{k});
    if options.save_eps
        saveas(fg{k}, save_str, 'epsc');
    end

    if options.save_fig
        savefig(fg{k}, save_str, 'compact');
    end
end

%% Generate TeX report
if options.tex_report && options.save_eps
    tex_report(run_data);
end
end
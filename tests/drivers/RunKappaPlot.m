function run_data = RunKappaPlot(mat_type, options, config_file, var_eps)
% run_data = RUNKAPPAPLOT(mat_type, options, config_file, var_eps) is a
% test manager for all KappaPlot tests.
%
% INPUTS:
% - mat_type: 'default', 'glued', 'laeuchli', 'monomial', 'piled';
%   trial matrices are stored in the cell XX throughout this routine
%   default: 'default'
%
% - options struct with the following fields:
%   - .scale: vector determining how condition numbers of specified
%      mat_type vary
%      default: built-in for each mat_type; see OPTIONS_INIT
%   - .num_rows: (m) scalar denoting the number of rows in each trial
%      matrix
%      default: 100; 120 for 'monomial'
%   - .num_partitions: (p) scalar denoting the number of block partitions
%      in each trial matrix
%      default: 10; 20 for 'monomial'
%   - .block_size: (s) scalar denoting size of block partitions in each
%      trial matrix
%      default: 2; 6 for 'monomial'
%   - .save_eps: boolean for whether to save figures as .eps files
%      default: false
%   - .save_fig: boolean for whether to save figures as .fig files
%      default: false
%   - .tex_report: boolean for whether to generate a TeX report
%      default: true
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
%   default: see demo.json
%
% - var_eps: float specifying the desired machine precision; determines
%   how comparison lines are plotted for O(eps)kappa(X);
%   default: eps (2^(-52))
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
% Part of [BlockStab](https://github.com/katlund) package.  Check README
% for how to properly cite and reuse this file.

%%
% Set time and date
dtnow = datetime('now');

% Turn off specific warnings
warning('off','MATLAB:nearlySingularMatrix')

% Defaults
if nargin == 0
    mat_type = 'default';
    options = options_init(mat_type);
    config_file = 'demo.json';
    var_eps = eps;
elseif nargin == 1
    if isempty(mat_type)
        mat_type = 'default';
    end
    options = options_init(mat_type);
    config_file = 'demo.json';
    var_eps = eps;
elseif nargin == 2
    if isempty(mat_type)
        mat_type = 'default';
    end
    if isempty(options)
        options = options_init(mat_type);
    else
        options = options_init(mat_type, options);
    end
    config_file = 'demo.json';
    var_eps = eps;
elseif nargin == 3
    if isempty(mat_type)
        mat_type = 'default';
    end
    if isempty(options)
        options = options_init(mat_type);
    else
        options = options_init(mat_type, options);
    end
    if isempty(config_file)
        config_file = 'demo.json';
    end
    var_eps = eps;
elseif nargin == 4
    if isempty(mat_type)
        mat_type = 'default';
    end
    if isempty(options)
        options = options_init(mat_type);
    else
        options = options_init(mat_type, options);
    end
    if isempty(config_file)
        config_file = 'demo.json';
    end
end

% Extract dimensions
m = options.num_rows;
p = options.num_partitions;
s = options.block_size;
n = p*s;
I = eye(n);

%% Build algorithm configurations
[skel, musc, param] = alg_config(config_file);

% Set up storage for musc IDs
musc_id = cell(size(musc));

%% Set up matrices
n_mat = length(options.scale);
XX = cell(n_mat,1);
switch lower(mat_type)

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
            rng(4); Y = rand(m, mat_p*mat_s);
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

    case 'piled'
        for i = 1:n_mat
            XX{i} = create_piled_matrix(m, p, s, options.scale(i));
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
        % Extra musc_id from multiIOs
        if isstruct(musc{j})
            musc_id{j} = musc{j}.id;
        else
            musc_id{j} = musc{j};
        end
        if isempty(skel{j})
            % Call IntraOrtho muscle
            [QQ, RR] = IntraOrtho(XX{i}, musc{j}, param{j});

            % Display algorithm configuration
            fprintf('(%d) musc: %s\n', j, musc_id{j});
            if ~isempty(param{j})
                disp(param{j})
            else
                fprintf('\n')
            end
            
        else
            % Call BGS skeleton-muscle configuration
            [QQ, RR] = BGS(XX{i}, s, skel{j}, musc{j}, param{j});

            % Display algorithm configuration
            fprintf('(%d) skel: %s, ', j, skel{j})
            if isempty(musc{j})
                fprintf('musc: default\n');
            else
                fprintf('musc: %s\n', musc_id{j});
            end
                
            if ~isempty(param{j})
                disp(param{j})
            else
                fprintf('\n')
            end
        end
    
        % Compute loss of orthonormality
        loss_ortho(i,j) = norm(I - InnerProd(QQ, QQ, musc{j}), 2);
    
        % Compute relative residual
        rel_res(i,j) = norm(XX{i} - QQ * RR, 2) / normXX(i);
        
        % Compute relative residual for Cholesky relation
        rel_chol_res(i,j) = norm(InnerProd(XX{i}, XX{i}, musc{j})...
            - RR' * RR, 2) / normXX(i)^2;
        
        % Clear computed variables before next run
        clear QQ RR
    end
    fprintf('\n')
end

%% Save run data and plot color and symbol order
% Create directory
dir_str = sprintf('results/%s/%s_m%d_p%d_s%d', ...
    config_file(1:end-5), mat_type, m, p, s);
mkdir(dir_str)

% Build legend
lgd = alg_string(skel, musc);

% Line colors play a big role in whether the plots are legible
if n_alg <= 7
    % Start with built-in lines
    alg_cmap = lines(n_alg);
else
    % Use LINSPECER for beautiful, distinguishable colors
    alg_cmap = [lines(7); linspecer(n_alg-7, 'qualitative')];
end

% Randomly permute symbols to avoid matching colors for repeated symbols;
% save fixed symbol order for gen_plots (necessary for selective plotting
% as in ROADMAP)
symb = {'s-', 'o-', '*-', '^-', 'p-', 'h-', 'd-', ...
    's:', 'o:', '*:', '^:', 'p:', 'h:', 'd:'};
n_symb = length(symb);
rng(4);
if n_alg <= n_symb
    symb = symb(randperm(n_alg));
else
    n_loops = floor(n_alg / n_symb);
    symb = symb(randperm(n_symb));
    for i = 1:n_loops-1
        symb = [symb symb(randperm(n_symb))]; %#ok<AGROW> 
    end
    symb = [symb symb(randperm(rem(n_alg, n_symb)))];
end

% Build struct
run_data = struct( ...
    'condXX', {condXX}, ...
    'config_file', config_file, ...
    'dtnow', dtnow, ...
    'dir_str', dir_str, ...
    'lgd', {lgd}, ...
    'loss_ortho', {loss_ortho}, ...
    'mat_type', mat_type, ...
    'musc', {musc}, ...
    'musc_id', {musc_id}, ...
    'normXX', {normXX}, ...
    'options', options, ...
    'param', {param}, ...
    'rel_chol_res', {rel_chol_res}, ...
    'rel_res', {rel_res}, ...
    'skel', {skel}, ...
    'alg_cmap', {alg_cmap}, ...
    'symb', {symb});

% Save run data
save_str = sprintf('%s/run_data', dir_str);
save(save_str, 'run_data');
fprintf('MAT file saved in %s\n', dir_str);

%% Generate plots
gen_plots(run_data, [], var_eps);

%% Generate TeX report
if options.tex_report && (options.save_pdf || options.save_eps)
    tex_report(run_data);
    fprintf('TEX file saved in %s\n', dir_str);
end
end
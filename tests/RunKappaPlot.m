function alg_list = RunKappaPlot(options, alg_config)
% alg_lsit = RUNKAPPAPLOT(options, alg_config) is a test manager for all
% KappaPlot tests.
%
% INPUTS:
% - options struct with the following fields:
%   - .matrix_type: 'default', 'glued', 'laeuchli', and 'monomial'
%      default: 'default'
%   - .scale: vector determining how condition numbers of specified
%      matrix_type vary
%      default: built-in for each matrix_type
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
% - alg_config: string specifying a .json file encoding algorithm
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
%   default: see sample_config.json
%
% Note on problem dimensions: trial matrices are constructed with dimension
% m x ps.  A non-block problem can be constructed from either p = 1 or s =
% 1.
%
% OUTPUT: alg_list struct corresponding to alg_config
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
% Defaults
options = options_init(options);

% Build algorithm configurations
alg_list = jsondecode(fileread(alg_config));

% Kappa plots
switch lower(options.matrix_type)

    case 'default'

        
    case 'glued'
        
        
    case 'laeuchli'
        
        
    case 'monomial'
        
        
end

%% Save data
folder_str = sprintf('results/%s_m%d_p%d_s%d', fstr, m, p, s);
mkdir(folder_str)
savestr = sprintf('%s/out', folder_str);
save(savestr, 'x', 'loss_ortho','res');

%% Generate plots
skel_cmap = lines(nskel);
musc_lbl = {'s-', 'o-', '*-', '^-', 'p-', '.-', 'h-', 'd-'};

x = XXcond;
lgd_str = {};

% Initialize figures and axes
fg = cell(1,3); ax = cell(1,3);
for i = 1:3
    fg{i} = figure;
    ax{i} = gca;
    hold on;
end

% Plot data
for j = 1:nskel
    for k = 1:nmusc
        plot(ax{1}, x, loss_ortho(:,j,k),...
            musc_lbl{k}, 'Color', skel_cmap(j,:),'MarkerSize',10,'LineWidth',1);
        plot(ax{2}, x, res(:,j,k),... 
            musc_lbl{k}, 'Color', skel_cmap(j,:),'MarkerSize',10,'LineWidth',1);
        plot(ax{3}, x, res_chol(:,j,k),...
            musc_lbl{k}, 'Color', skel_cmap(j,:),'MarkerSize',10,'LineWidth',1);
        lgd_str{end+1} = sprintf('%s $\\circ$ %s', skel_str{j}, musc_str{k});
    end
end
plot(ax{1}, x, eps*x, 'k--', x, eps*(x.^2), 'k-')

% Make plots pretty and save figures
pretty_kappa_plot(fg, ax, lgd_str, folder_str);

% Save plots
if options.save_eps


% Generate report


end
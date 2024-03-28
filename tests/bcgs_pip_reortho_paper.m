% Script to reproduce plots in paper.  Note that more tests are generated
% than are presented in the paper.  The subroutine GEN_PLOTS extracts a
% subset of results to display.
%
% Algorithm configurations should be indexed as follows:
% (1) BCGS-PIP(CholQR)
% (2) BCGS-PIP(HouseQR)
% (3) BCGS-PIP+(CholQR)
% (4) BCGS-PIP+(HouseQR)
% (5) BCGS-PIPI+(CholQR)
% (6) BCGS-PIPI+(HouseQR)
% (7) BCGS-PIP^{MP}(CholQR)
% (8) BCGS-PIP^{MP}(HouseQR)
% (9) BCGS-PIP+^{MP}(CholQR)
% (10) BCGS-PIP+^{MP}(HouseQR)
% (11) BCGS-PIPI+^{MP}(CholQR)
% (12) BCGS-PIPI+^{MP}(HouseQR)

%% SET UP
% Specify algorithm configuration file
config_file = [mfilename '.json'];

% Set up options struct to be reused by GEN_PLOTS
options = [];
options.save_eps = true;
options.save_pdf = true;
options.save_fig = false;

%% Kappa plots
% Default kappa plot
mat_type = 'default';
options_default.num_rows = 100;
options_default.num_partitions = 10;
options_default.block_size = 2;
run_data_default = RunKappaPlot(mat_type, options_default, config_file);
run_data_default.options = options;
close all;

% Glued kappa plot
mat_type = 'glued';
options_glued.num_rows = 100;
options_glued.num_partitions = 10;
options_glued.block_size = 2;
options_glued.scale = 1:12;
run_data_glued = RunKappaPlot(mat_type, options_glued, config_file);
run_data_glued.options = options;
close all;

% Monomial kappa plot
mat_type = 'monomial';
options_monomial.num_rows = 2000;
options_monomial.num_partitions = 120;
options_monomial.block_size = 8;
run_data_monomial = RunKappaPlot(mat_type, options_monomial, config_file);
run_data_monomial.options = options;
close all;

% Piled kappa plot
mat_type = 'piled';
options_piled.num_rows = 1000;
options_piled.num_partitions = 10;
options_piled.block_size = 10;
run_data_piled = RunKappaPlot(mat_type, options_piled, config_file);
run_data_piled.options = options;
close all;

%% GEN_PLOTS
% Figure 1 (glued)
new_dir_str = sprintf('pip_vs_reorth/%s', run_data.dir_str);
mkdir(new_dir_str);
ind = 1:6;
gen_plots(mod_run_data(run_data_glued, ind), new_dir_str);
close all;

% Figure 2 (monomial)
new_dir_str = sprintf('reorth/%s', run_data.dir_str);
mkdir(new_dir_str);
ind = 3:6;
gen_plots(mod_run_data(run_data_monomial, ind), new_dir_str);
close all;

% Figure 3 (default)
new_dir_str = sprintf('reorth_vs_mp/%s', run_data.dir_str);
mkdir(new_dir_str);
ind = [3:6 9:12];
gen_plots(mod_run_data(run_data_default, ind), new_dir_str);
close all;

% Figure 4 (glued)
new_dir_str = sprintf('reorth_vs_mp/%s', run_data.dir_str);
mkdir(new_dir_str);
ind = [3:6 9:12];
gen_plots(mod_run_data(run_data_glued, ind), new_dir_str);
close all;

% Figure 5 (piled)
new_dir_str = sprintf('reorth_vs_mp/%s', run_data.dir_str);
mkdir(new_dir_str);
ind = [3:6 9:12];
gen_plots(mod_run_data(run_data_piled, ind), new_dir_str);
close all;
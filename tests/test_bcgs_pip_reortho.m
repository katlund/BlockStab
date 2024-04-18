% This is a script to reproduce plots in bcgs_pip_reortho paper.  Note
% that more tests are run than are presented in the paper.  The
% subroutine GEN_PLOTS extracts a subset of results to format for the
% paper.
%
% If you have the Symbolic Math Toolbox but not Advanpix, you can still
% run these tests, but you will need to edit configs/bcgs_pip_reortho.json
% and replace "advanpix" everywhere by "vpa".  If you have neither, then
% you could qualitatively reproduce the multiprecision (MP) tests with
% "mp_pair": ["single", "double"] instead of "mp_pair": ["double", "quad"],
% which will force MP routines to use the built-in SINGLE function.
%
% This script may take a while to run on a personal laptop, even if
% Advanpix is used.  Using the Symbolic Math Toolbox instead will take
% even more time. The slowest tests are the 'monomial' matrices.  You can
% comment out any matrix class you are not interested in running.
%
% If you are running this file to compare with or reproduce results from
% the associated paper, note that a TeX report is auto-generated for each
% matrix class.  You can either compile the report, or look directly at the
% .pdf or .fig plots that are saved from each test.  The compiled report
% will have a timestamp and a summary of the all the algorithm
% configurations, which can be useful if you want to share results.
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
%
% Part of [BlockStab](https://github.com/katlund) package.  Check README
% for how to properly cite and reuse this file.

%% SET UP
% Specify algorithm configuration file
config_file = 'bcgs_pip_reortho.json';

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
options_monomial.block_size = 10;
run_data_monomial = RunKappaPlot(mat_type, options_monomial, config_file);
run_data_monomial.options = options;
close all;

% Piled kappa plot
mat_type = 'piled';
options_piled.num_rows = 100;
options_piled.num_partitions = 10;
options_piled.block_size = 5;
run_data_piled = RunKappaPlot(mat_type, options_piled, config_file);
run_data_piled.options = options;
close all;

%% GEN_PLOTS
% Figure 1 (glued)
new_dir_str = sprintf('%s/pip_vs_ro', run_data_glued.dir_str);
mkdir(new_dir_str);
ind = 1:6;
gen_plots(mod_run_data(run_data_glued, ind), new_dir_str);
close all;

% Figure 2 (monomial)
new_dir_str = sprintf('%s/ro_only', run_data_monomial.dir_str);
mkdir(new_dir_str);
ind = 3:6;
gen_plots(mod_run_data(run_data_monomial, ind), new_dir_str);
close all;

% Figure 3 (default)
new_dir_str = sprintf('%s/ro_vs_mp', run_data_default.dir_str);
mkdir(new_dir_str);
ind = [3:6 9:12];
gen_plots(mod_run_data(run_data_default, ind), new_dir_str);
close all;

% Figure 4 (glued)
new_dir_str = sprintf('%s/ro_vs_mp', run_data_glued.dir_str);
mkdir(new_dir_str);
ind = [3:6 9:12];
gen_plots(mod_run_data(run_data_glued, ind), new_dir_str);
close all;

% Figure 5 (piled)
new_dir_str = sprintf('%s/ro_vs_mp/%s', run_data_piled.dir_str);
mkdir(new_dir_str);
ind = [3:6 9:12];
gen_plots(mod_run_data(run_data_piled, ind), new_dir_str);
close all;

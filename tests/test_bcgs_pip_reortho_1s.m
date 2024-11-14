% This is a script to reproduce plots in stable one-sync reortho BCGS paper.
% Note that more tests are run than are presented in the paper.
%
% Algorithm configurations should be indexed as follows:
% (1) BCGS-PIPI+(HouseQR)
% (2) BCGSI+P-1s2s(HouseQR)
% (3) BCGSI+A(HouseQR)
% (4) BCGSI+P-2s(HouseQR)
% (5) BCGSI+P-1s(HouseQR)
% (6) BCGSI+A-1s(HouseQR)
%
% Part of [BlockStab](https://github.com/katlund) package.  Check README
% for how to properly cite and reuse this file.

%% SET UP
% Specify algorithm configuration file
config_file = 'bcgs_pip_iro_1s.json';

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


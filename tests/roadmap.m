%% ROADMAP
% Reproduce plot from "roadmap" discussion, focusing on methods with O(eps)
% first vector
mat_type = 'monomial';
options.num_rows = 1000;
options.num_partitions = 120;
options.block_size = 8;
config_file = 'roadmap.json';
run_data = RunKappaPlot(mat_type, options, config_file);
close all;

%% Stepwise plots for paper
options.save_eps = true; options.save_fig = false;

% BCGS(HouseQR, CholQR) up to BCGSI+(HouseQR, CholQR)
new_dir_str = sprintf('%s/1_bcgs_to_iro', run_data.dir_str);
mkdir(new_dir_str);
ind = 1:4;
gen_plots(mod_run_data(run_data, ind), options, new_dir_str);
close all;

% BCGSI+(HouseQR, CholQR) up to BCGSI+A(HouseQR, CholQR)
new_dir_str = sprintf('%s/2_iro_to_iro_a', run_data.dir_str);
mkdir(new_dir_str);
ind = 3:6;
gen_plots(mod_run_data(run_data, ind), options, new_dir_str);
close all;

% BCGSI+A(HouseQR, CholQR) up to BCGSI+A-3S(HouseQR, CholQR)
new_dir_str = sprintf('%s/3_a4s_to_a3s', run_data.dir_str);
mkdir(new_dir_str);
ind = 5:8;
gen_plots(mod_run_data(run_data, ind), options, new_dir_str);
close all;

% BCGSI+A-3S(CholQR) up to BCGSI+A-1S(CholQR)
new_dir_str = sprintf('%s/4_a3s_to_a1s', run_data.dir_str);
mkdir(new_dir_str);
ind = 8:10;
gen_plots(mod_run_data(run_data, ind), options, new_dir_str);
close all;

%% Better behaved example
mat_type = 'default';
config_file = 'roadmap.json';
run_data = RunKappaPlot(mat_type, [], config_file);
close all;
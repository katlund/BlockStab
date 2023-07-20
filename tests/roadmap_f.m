%% ROADMAP_F
% Reproduce plot from "roadmap" discussion, focusing on methods with O(eps)
% first vector
mat_type = 'monomial';
options.num_rows = 1000;
options.num_partitions = 120;
options.block_size = 8;
config_file = 'roadmap_f.json';
run_data = RunKappaPlot(mat_type, options, config_file);

%% File strings
dir_str = sprintf('results/%s/%s_m%d_p%d_s%d', ...
    mfilename, mat_type, options.num_rows, options.num_partitions, options.block_size);

%% Stepwise plots for paper
% BCGS(HouseQR, CholQR) up to BCGSI+(HouseQR, CholQR)
new_dir_str = sprintf('%s/1_bcgs_to_iro', dir_str);
new_run_data = run_data;
ind = 1:4;
new_run_data.loss_ortho = run_data.loss_ortho(ind);
new_run_data.rel_res = run_data.rel_res(ind);
new_run_data.rel_chol_res = run_data.rel_chol_res(ind);
new_run_data.lgd = run_data.lgd(ind);
new_run_data.alg_cmap = run_data.alg_cmap(ind,:);
new_run_data.symb = run_data.symb(ind);
gen_plots(run_data, new_dir_str);

% BCGSI+(HouseQR, CholQR) up to BCGSI+F(HouseQR, CholQR)
new_dir_str = sprintf('%s/2_iro_to_iro_f', dir_str);
ind = 3:6;
new_run_data.loss_ortho = run_data.loss_ortho(ind);
new_run_data.rel_res = run_data.rel_res(ind);
new_run_data.rel_chol_res = run_data.rel_chol_res(ind);
new_run_data.lgd = run_data.lgd(ind);
new_run_data.alg_cmap = run_data.alg_cmap(ind,:);
new_run_data.symb = run_data.symb(ind);
gen_plots(run_data, new_dir_str);

% BCGSI+F(HouseQR, CholQR) up to BCGSI+F-3S(HouseQR, CholQR)
new_dir_str = sprintf('%s/3_f4s_to_f3s', dir_str);
ind = 5:8;
new_run_data.loss_ortho = run_data.loss_ortho(ind);
new_run_data.rel_res = run_data.rel_res(ind);
new_run_data.rel_chol_res = run_data.rel_chol_res(ind);
new_run_data.lgd = run_data.lgd(ind);
new_run_data.alg_cmap = run_data.alg_cmap(ind,:);
new_run_data.symb = run_data.symb(ind);
gen_plots(run_data, new_dir_str);

% BCGSI+F-3S(CholQR) up to BCGSI+F-1S(CholQR)
new_dir_str = sprintf('%s/4_f3s_to_f1s', dir_str);
ind = 8:10;
new_run_data.loss_ortho = run_data.loss_ortho(ind);
new_run_data.rel_res = run_data.rel_res(ind);
new_run_data.rel_chol_res = run_data.rel_chol_res(ind);
new_run_data.lgd = run_data.lgd(ind);
new_run_data.alg_cmap = run_data.alg_cmap(ind,:);
new_run_data.symb = run_data.symb(ind);
gen_plots(run_data, new_dir_str);
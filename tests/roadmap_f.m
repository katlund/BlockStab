%% ROADMAP_F
% Reproduce plot from "roadmap" discussion, focusing on methods with O(eps)
% first vector
mat_type = 'monomial';
options.num_rows = 1000;
options.num_partitions = 120;
options.block_size = 8;
config_file = 'roadmap_f.json';
RunKappaPlot(mat_type, options, config_file);

%% Reconstruct plots stepwise for paper

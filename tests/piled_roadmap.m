%% ROADMAP
% Reproduce plot from "roadmap" discussion, focusing on methods with O(eps)
% first vector
mat_type = 'piled';
options.num_rows = 1000;
options.num_partitions = 10;
options.block_size = 10;
config_file = 'roadmap.json';
run_data = RunKappaPlot(mat_type, options, config_file);
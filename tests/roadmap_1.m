% Reproduce plot from "roadmap" discussion
mat_type = 'monomial';
options.num_rows = 1000;
options.num_partitions = 120;
options.block_size = 8;
config_file = 'roadmap_1.json';
RunKappaPlot(mat_type, options, config_file);
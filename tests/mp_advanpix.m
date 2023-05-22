% Look at BCGS-PIP and variants
mat_type = 'glued';
options.num_rows = 1000;
options.num_partitions = 120;
options.block_size = 8;
config_file = 'mp.json';
RunKappaPlot(mat_type, options, config_file);
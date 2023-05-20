% Look at BCGS-PIP and variants
mat_type = 'monomial';
options.num_rows = 1000;
options.num_partitions = 120;
options.block_size = 8;
config_file = 'pipi.json';
RunKappaPlot(mat_type, options, config_file);
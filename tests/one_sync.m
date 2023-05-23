% All one-sync methods in comparison with BMGS and HouseQR
mat_type = 'monomial';
options.num_rows = 1000;
options.num_partitions = 120;
options.block_size = 8;
config_file = 'one_sync.json';
RunKappaPlot(mat_type, options, config_file);
% All one-sync methods in comparison with BMGS and HouseQR
config_file = 'one_sync.json';

mat_type = {'glued', 'laeuchli', 'default'};
options.num_rows = 1000;
for i = 1:3
    RunKappaPlot(mat_type{i}, options, config_file);
    close all;
end

mat_type = 'monomial';
options.num_partitions = 120;
options.block_size = 8;
config_file = 'one_sync.json';
RunKappaPlot(mat_type, options, config_file);
close all;
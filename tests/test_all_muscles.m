% All muscles
config_file = 'all_muscles.json';

mat_type = {'glued', 'laeuchli', 'default'};
options.num_rows = 1000;
for i = 1:3
    RunKappaPlot(mat_type{i}, options, config_file);
    close all;
end

mat_type = 'monomial';
options.num_partitions = 120;
options.block_size = 8;
RunKappaPlot(mat_type, options, config_file);
close all;
% Test for all one-sync methods in comparison with BMGS and HouseQR
%
% Part of [BlockStab](https://github.com/katlund/BlockStab) package.  Check README
% for how to properly cite and reuse this file.

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
RunKappaPlot(mat_type, options, config_file);
close all;
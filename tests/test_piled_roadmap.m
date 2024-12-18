% Test for reproducing "roadmap" plots from bcgs_iro_ls paper, but with
% piled matrices.
%
% Part of [BlockStab](https://github.com/katlund/BlockStab) package.  Check README
% for how to properly cite and reuse this file.

mat_type = 'piled';
options.num_rows = 1000;
options.num_partitions = 10;
options.block_size = 10;
config_file = 'roadmap.json';
run_data = RunKappaPlot(mat_type, options, config_file);
close all;
clc
close all

mat_type = 'glued';
config_file = 'pip_mp_only.json';
RunKappaPlot(mat_type, [], config_file);

mat_type = 'glued';
config_file = 'pip_vs_reorth.json';
RunKappaPlot(mat_type, [], config_file);

mat_type = 'default';
config_file = 'reorth_vs_mp.json';
RunKappaPlot(mat_type, [], config_file);

mat_type = 'monomial';
options.num_rows = 1000;
options.num_partitions = 120;
options.block_size = 8;
config_file = 'pip_vs_reorth.json';
RunKappaPlot(mat_type, options, config_file);

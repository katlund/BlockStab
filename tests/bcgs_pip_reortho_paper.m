clc
close all

mat_type = 'glued';
config_file = 'pip_vs_reorth.json';
RunKappaPlot(mat_type, [], config_file);

mat_type = 'monomial';
options.num_rows = 2000;
options.num_partitions = 120;
options.block_size = 10;
config_file = 'reorth.json';
RunKappaPlot(mat_type, options, config_file);

mat_type = 'default';
config_file = 'reorth_vs_mp.json';
RunKappaPlot(mat_type, [], config_file);

mat_type = 'glued';
config_file = 'reorth_vs_mp.json';
RunKappaPlot(mat_type, [], config_file);
%%
mat_type = 'piled';
options.num_rows = 100;%1000;
options.num_partitions = 8;%10;
options.block_size = 2;%10;
config_file = 'pip_pipi_mp_only_2_prec.json';
run_data = RunKappaPlot(mat_type, options, config_file);






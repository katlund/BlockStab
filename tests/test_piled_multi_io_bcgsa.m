% Tests behavior of multi_io feature for bcgsa in piled matrices with
% various block sizes
mat_type = 'piled';
options.num_rows = 1000;
options.num_partitions = 100;
config_file = 'multi_io_bcgs_a_eda.json';
p = [1,100];
for i=1:length(p)
    options.block_size = p(i);
    run_data = RunKappaPlot(mat_type, options, config_file);
    close all
end

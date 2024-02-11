% Subset for presentation
mat_type = {'glued', 'monomial'};
options.num_rows = 1000;
options.num_partitions = 120;
options.block_size = 2;
for i = 1:length(mat_type)
    RunKappaPlot(mat_type{i}, options, 'pres_pip.json');
    close all;
end

mat_type = {'glued', 'monomial'};
options.num_rows = 1000;
options.num_partitions = 120;
options.block_size = 2;
for i = 1:length(mat_type)
    RunKappaPlot(mat_type{i}, options, 'pres_ls.json');
    close all;
end
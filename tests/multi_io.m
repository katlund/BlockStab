% Tests behavior of multi_io feature
mat_type = 'default';
config_file = {'multi_io_bcgs_a.json', ...
    'multi_io_bcgs_iro_a.json',...
    'multi_io_bcgs_iro_3-1s.json'};
for i = 1:3
    run_data = RunKappaPlot(mat_type, [], config_file{i});
    close all;
end
% Tests behavior of multi_io feature
mat_type = 'default';
config_file = {'multi_io_bcgs_a.json', ...
    'multi_io_bcgs_iro_a.json',...
    'multi_io_bcgs_iro_a_3s.json'};
for i = 1:3
    run_data = RunKappaPlot(mat_type, [], config_file{i});
    close all;

    % Check that the run_data matches exactly
    assert(norm(run_data.loss_ortho(:,1) - run_data.loss_ortho(:,3)) == 0); % HouseQR
    assert(norm(run_data.loss_ortho(:,2) - run_data.loss_ortho(:,4)) == 0); % CholQR
end
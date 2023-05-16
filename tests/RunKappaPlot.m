function run_data = RunKappaPlot(options)
% run_data = RUNKAPPAPLOT(options) is a test manager for all KappaPlot tests.
%
% INPUT: options struct with the following fields:
% - .matrix_type: 'default', 'glued', 'laeuchli', and 'monomial'
%    default: 'default'
% - .musc: cell array of muscles
%    default: {'HouseQR'}
% - .scale: vector determining how condition numbers of specified
%    matrix_type vary
%    default: depends on matrix_type
% - .skel: cell array of skeletons; leads to Cartesian product with muscles
%    (excluding duplicates and infeasible combinations)
%    default: {'BCGS', 'BMGS', 'BCGS_IRO'}
% - .saveas: 'eps', 'fig', etc; if unspecified or empty then plots are not
%    saved
%    default: []
% - .tex_report: true or false to generate a TeX report
%    default: false
%
% OUTPUT: run_data struct with the following fields:
% - .alg_config: cell of algorithm configurations
% - .prob_size: vector of XX matrix dimensions
% - .block_size: scalar denoting number of columns per block in XX
%%


end
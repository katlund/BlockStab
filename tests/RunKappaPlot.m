function RunKappaPlot(options)
% RUNKAPPAPLOT(options) is a test manager for all KappaPlot tests.
% options struct
% - .matrix_type: 'default', 'glued', 'laeuchli', and 'monomial'
%    default: 'default'
% - .musc: cell array of muscles
%    default: {'HouseQR'}
% - .scale: vector determining how condition numbers of specified
%    matrix_type vary; see test scripts for defaults
% - .skel: cell array of skeletons; leads to Cartesian product with muscles
%    (excluding duplicates and infeasible combinations)
%    default: {'BCGS', 'BMGS', 'BCGS_IRO'}
% - .saveas: 'eps', 'fig', etc; if unspecified or empty then plots are not
%    saved
%    default: []
% - .tex_report: true or false to generate a TeX report
%    default: false

%%


end
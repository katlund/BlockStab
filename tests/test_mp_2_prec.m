% Test for comparing multiprecision methods with two precisions.
%
% Part of [BlockStab](https://github.com/katlund) package.  Check README
% for how to properly cite and reuse this file.

mat_type = 'glued';

config_file = 'mp_2_prec.json';
RunKappaPlot(mat_type, [], config_file);

config_file = 'mp_2_prec_vpa.json';
RunKappaPlot(mat_type, [], config_file);

config_file = 'mp_2_prec_advanpix.json';
RunKappaPlot(mat_type, [], config_file);
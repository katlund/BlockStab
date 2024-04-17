function [QQ, RR] = bcgs_pip_ro(XX, s, musc, param)
% [QQ, RR] = BCGS_PIP_RO(XX, s, musc, param) performs Block Classical
% Gram-Schmidt with Pythagorean Inner Product modification and (outer)
% ReOrthogonalization on the m x n matrix XX with p = n/s block partitions
% each of size s with intra-orthogonalization procedure determined by musc.
%
% See BGS for more details about the parameters, and INTRAORTHO for musc
% options.
%
% Part of [BlockStab](https://github.com/katlund) package.  Check README
% for how to properly cite and reuse this file.

%%
% Default: debugging off
if nargin < 4
    param.verbose = 0;
end

% Run BCGS_PIP twice
[QQ, RR] = bcgs_pip(XX, s, musc, param);
[QQ, RR2] = bcgs_pip(QQ, s, musc, param);
RR = RR2 * RR;
end
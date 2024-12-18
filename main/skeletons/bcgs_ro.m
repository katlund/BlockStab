function [QQ, RR] = bcgs_ro(XX, s, musc, param)
% [QQ, RR] = BCGS_RO(XX, s, musc, param) performs Block Classical
% Gram-Schmidt with (outer) ReOrthogonalization on the m x n matrix XX with
% p = n/s block partitions each of size s with inner orthogonalization
% procedure determined by musc.  BCGS_RO is the block generalization of
% CGS_RO.
%
% See BGS for more details about the parameters, and INTRAORTHO for musc
% options.  See BCGS for the base routine that is called twice.
%
% Part of [BlockStab](https://github.com/katlund/BlockStab) package.  Check README
% for how to properly cite and reuse this file.

%%
% Default: debugging off
if nargin < 4
    param.verbose = 0;
end

% Run BCGS twice
[QQ, RR] = bcgs(XX, s, musc, param);
[QQ, RR2] = bcgs(QQ, s, musc, param);
RR = RR2 * RR;
end
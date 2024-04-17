function [QQ, RR] = bcgs_pip_ro_mp(XX, s, musc, param)
% [QQ, RR] = BCGS_PIP_RO_MP(XX, s, musc, param) performs Block Classical
% Gram-Schmidt with Pythagorean Inner Product modification and (outer)
% ReOrthogonalization on the m x n matrix XX with p = n/s block partitions
% each of size s with intra-orthogonalization procedure determined by musc.
%
% This multiprecision version computes the inputs to Cholesky and the
% Cholesky factorization itself in simulated quadruple (or other,
% user-specified precision) precision.  See MP_SWITCH for details on the
% param struct.
%
% See BGS and MP_SWITCH for more details about the parameters, and
% INTRAORTHO for musc options.
%
% Part of [BlockStab](https://github.com/katlund) package.  Check README
% for how to properly cite and reuse this file.

%%
% Default: debugging off
if nargin < 4
    param.verbose = 0;
end

% Run BCGS_PIP_MP twice
[QQ, RR] = bcgs_pip_mp(XX, s, musc, param);
[QQ, RR2] = bcgs_pip_mp(QQ, s, musc, param);
RR = RR2 * RR;
end
function [QQ, RR] = bcgs_pio_ro(XX, s, musc, param)
% [QQ, RR] = BCGS_PIP_RO(XX, s, musc, param) Block Classical Gram-Schmidt
% with Pythagorean Intra-Orthogonalization modification and (outer)
% ReOrthogonalization on the m x n matrix XX with p = n/s block partitions
% each of size s and with intra-orthogonalization procedure determined by
% musc.
%
% See BGS for more details about the parameters, and INTRAORTHO for musc
% options.
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
% Default: debugging off
if nargin < 4
    param.verbose = 0;
end

% Run BCGS_PIO twice
[QQ, RR] = bcgs_pio(XX, s, musc, param);
[QQ, RR2] = bcgs_pio(QQ, s, musc, param);
RR = RR2 * RR;
end
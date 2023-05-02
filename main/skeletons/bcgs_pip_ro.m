function [QQ, RR] = bcgs_pip_ro(XX, s, musc, verbose)
% [QQ, RR] = BCGS_PIP_IRO(XX, s, musc, verbose) performs Block Classical
% Gram-Schmidt with Pythagorean Inner Product modification and (outer)
% ReOrthogonalization on the m x n matrix XX with p = n/s block partitions
% each of size s with intra-orthogonalization procedure determined by musc.
%
% See BGS for more details about the parameters, and INTRAORTHO for musc
% options.
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
% Default: debugging off
if nargin < 4
    verbose = 0;
end

% Run BCGS_PIP twice
[QQ, RR] = bcgs_pip(XX, s, musc, verbose);
[QQ, RR2] = bcgs_pip(QQ, s, musc, verbose);
RR = RR2 * RR;
end
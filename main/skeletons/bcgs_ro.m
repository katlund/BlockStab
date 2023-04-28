function [QQ, RR] = bcgs_ro(XX, s, musc, verbose)
% [QQ, RR] = BCGS_RO(XX, s, musc, verbose) performs Block Classical
% Gram-Schmidt with (outer) ReOrthogonalization on the m x n matrix XX with
% p = n/s block partitions each of size s with inner orthogonalization
% procedure determined by musc.  BCGS_RO is the block generalization of
% CGS_RO.
%
% See BGS for more details about the parameters, and INTRAORTHO for musc
% options.  See BCGS for the base routine that is called twice.

%%
% Default: debugging off
if nargin < 4
    verbose = 0;
end

% Run BCGS twice
[QQ, RR] = bcgs(XX, s, musc, verbose);
[QQ, RR2] = bcgs(QQ, s, musc, verbose);
RR = RR2 * RR;
end
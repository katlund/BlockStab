function [Q, R] = mgs_ro(X, verbose)
% [Q, R] = MGS_RO(X, verbose) performs Modified Gram-Schmidt with (outer)
% ReOrthogonalization on the m x s matrix X.
%
% See INTRAORTHO for more details about the parameters.
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
% Default: debugging off
if nargin < 2
    verbose = 0;
end

[Q1, R1] = mgs(X, verbose);
[Q, R] = mgs(Q1, verbose);
R = R*R1;
end
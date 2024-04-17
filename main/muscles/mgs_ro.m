function [Q, R] = mgs_ro(X, verbose)
% [Q, R] = MGS_RO(X, verbose) performs Modified Gram-Schmidt with (outer)
% ReOrthogonalization on the m x s matrix X.
%
% See INTRAORTHO for more details about the parameters.

%%
% Default: debugging off
if nargin < 2
    verbose = 0;
end

[Q1, R1] = mgs(X, verbose);
[Q, R] = mgs(Q1, verbose);
R = R*R1;
end
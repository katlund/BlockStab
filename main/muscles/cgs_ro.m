function [Q, R] = cgs_ro(X, verbose)
% [Q, R] = CGS_RO(X, verbose) performs Classical Gram-Schmidt with (outer)
% ReOrthogonalization on the m x s matrix X.  See CGS for the base version.
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

% Run CGS_P twice
[Q, R] = cgs(X, verbose);
[Q, R2] = cgs(Q, verbose);
R = R2 * R;

end
function [Q, R] = cgs_p_ro(X, verbose)
% [Q, R] = CGS_P_RO(X, verbose) performs Classical Gram-Schmidt with
% Pythagorean modification and (outer) ReOrthogonalization on the m x s
% matrix X.  See CGS_P for the base version.
%
% See INTRAORTHO for more details about the parameters.
%
% Part of [BlockStab](https://github.com/katlund/BlockStab) package.  Check README
% for how to properly cite and reuse this file.

%%
% Default: debugging off
if nargin < 2
    verbose = 0;
end

% Run CGS_P twice
[Q, R] = cgs_p(X, verbose);
[Q, R2] = cgs_p(Q, verbose);
R = R2 * R;

end
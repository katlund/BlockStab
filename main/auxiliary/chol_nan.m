function [R, flag] = chol_nan(R)
% [R, nan_flag] = CHOL_NAN(R) is a simple script for catching when MATLAB's
% chol refuses to perform and returns a NaN matrix and flag when R cannot
% be computed.
%
% Part of [BlockStab](https://github.com/katlund/BlockStab) package.  Check README
% for how to properly cite and reuse this file.

%%
[~, flag] = chol(R);
if flag == 0
    R = chol(R);
else
    R = NaN(size(R));
end
end
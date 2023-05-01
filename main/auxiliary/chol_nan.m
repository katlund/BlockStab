function [R, nan_flag] = chol_nan(R)
% A simple script for catching when MATLAB's chol refuses to perform and
% returns a NaN matrix and flag when R cannot be computed.

%%
[~, nan_flag] = chol(R);
if nan_flag == 0
    R = chol(R);
else
    R = NaN(size(R));
end
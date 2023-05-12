function XY = InnerProd(X, Y, musc)
% S = InnerProd(X, Y, musc) is a wrapper function for switching between
% different inner product paradigms, determined by the char musc.
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
% Defaults
if nargin == 2
    musc = [];
end

% Switch
musc = lower(musc);
switch musc
    case {'global'}
        s = size(X,2);
        XY = (trace(X'*Y)/s) * eye(s);  % with s-scaling
        
    case {'global-no-scale'}
        s = size(X,2);
        XY = trace(X'*Y) * eye(s);
        
    otherwise
        XY = X'*Y;
        
end
end
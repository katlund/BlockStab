function XY = InnerProd(X, Y, musc)
% S = InnerProd(X, Y, musc) is a wrapper function for switching between
% different inner product paradigms, determined by the char musc.

%%
% Defaults
if nargin <= 1
    error('Two arguments must be provided.')
elseif nargin == 2
    musc = '';
elseif nargin == 3
    if isempty(musc)
        musc = '';
    end
end

% Switch
musc = lower(musc);
switch musc
    case {'global'}
        s = size(X,2);
        XY = (trace(X'*Y)/s) * eye(s);
        
    otherwise
        XY = X'*Y;
        
end
end
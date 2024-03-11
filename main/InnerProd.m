function XY = InnerProd(X, Y, musc)
% S = InnerProd(X, Y, musc) is a wrapper function for switching between
% different inner product paradigms, determined by the char musc.
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

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

% Handle raw multiIO struct
if isstruct(musc)
    musc = unpack_multi_io(musc);
end

% Handle multiIO processed by UNPACK_MULTI_IO
if iscell(musc)
    if any(strcmp(musc, 'global'))
        musc_id = 'global';
    else
        musc_id = '';
    end
end

% Handle regular char musc
if ischar(musc)
    musc_id = lower(musc);
end

% Switch
switch musc_id
    case {'global'}
        s = size(X,2);
        XY = (trace(X'*Y)/s) * eye(s);
        
    otherwise
        XY = X'*Y;
        
end
end
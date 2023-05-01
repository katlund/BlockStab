function XY = InnerProd(X, Y, str)
% S = InnerProd(X,Y) is a wrapper function for switching between different
% inner product paradigms, determined by str.
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
% Defaults
if nargin == 2
    str = 'default';
end
if nargin == 3
    if isempty(str)
        str = 'default';
    end
end

% Switch
str = lower(str);
switch str
    case {'default', 'cgs', 'cgs_p',...
            'cgs_ro', 'cgs_iro', 'cgs_iro_ls', 'cgs_sro', 'cgs_sror',...
            'mgs', 'mgs_iro', 'mgs_ro',...
            'mgs_svl', 'mgs_lts', 'mgs_icwy', 'mgs_cwy',...
            'houseqr',...
            'cholqr', 'cholqr_ro', 'sh_cholqr_roro', 'iter_cholqr',...
            'cholqr_pinv', 'cholqr_vpa', 'cholqr_mp'}
        XY = X'*Y;
        
    case {'global'}
        s = size(X,2);
        XY = (trace(X'*Y)/s) * eye(s);  % with s-scaling
        
    case {'global-no-scale'}
        s = size(X,2);
        XY = trace(X'*Y) * eye(s);
        
    otherwise
        error('%s is not a viable inner product option', str);
end
end
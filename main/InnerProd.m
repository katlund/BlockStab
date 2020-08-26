function XY = InnerProd(X, Y, str)
% S = InnerProd(X,Y) is a wrapper function for switching between different
% inner product paradigms, determined by str.

%%
str = lower(str);
switch str
    case {'cgs', 'cgs_ro', 'cgs_iro', 'cgs_iro_ls', 'cgs_sro', 'cgs_sror',...
            'mgs', 'mgs_iro', 'mgs_ro', 'mgs_svl', 'houseqr', 'cholqr',...
            'cholqr_ro', 'sh_cholqr_roro', 'svqb', 'iter_cholqr',...
            'mgs_lts', 'mgs_icwy', 'mgs_cwy', 'cgs_ls_ro',...
            'mgs_lts_trans', 'mgs_svl_trans', 'mgs_icwy_trans', 'mgs_cwy_trans'}
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
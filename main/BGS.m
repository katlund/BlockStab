function [QQ, RR, TT, TotTime] = BGS(XX, s, skel, musc, rpltol, verbose)
% [QQ, RR, TT, TotTime] = BGS(XX, s, skel, musc, rpltol, verbose) is a
% wrapper function for calling different Block Gram-Schmidt skeleton-muscle
% configurations. XX must be a matrix, s a scalar specifying the size of
% block vectors, and skel and musc char arrays.  s must divide the number
% of columns of XX. rpltol is an optional argument for BCGS_SROR, a scalar
% which determines the replacement tolerance; it is ignored otherwise.
% verbose is a Boolean for whether to print intermediate loss of
% orthogonality (LOO) or relative residual (RelRes) per iteration.
%
% For all possible muscle options, see INTRAORTHO.

%%
skel = lower(skel);
switch skel
    case {'bcgs'}
        tic;
        [QQ, RR] = bcgs(XX, s, musc, verbose);
        TotTime = toc;
        TT = eye(size(RR));
        
    case {'bcgs_iro'}
        tic;
        [QQ, RR] = bcgs_iro(XX, s, musc, verbose);
        TotTime = toc;
        TT = eye(size(RR));
        
    case {'bcgs_sror'}
        if nargin == 5
            tic;
            [QQ, RR] = bcgs_sror(XX, s, rpltol, verbose);
            TotTime = toc;
        else
            tic;
            [QQ, RR] = bcgs_sror(XX, s, [], verbose);
            TotTime = toc;
        end
        TT = eye(size(RR));
        
    case {'bcgs_iro_ls'}
        tic;
        [QQ, RR] = bcgs_iro_ls(XX, s, musc, verbose);
        TotTime = toc;
        
    case {'bmgs'}
        tic;
        [QQ, RR] = bmgs(XX, s, musc, verbose);
        TotTime = toc;
        TT = eye(size(RR));
        
    case {'bmgs_svl'}
        tic;
        [QQ, RR, TT] = bmgs_svl(XX, s, musc, verbose);
        TotTime = toc;
        
    case {'bmgs_lts'}
        tic;
        [QQ, RR, TT] = bmgs_lts(XX, s, musc, verbose);
        TotTime = toc;
        
    case {'bmgs_cwy'}
        tic;
        [QQ, RR, TT] = bmgs_cwy(XX, s, musc, verbose);
        TotTime = toc;
        
    case {'bmgs_icwy'}
        tic;
        [QQ, RR, TT] = bmgs_icwy(XX, s, musc, verbose);
        TotTime = toc;
        
% T-variants --------------------------------------------------------------
    case {'bcgs_iro_t'}
        tic;
        [QQ, RR, TT] = bcgs_iro_t(XX, s, musc, verbose);
        TotTime = toc;
        
    case {'bmgs_t'}
        tic;
        [QQ, RR, TT] = bmgs_t(XX, s, musc, verbose);
        TotTime = toc;
        
% Reorthogonalize first step ----------------------------------------------
    case {'bcgs_iro_1'}
        tic;
        [QQ, RR] = bcgs_iro_1(XX, s, musc, verbose);
        TotTime = toc;
        TT = eye(size(RR));
        
% P-variants --------------------------------------------------------------
    case {'bcgs_pio'}
        tic;
        [QQ, RR] = bcgs_pio(XX, s, musc, verbose);
        TotTime = toc;
        TT = eye(size(RR));
        
    case {'bcgs_pip'}
        tic;
        [QQ, RR] = bcgs_pip(XX, s, musc, verbose);
        TotTime = toc;
        TT = eye(size(RR));
        
    otherwise
        error('%s is not a viable skeleton option', skel);
end
function [QQ, RR, TT, TotTime] = BGS(XX, s, skel, musc, rpltol, verbose)
% [QQ, RR, TT, TotTime] = BGS(XX, s, skel, musc, rpltol, verbose) is a
% wrapper function for calling different Block Gram-Schmidt skeleton-muscle
% configurations. XX must be a matrix, s a scalar specifying the size of
% block vectors, and skel and musc char arrays.  s must cleanly divide the
% number of columns of XX. rpltol is an optional scalar argument for
% BCGS_SROR that determines the replacement tolerance; it is ignored
% otherwise. verbose is a boolean for whether to print intermediate loss of
% orthogonality (LOO) or relative residual (RelRes) per iteration.
%
% For all possible muscle options, see INTRAORTHO.
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
% Defaults
if nargin == 4
    rpltol = [];
    verbose = 0;
elseif nargin == 5
    verbose = 0;
end

skel = lower(skel);
switch skel
% Standard BCGS and reorthogonalized variants -----------------------------
    case {'bcgs'}
        tic;
        [QQ, RR] = bcgs(XX, s, musc, verbose);
        TotTime = toc;

    case {'bcgs_ro'}
        tic;
        [QQ, RR] = bcgs(XX, s, musc, verbose);
        TotTime = toc;
        
    case {'bcgs_iro'}
        tic;
        [QQ, RR] = bcgs_iro(XX, s, musc, verbose);
        TotTime = toc;

% [Stewart 2008] variant --------------------------------------------------
    case {'bcgs_sror'}
        if nargin >= 5
            tic;
            [QQ, RR] = bcgs_sror(XX, s, rpltol, verbose);
            TotTime = toc;
        else
            tic;
            [QQ, RR] = bcgs_sror(XX, s, [], verbose);
            TotTime = toc;
        end

% [Swirydowicz et al. 2020]/[Bielich et al. 2022] variant -----------------
    case {'bcgs_iro_ls'}
        tic;
        [QQ, RR] = bcgs_iro_ls(XX, s, musc, verbose);
        TotTime = toc;

% "Roadmap" variants ------------------------------------------------------
    case {'bcgs_iro_3s'}
        tic;
        [QQ, RR] = bcgs_iro_3s(XX, s, musc, verbose);
        TotTime = toc;

    case {'bcgs_iro_2s'}
        tic;
        [QQ, RR] = bcgs_iro_2s(XX, s, musc, verbose);
        TotTime = toc;
        
    case {'bcgs_iro_1s'}
        tic;
        [QQ, RR] = bcgs_iro_1s(XX, s, musc, verbose);
        TotTime = toc;

% Standard BMGS -----------------------------------------------------------
    case {'bmgs'}
        tic;
        [QQ, RR] = bmgs(XX, s, musc, verbose);
        TotTime = toc;
        
% 3-sync BMGS variants ----------------------------------------------------
    case {'bmgs_svl'}
        tic;
        [QQ, RR, TT] = bmgs_svl(XX, s, musc, verbose);
        TotTime = toc;
        
    case {'bmgs_lts'}
        tic;
        [QQ, RR, TT] = bmgs_lts(XX, s, musc, verbose);
        TotTime = toc;
        
% 1-sync BMGS variants ----------------------------------------------------
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
        
    case {'bcgs_pio_iro_1'}
        tic;
        [QQ, RR] = bcgs_pio_iro_1(XX, s, musc, verbose);
        TotTime = toc;   
        
    case {'bcgs_pip_iro_1'}
        tic;
        [QQ, RR] = bcgs_pip_iro_1(XX, s, musc, verbose);
        TotTime = toc;
        
% P-variants --------------------------------------------------------------
    case {'bcgs_pio'}
        tic;
        [QQ, RR] = bcgs_pio(XX, s, musc, verbose);
        TotTime = toc;
        
    case {'bcgs_pip'}
        tic;
        [QQ, RR] = bcgs_pip(XX, s, musc, verbose);
        TotTime = toc;

% Reorthogonalized P-variants ---------------------------------------------
    case {'bcgs_pip_ro'}
        tic;
        [QQ1, RR1] = bcgs_pip(XX, s, musc, verbose);
        [QQ, RR2] = bcgs_pip(QQ1, s, musc, verbose);
        RR = RR2 * RR1;
        TotTime = toc;

    case {'bcgs_pip_iro'}
        tic;
        [QQ, RR] = bcgs_pip_iro(XX, s, musc, verbose);
        TotTime = toc;
		
% Multiprecision variants --------------------------------------------------
    case {'bcgs_iro_ls_mp'}
        tic;
        [QQ, RR] = bcgs_iro_ls_mp(XX, s, musc, verbose);
        TotTime = toc;    

    case {'bcgs_iro_ls_vpa'}
        tic;
        [QQ, RR] = bcgs_iro_ls_vpa(XX, s, musc, verbose);
        TotTime = toc; 
    case {'bcgs_pio_mp'}
        tic;
        [QQ, RR] = bcgs_pio_mp(XX, s, musc, verbose);
        TotTime = toc;
        
    case {'bcgs_pip_mp'}
        tic;
        [QQ, RR] = bcgs_pip_mp(XX, s, musc, verbose);
        TotTime = toc;  

    case {'bcgs_pip_ro_mp'}
        tic;
        [QQ1, RR1] = bcgs_pip_mp(XX, s, musc, verbose);
        [QQ, RR2] = bcgs_pip_mp(QQ1, s, musc, verbose);
        RR = RR2 * RR1;
        TotTime = toc;  

    case {'bcgs_pip_ro'}
        tic;
        [QQ1, RR1] = bcgs_pip(XX, s, musc, verbose);
        [QQ, RR2] = bcgs_pip(QQ1, s, musc, verbose);
        RR = RR2 * RR1;
        TotTime = toc;  
        
    case {'bcgs_pio_vpa'}
        tic;
        [QQ, RR] = bcgs_pio_vpa(XX, s, musc, verbose);
        TotTime = toc;
        
    case {'bcgs_pip_vpa'}
        tic;
        [QQ, RR] = bcgs_pip_vpa(XX, s, musc, verbose);
        TotTime = toc; 
        
    case {'bcgs_pip_iro'}
        tic;
        [QQ, RR] = bcgs_pip_iro(XX, s, musc, verbose);
        TotTime = toc;
        
    case {'bcgs_pio_iro'}
        tic;
        [QQ, RR] = bcgs_pio_iro(XX, s, musc, verbose);
        TotTime = toc;
    
	case {'bcgs_pip_iro_mp'}
        tic;
        [QQ, RR] = bcgs_pip_iro_mp(XX, s, musc, verbose);
        TotTime = toc;
        
    case {'bcgs_pio_iro_mp'}
        tic;
        [QQ, RR] = bcgs_pio_iro_mp(XX, s, musc, verbose);
        TotTime = toc;

    case {'bcgs_pip_iro_vpa'}
        tic;
        [QQ, RR] = bcgs_pip_iro_vpa(XX, s, musc, verbose);
        TotTime = toc;
        
    case {'bcgs_pio_iro_vpa'}
        tic;
        [QQ, RR] = bcgs_pio_iro_vpa(XX, s, musc, verbose);
        TotTime = toc;
        
    otherwise
        error('%s is not a viable skeleton option', skel);
end
if ~exist('TT', 'var')
    TT = eye(size(RR));
end
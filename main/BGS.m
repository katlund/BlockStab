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
addpath(genpath('skeletons/'))
addpath(genpath('muscles/'))

% Defaults
if nargin == 4
    rpltol = [];
    verbose = 0;
elseif nargin == 5
    verbose = 0;
end

skel = lower(skel);
switch skel
    case {'bcgs'}
        tic;
        [QQ, RR] = bcgs(XX, s, musc, verbose);
        TotTime = toc;

    case {'bcgs_ro'}
        tic;
        [QQ1, RR1] = bcgs(XX, s, musc, verbose);
        [QQ, RR2] = bcgs(QQ1, s, musc, verbose);
        RR = RR2 * RR1;
        TotTime = toc;
        
    case {'bcgs_iro'}
        tic;
        [QQ, RR] = bcgs_iro(XX, s, musc, verbose);
        TotTime = toc;
        
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
        
    case {'bcgs_iro_ls'}
        tic;
        [QQ, RR] = bcgs_iro_ls(XX, s, musc, verbose);
        TotTime = toc;

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

    case {'bcgs_iro_ls_play'}
        tic;
        [QQ, RR] = bcgs_iro_ls_play(XX, s, musc, verbose);
        TotTime = toc;

    case {'bcgs_iro_ls_1'}
        tic;
        [QQ, RR] = bcgs_iro_ls_1(XX, s, musc, verbose);
        TotTime = toc;

    case {'bcgs_iro_ls_2'}
        tic;
        [QQ, RR] = bcgs_iro_ls_2(XX, s, musc, verbose);
        TotTime = toc;

    case {'bcgs_iro_1s'}
        tic;
        [QQ, RR] = bcgs_iro_1s(XX, s, musc, verbose);
        TotTime = toc;

    case {'bcgs_iro_2s'}
        tic;
        [QQ, RR] = bcgs_iro_2s(XX, s, musc, verbose);
        TotTime = toc;

    case {'bcgs_iro_3s'}
        tic;
        [QQ, RR] = bcgs_iro_3s(XX, s, musc, verbose);
        TotTime = toc;

    case {'bcgs_iro_bl'}
        tic;
        [QQ, RR] = bcgs_iro_bl(XX, s, musc, verbose);
        TotTime = toc;

    case {'bcgs_iro_kl'}
        tic;
        [QQ, RR] = bcgs_iro_kl(XX, s, musc, verbose);
        TotTime = toc;

% MGS variants -------------------------------------------------------------
    case {'bmgs'}
        tic;
        [QQ, RR] = bmgs(XX, s, musc, verbose);
        TotTime = toc;
        
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
        
% P-variants --------------------------------------------------------------
    case {'bcgs_pio'}
        tic;
        [QQ, RR] = bcgs_pio(XX, s, musc, verbose);
        TotTime = toc;
        
    case {'bcgs_pip'}
        tic;
        [QQ, RR] = bcgs_pip(XX, s, musc, verbose);
        TotTime = toc;
        
    case {'bcgs_pio_free'}
        tic;
        [QQ, RR] = bcgs_pio_free(XX, s, musc, verbose);
        TotTime = toc;
        
    case {'bcgs_pip_free'}
        tic;
        [QQ, RR] = bcgs_pip_free(XX, s, musc, verbose);
        TotTime = toc;
        
    otherwise
        error('%s is not a viable skeleton option', skel);
end
if ~exist('TT', 'var')
    TT = eye(size(RR));
end
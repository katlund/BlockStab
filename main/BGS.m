function [QQ, RR, TT, run_time] = BGS(XX, s, skel, musc, param)
% [QQ, RR, TT, run_time] = BGS(XX, s, skel, musc, param) is a wrapper
% function for calling different Block Gram-Schmidt skeleton-muscle
% configurations.
% 
% INPUTS:
% - XX: n x m matrix
% - s: a scalar specifying the size of block vectors (must cleanly divide
%   integer m for p = m/s block partitions)
% - skel: char array specifying inter-orthogonalation (between blocks)
%   routine
% - musc: char array specifying intra-orthogonalization (within block)
%   routine; see INTRAORTHO for all possible options
%   default: 'houseqr'
% - param: a struct with the following optional fields:
%   - .chol: char specifying what type of Cholesky subroutine to call for
%      skeletons that hard-code Cholesky via a block Pythagorean trick
%      (e.g., BCGS_PIP, BCGS_PIO, BCGS_IRO_LS, BMGS_CWY, BMGS_ICWY, and
%      their reorthogonalized and multi-precision versions)
%      default: 'chol_nan'
%   - .mp_package: char specifying either 'advanpix' or 'symbolic toolbox'
%      as the mixed precision package for routines with *_MP
%      default: []
%   - .mp_digits: int specifiying number of precision digits, e.g., 34 for
%      quadruple precision (in Advanpix) or 32 for quadruple precision in
%      Symbolic Toolbox
%      default: []
%   - .rpltol: scalar argument for BCGS_SROR that determines the
%      replacement tolerance
%      default: 1
%   - .verbose: boolean for whether to print intermediate loss of
%      orthogonality (LOO) or relative residual (RelRes) per iteration
%      default: 0
%
% OUTPUTS:
% - QQ: n x m orthonormal matrix
% - RR: m x m upper triangular matrix
% - TT: m x m triangular matrix returned by routines such as BMGS_SVL and
%   BMGS_LTS; regardless of TT, it should hold that
%                      XX = QQ * RR.
% - run_time: total time elapsed measured by tic and toc; does not exclude
%   screen outputs from verbose
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
% Defaults
if nargin == 4
    param = param_init;
elseif nargin == 5
    param = param_init(param);
end

if isempty(musc)
    musc = 'houseqr';
end

switch lower(skel)
% Standard BCGS and reorthogonalized variants -----------------------------
    case {'bcgs'}
        tic;
        [QQ, RR] = bcgs(XX, s, musc, param);
        run_time = toc;

    case {'bcgs_ro'}
        tic;
        [QQ, RR] = bcgs_ro(XX, s, musc, param);
        run_time = toc;
        
    case {'bcgs_iro'}
        tic;
        [QQ, RR] = bcgs_iro(XX, s, musc, param);
        run_time = toc;

% [Stewart 2008] variant --------------------------------------------------
    case {'bcgs_sror'}
        tic;
        [QQ, RR] = bcgs_sror(XX, s, param.rpltol, param.verbose);
        run_time = toc;

% [Swirydowicz et al. 2020]/[Bielich et al. 2022] variant -----------------
    case {'bcgs_iro_ls'}
        tic;
        [QQ, RR] = bcgs_iro_ls(XX, s, musc, param);
        run_time = toc;

% "Roadmap" variants ------------------------------------------------------
    case {'bcgs_iro_3s'}
        tic;
        [QQ, RR] = bcgs_iro_3s(XX, s, musc, param);
        run_time = toc;

    case {'bcgs_iro_2s'}
        tic;
        [QQ, RR] = bcgs_iro_2s(XX, s, musc, param);
        run_time = toc;
        
    case {'bcgs_iro_1s'}
        tic;
        [QQ, RR] = bcgs_iro_1s(XX, s, musc, param);
        run_time = toc;

% Standard BMGS -----------------------------------------------------------
    case {'bmgs'}
        tic;
        [QQ, RR] = bmgs(XX, s, musc, param);
        run_time = toc;
        
% 3-sync BMGS variants ----------------------------------------------------
    case {'bmgs_svl'}
        tic;
        [QQ, RR, TT] = bmgs_svl(XX, s, musc, param);
        run_time = toc;
        
    case {'bmgs_lts'}
        tic;
        [QQ, RR, TT] = bmgs_lts(XX, s, musc, param);
        run_time = toc;
        
% 1-sync BMGS variants ----------------------------------------------------
    case {'bmgs_cwy'}
        tic;
        [QQ, RR, TT] = bmgs_cwy(XX, s, musc, param);
        run_time = toc;
        
    case {'bmgs_icwy'}
        tic;
        [QQ, RR, TT] = bmgs_icwy(XX, s, musc, param);
        run_time = toc;
        
% T-variants --------------------------------------------------------------
    case {'bcgs_iro_t'}
        tic;
        [QQ, RR, TT] = bcgs_iro_t(XX, s, musc, param);
        run_time = toc;
        
    case {'bmgs_t'}
        tic;
        [QQ, RR, TT] = bmgs_t(XX, s, musc, param);
        run_time = toc;
        
% Reorthogonalize first step (_f) -----------------------------------------
    case {'bcgs_iro_a'}
        tic;
        [QQ, RR] = bcgs_iro_a(XX, s, musc, param);
        run_time = toc;

    case {'bcgs_iro_a_3s'}
        tic;
        [QQ, RR] = bcgs_iro_a_3s(XX, s, musc, param);
        run_time = toc;

    case {'bcgs_iro_a_2s'}
        tic;
        [QQ, RR] = bcgs_iro_a_2s(XX, s, musc, param);
        run_time = toc;
        
    case {'bcgs_iro_a_1s'}
        tic;
        [QQ, RR] = bcgs_iro_a_1s(XX, s, musc, param);
        run_time = toc;
        
% P-variants --------------------------------------------------------------
    case {'bcgs_pio'}
        tic;
        [QQ, RR] = bcgs_pio(XX, s, musc, param);
        run_time = toc;
        
    case {'bcgs_pip'}
        tic;
        [QQ, RR] = bcgs_pip(XX, s, musc, param);
        run_time = toc;

% Reorthogonalized P-variants ---------------------------------------------
    case {'bcgs_pio_ro'}
        tic;
        [QQ, RR] = bcgs_pio_ro(XX, s, musc, param);
        run_time = toc;

    case {'bcgs_pio_iro'}
        tic;
        [QQ, RR] = bcgs_pio_iro(XX, s, musc, param);
        run_time = toc;

    case {'bcgs_pip_ro'}
        tic;
        [QQ, RR] = bcgs_pip_ro(XX, s, musc, param);
        run_time = toc;

    case {'bcgs_pip_iro'}
        tic;
        [QQ, RR] = bcgs_pip_iro(XX, s, musc, param);
        run_time = toc;
		
% Mixed Precision variants ------------------------------------------------
    case {'bcgs_iro_ls_mp'}
        tic;
        [QQ, RR] = bcgs_iro_ls_mp(XX, s, musc, param);
        run_time = toc;    
        
    case {'bcgs_pio_mp'}
        tic;
        [QQ, RR] = bcgs_pio_mp(XX, s, musc, param);
        run_time = toc;

    case {'bcgs_pio_ro_mp'}
        tic;
        [QQ, RR] = bcgs_pio_iro_mp(XX, s, musc, param);
        run_time = toc;
        
    case {'bcgs_pio_iro_mp'}
        tic;
        [QQ, RR] = bcgs_pio_iro_mp(XX, s, musc, param);
        run_time = toc;
        
    case {'bcgs_pip_mp'}
        tic;
        [QQ, RR] = bcgs_pip_mp(XX, s, musc, param);
        run_time = toc;  

    case {'bcgs_pip_ro_mp'}
        tic;
        [QQ, RR] = bcgs_pip_ro_mp(XX, s, musc, param);
        run_time = toc;
    
	case {'bcgs_pip_iro_mp'}
        tic;
        [QQ, RR] = bcgs_pip_iro_mp(XX, s, musc, param);
        run_time = toc;
        
    otherwise
        error('%s is not a viable skeleton option', skel);
end
if ~exist('TT', 'var')
    TT = eye(size(RR));
end
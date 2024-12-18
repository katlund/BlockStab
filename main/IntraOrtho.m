function [Q, R, T] = IntraOrtho(X, musc, param)
% [Q, R, T] = INTRAORTHO(X, musc, param) is a wrapper function for
% switching between intra-orthonormalization (column-vector-wise)
% algorithms.
% 
% INPUTS:
% - X: an n x s matrix, s <=n
% - musc is a char specifying the algorithm
% - param: a struct with the following optional fields:
%   - .chol: char specifying what type of Cholesky subroutine to call for
%      skeletons that hard-code Cholesky via a block Pythagorean trick
%      (e.g., BCGS_PIP, BCGS_PIO, BCGS_IRO_LS, BMGS_CWY, BMGS_ICWY, and
%      their reorthogonalized and multi-precision versions)
%      default: 'chol_nan'
%   - .rpltol: scalar argument for CGS_SROR that determines the
%      replacement tolerance
%      default: 1
%   - .verbose: boolean for whether to print intermediate loss of
%      orthogonality (LOO) or relative residual (RelRes) per iteration
%      default: 0
%
% OUTPUTS:
% - Q: m x s orthonormal matrix
% - R: s x s upper triangular matrix
% - T: s x s triangular matrix returned by routines such as MGS_SVL and
%   MGS_LTS; regardless of T, it should hold that
%                      X = Q * R.
%
% Part of [BlockStab](https://github.com/katlund/BlockStab) package.  Check README
% for how to properly cite and reuse this file.

%%
% Defaults
if nargin == 2
    param = param_init;
elseif nargin == 3
    param = param_init(param);
end

switch lower(musc)
    % CGS -----------------------------------------------------------------
    case {'cgs'}
        [Q, R] = cgs(X, param.verbose);
        
    % CGS with reorthogonalization ----------------------------------------
    case {'cgs_ro'}
        [Q, R] = cgs_ro(X, param.verbose);
        
    case {'cgs_iro'}
        [Q, R] = cgs_iro(X, param.verbose);
        
    case {'cgs_sro'}
        [Q, R] = cgs_sror(X, 0, param.verbose);    % rpltol = 0 ensures no replacement
        
    case {'cgs_sror'}
        if ~isfield(param, 'rpltol')
            param.rpltol = 10;
        end
        if isempty(param.rpltol)
            param.rpltol = 10;
        end
        [Q, R] = cgs_sror(X, param.rpltol, param.verbose);
        
    case {'cgs_iro_ls'}
        [Q, R] = cgs_iro_ls(X, param.verbose);

% CGS (P-variants) --------------------------------------------------------
    case {'cgs_p'}
        [Q, R] = cgs_p(X, param.verbose);
        
    case {'cgs_p_ro'}
        [Q, R] = cgs_p_ro(X, param.verbose);
        
    case {'cgs_p_iro'}
        [Q, R] = cgs_p_iro(X, param.verbose);
        
% MGS ---------------------------------------------------------------------
    case {'mgs'}
        [Q, R] = mgs(X, param.verbose);

    case {'mgs_ro'}
        [Q, R] = mgs_ro(X, param.verbose);
        
    case {'mgs_iro'}
        [Q, R] = mgs_iro(X, param.verbose);
        
% MGS (3-sync) ------------------------------------------------------------
    case {'mgs_svl'}
        [Q, R, T] = mgs_svl(X, param.verbose);
        
    case {'mgs_lts'}
        [Q, R, T] = mgs_lts(X, param.verbose);
        
% MGS (1-sync) ------------------------------------------------------------
    case {'mgs_icwy'}
        [Q, R, T] = mgs_icwy(X, param.verbose);
        
    case {'mgs_cwy'}
        [Q, R, T] = mgs_cwy(X, param.verbose);
        
% HouseQR -----------------------------------------------------------------
    case {'houseqr'}
        [Q, R] = qr(X,0);
        
% CholQR ------------------------------------------------------------------
    case {'cholqr'}
        [Q, R] = cholqr(X, param);
        
    case {'cholqr_ro'}
        [Q, R] = cholqr_ro(X, param);
        
    case {'iter_cholqr'}
        [Q, R] = iter_cholqr(X);
        
    case {'sh_cholqr_roro'}
        [Q, R] = sh_cholqr_roro(X, param);

    case {'cholqr_pinv'}
        [Q, R] = cholqr_pinv(X);

    case {'global', 'globalqr'}
        [Q, R] = globalqr(X);
        
    otherwise
        error('%s is not a viable muscle option', musc);
end
if ~exist('T', 'var')
    T = eye(size(R));
end
end
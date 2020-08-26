function [Q, R, T] = IntraOrtho(X, musc, rpltol)
% [Q, R, T] = INTRAORTHO(X, musc, rpltol) is a wrapper function for
% switching between intra-orthonormalization algorithms. musc is a char
% specifying the algorithm, and rpltol is an optional argument determining
% the replacement tolerance in CGS_SROR.
%
% For all possible muscle options, see the header of RUNTEST.

%%
addpath(genpath('muscles'))
musc = lower(musc);
switch musc        
    case {'cgs'}
        [Q, R] = cgs(X);
        
    case {'cgs_ro'}
        [Q1, R1] = cgs(X);
        [Q, R] = cgs(Q1);
        R = R*R1;
        
    case {'cgs_iro'}
        [Q, R] = cgs_iro(X);
        
    case {'cgs_sro'}
        [Q, R] = cgs_sror(X, 0);    % rpltol = 0 ensures no replacement
        
    case {'cgs_sror'}
        if nargin == 2
            rpltol = 10;
        end
        [Q, R] = cgs_sror(X, rpltol);
        
    case {'cgs_iro_ls'}
        [Q, R] = cgs_iro_ls(X);
        
    case {'cgs_iro_ls_v2'}
        [Q, R] = cgs_iro_ls_v2(X);
        
    case {'mgs'}
        [Q, R] = mgs(X);
        
    case {'mgs_ro'}
        [Q1, R1] = mgs(X);
        [Q, R] = mgs(Q1);
        R = R*R1;
        
    case {'mgs_iro'}
        [Q, R] = mgs_iro(X);
        
    case {'mgs_svl'}
        [Q, R, T] = mgs_svl(X);
        
    case {'mgs_lts'}
        [Q, R, T] = mgs_lts(X);
        
    case {'mgs_icwy'}
        [Q, R, T] = mgs_icwy(X);
        
    case {'mgs_cwy'}
        [Q, R, T] = mgs_cwy(X);
        
%--------------------------------------------------------------------------

    case {'mgs_svl_trans'}
        [Q, R, T] = mgs_svl_trans(X);
        
    case {'mgs_lts_trans'}
        [Q, R, T] = mgs_lts_trans(X);
        
    case {'mgs_icwy_trans'}
        [Q, R, T] = mgs_icwy_trans(X);  % this is identical to UTS!
        
    case {'mgs_cwy_trans'}
        [Q, R, T] = mgs_cwy_trans(X);
        
%--------------------------------------------------------------------------        
    case {'houseqr'}
        [Q, R] = qr(X,0);
        
    case {'cholqr'}
        [Q, R] = cholqr(X);
        
    case {'cholqr_ro'}
        [Q1, R1] = cholqr(X);
        [Q, R] = cholqr(Q1);
        R = R*R1;
        
    case {'iter_cholqr'}
        [Q, R] = iter_cholqr(X);
        
    case {'sh_cholqr_roro'}
        [Q, R] = sh_cholqr_roro(X);
        
    case {'svqb'}
        [Q, R] = svqb(X);
        
%--------------------------------------------------------------------------
    case {'global'}
        s = size(X,2);
        R = norm(X,'fro')/sqrt(s);  % with s-scaling
        Q = X/R;
        R = R*eye(s);
        
    case {'global-no-scale'}
        s = size(X,2);
        R = norm(X,'fro');
        Q = X/R;
        R = R*eye(s);
        
    otherwise
        error('%s is not a viable muscle option', musc);
end
if ~exist('T', 'var')
    T = eye(size(R));
end
end
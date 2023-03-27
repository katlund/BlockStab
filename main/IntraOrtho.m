function [Q, R, T] = IntraOrtho(X, musc, rpltol, verbose)
% [Q, R, T] = INTRAORTHO(X, musc, rpltol, verbose) is a wrapper function
% for switching between intra-orthonormalization algorithms. musc is a char
% specifying the algorithm, and rpltol is an optional argument determining
% the replacement tolerance in CGS_SROR. verbose is a Boolean for whether
% to print intermediate loss of orthogonality (LOO) or relative residual
% (RelRes) computations.

%%
% Defaults
if nargin == 2
    rpltol = [];
    verbose = 0;
elseif nargin == 3
    verbose = 0;
end

addpath(genpath('muscles'))
musc = lower(musc);
switch musc        
    case {'cgs'}
        [Q, R] = cgs(X, verbose);
        
    case {'cgs_p'}
        [Q, R] = cgs_p(X, verbose);
        
    case {'cgs_ro'}
        [Q1, R1] = cgs(X, verbose);
        [Q, R] = cgs(Q1, verbose);
        R = R*R1;
        
    case {'cgs_iro'}
        [Q, R] = cgs_iro(X, verbose);
        
    case {'cgs_sro'}
        [Q, R] = cgs_sror(X, 0, verbose);    % rpltol = 0 ensures no replacement
        
    case {'cgs_sror'}
        if nargin == 2
            rpltol = 10;
        end
        [Q, R] = cgs_sror(X, rpltol, verbose);
        
    case {'cgs_iro_ls'}
        [Q, R] = cgs_iro_ls(X, verbose);

    case {'cgs_iro_bl'}
        [Q, R] = cgs_iro_bl(X, verbose);
        
%--------------------------------------------------------------------------
    case {'mgs'}
        [Q, R] = mgs(X, verbose);
        
    case {'mgs_ro'}
        [Q1, R1] = mgs(X, verbose);
        [Q, R] = mgs(Q1, verbose);
        R = R*R1;
        
    case {'mgs_iro'}
        [Q, R] = mgs_iro(X, verbose);
        
%--------------------------------------------------------------------------
    case {'mgs_svl'}
        [Q, R, T] = mgs_svl(X, verbose);
        
    case {'mgs_lts'}
        [Q, R, T] = mgs_lts(X, verbose);
        
    case {'mgs_icwy'}
        [Q, R, T] = mgs_icwy(X, verbose);
        
    case {'mgs_cwy'}
        [Q, R, T] = mgs_cwy(X, verbose);
        
%--------------------------------------------------------------------------
    case {'houseqr'}
        [Q, R] = qr(X,0);
        
%--------------------------------------------------------------------------
    case {'cholqr'}
        [Q, R] = cholqr(X);
        
    case {'cholqr_pinv'}
        [Q, R] = cholqr_pinv(X);
        
    case {'cholqr_ro'}
        [Q1, R1] = cholqr(X);
        [Q, R] = cholqr(Q1);
        R = R*R1;
        
    case {'iter_cholqr'}
        [Q, R] = iter_cholqr(X);
        
    case {'sh_cholqr_roro'}
        [Q, R] = sh_cholqr_roro(X);
        
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
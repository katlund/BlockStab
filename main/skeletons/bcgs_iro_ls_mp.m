function [QQ, RR] = bcgs_iro_ls_mp(XX, s, musc, verbose, param)
% [QQ, RR] = BCGS_IRO_LS_MP(XX, s, musc, verbose) performs Block Classical
% Gram-Schmidt with Reorthogonalization and a Low-Sync formulation
% (actually, a one-sync) on the n x m matrix XX. It is a block
% generalization of CGS_IRO_LS, i.e., Algorithm 3 from [Swirydowicz, et.
% al. 2020] or Algorithm 2 from [Bielich, et al. 2022].  Note that no
% muscle is explicitly required, because CholQR is hard-coded for all
% intra-orthogonalizations; it can, however, be passed to InnerProd.
%
% This mixed precision version furthermore selectively uses simulated quad
% precision via Advanpix or the Symbolic Math Toolbox as specified by
% param.mp, with returned quantities cast to double.
%
% See BGS and MP_SWITCH for more details about the parameters.
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
% Default: debugging off
if nargin < 4
    verbose = 0;
    param.mp_package = 'advanpix';
elseif nargin < 5
    param.mp_package = 'advanpix';
end

% Set up quad-precision subroutine
switch param.mp_package
    case 'advanpix'
        param.mp_digits = 34;
    case 'symbolic math'
        param.mp_digits = 32;
end
qp = @(x) mp_switch(x, param);

% Pre-allocate memory for QQ and RR
[m, n] = size(XX);
QQ = qp(zeros(m,n));
RR = zeros(n,n);
p = n/s;

% Set up block indices
kk = 1:s;
s1 = 1:s;
s2 = s+1:2*s;

% Initialize Q_tmp
Q_tmp = qp(XX(:,kk));

if verbose
    fprintf('         LOO      |    RelRes\n');
    fprintf('-----------------------------------\n');
end

for k = 2:p
    % Increment block index
    kk = kk + s;
    sk = s*k;
    
    % Pull out block vector (keeps MATLAB from copying full X repeatedly)
    U = qp(XX(:,kk));
    
    % Compute temporary quantities (inner product) in quad precision.  Note
    % that for quantities stored in quad, MATLAB continues to compute in
    % quad until explicitly told otherwise.
    if k == 2
        R_tmp = InnerProd(Q_tmp, [Q_tmp U], musc);
    else
        tmp = InnerProd([QQ(:, 1:sk-2*s) Q_tmp], [Q_tmp U], musc);
        W = tmp(1:sk-2*s, s1);
        Z = tmp(1:sk-2*s, s2);
        R_tmp = tmp(end-s+1:end,:) - W' * [W Z];
    end
    
    % Pythagorean trick for RR diagonals; compute in quad
    R_diag = chol_free_mp(R_tmp(:, s1), param);
    
    % Assign finished entries of RR; cast back to double
    RR(kk-s, kk-s) = double(R_diag);
    RR(kk-s, kk) = double(R_diag'\ R_tmp(:, s2)); 
    
    if k == 2
        % Finish normalizing QQ(:,k-1)
        QQ(:,kk-s) = Q_tmp / R_diag; % quad precision
    else
        % Assign finished entries of RR and cast back to double
        RR(1:sk-2*s, kk-s) = RR(1:sk-2*s, kk-s) + double(W);
        RR(1:sk-2*s, kk) = double(Z);
        
        % Finish normalizing QQ(:,k-1)
        QQ(:,kk-s) = qp((Q_tmp - QQ(:, 1:sk-2*s) * W)) / R_diag;
    end
    
    % Set up temporary block vector for next iteration
    Q_tmp = qp(U - QQ(:, 1:sk-s) * RR(1:sk-s, kk));
    
    if verbose
        fprintf('%3.0d:', k-1);
        fprintf('  %2.4e  |',...
            norm( eye(sk-s) - InnerProd(double(QQ(:, 1:sk-s)), double(QQ(:, 1:sk-s)), musc) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk-s) - double(QQ(:,1:sk-s)) * RR(1:sk-s,1:sk-s) ) / norm(XX(:,1:sk-s)) );
    end
end

% Finish renormalizing last basis vector and assign last diagonal entry of
% RR.  Note that this requires just one more sync, no IntraOrtho needed.
tmp = InnerProd([QQ(:,1:n-s) Q_tmp], Q_tmp, musc);
W = tmp(1:n-s,:);
R_tmp = tmp(end-s+1:end,:) - W' * W;
R_diag = chol_free_mp(R_tmp);
RR(kk, kk) = double(R_diag);
RR(1:n-s, kk) = RR(1:n-s, kk) + double(W);
QQ(:,kk) = qp((Q_tmp - QQ(:,1:n-s) * W)) / R_diag;

% Cast QQ back to double
QQ = double(QQ);

if verbose
    fprintf('%3.0d:', k+1);
    fprintf('  %2.4e  |', norm( eye(n) - QQ' * QQ ) );
    fprintf('  %2.4e\n', norm( XX - QQ * RR ) / norm(XX) );
end
end

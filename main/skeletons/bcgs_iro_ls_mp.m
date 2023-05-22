function [QQ, RR] = bcgs_iro_ls_mp(XX, s, musc, param)
% [QQ, RR] = BCGS_IRO_LS_MP(XX, s, musc, param) performs Block Classical
% Gram-Schmidt with Reorthogonalization and a Low-Sync formulation
% (actually, a one-sync) on the n x m matrix XX. It is a block
% generalization of CGS_IRO_LS, i.e., Algorithm 3 from [Swirydowicz, et.
% al. 2020] or Algorithm 2 from [Bielich, et al. 2022].  Note that no
% muscle is explicitly required, because CholQR is hard-coded for all
% intra-orthogonalizations; it can, however, be passed to InnerProd.
%
% This mixed precision version computes the inputs to Cholesky and the
% Cholesky factorization itself in simulated quadruple (or other,
% user-specified precision) precision.  See MP_SWITCH for details on the
% param struct.
%
% See BGS and MP_SWITCH for more details about the parameters, and
% INTRAORTHO for musc options.
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
% Defaults
if nargin < 4
    param = mp_param_init;
elseif nargin == 4
    param = mp_param_init(param);
end

% Set up quad-precision subroutine
qp = @(x) mp_switch(x, param);

% Pre-allocate memory for QQ and RR
[m, n] = size(XX);
RR = zeros(n,n);
QQ = zeros(m,n);
p = n/s;

% Set up block indices
kk = 1:s;
s1 = 1:s;
s2 = s+1:2*s;

% Initialize
U = qp(XX(:,kk));

if param.verbose
    fprintf('         LOO      |    RelRes\n');
    fprintf('-----------------------------------\n');
end

for k = 2:p
    % Increment block index
    kk = kk + s;
    sk = s*k;
    
    % Pull out block vector (keeps MATLAB from copying full X repeatedly)
    W = XX(:,kk);
    
    % Compute temporary quantities (inner product) in qp.  Note that
    % standard operations are overloaded.
    if k == 2
        R_tmp = InnerProd(qp(U), qp([U W]), musc); % returned in qp
    else
        tmp = InnerProd(qp([QQ(:, 1:sk-2*s) U]), qp([U W]), musc);  % returned in qp
        Y = tmp(1:sk-2*s, s1);  % returned in qp
        Z = tmp(1:sk-2*s, s2);  % returned in qp
        R_tmp = tmp(end-s+1:end,:) - Y' * [Y Z];  % returned in qp
    end
    
    % Pythagorean trick for RR diagonals; R_tmp is already in qp from
    % previous step, and R_diag is returned in qp, so no need to recast it
    R_diag = chol_switch(R_tmp(:, s1), param);
    
    % Assign finished entries of RR; cast to double
    RR(kk-s, kk-s) = double(R_diag);
    RR(kk-s, kk) = double(R_diag'\ R_tmp(:, s2));
    
    if k == 2
        % Finish normalizing QQ(:,k-1) and cast to double
        QQ(:,kk-s) = double(qp(U) / R_diag);
    else
        % Assign finished entries of RR and cast to double
        RR(1:sk-2*s, kk-s) = RR(1:sk-2*s, kk-s) + double(Y);
        RR(1:sk-2*s, kk) = double(Z);
        
        % Finish normalizing QQ(:,k-1)
        QQ(:,kk-s) = double(qp(U - QQ(:, 1:sk-2*s) * Y) / R_diag);
    end
    
    % Set up temporary block vector for next iteration
    U = W - QQ(:, 1:sk-s) * RR(1:sk-s, kk);
    
    if param.verbose
        fprintf('%3.0d:', k-1);
        fprintf('  %2.4e  |',...
            norm( eye(sk-s) - InnerProd(QQ(:, 1:sk-s), QQ(:, 1:sk-s), musc) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk-s) - QQ(:,1:sk-s) * RR(1:sk-s,1:sk-s) ) / norm(XX(:,1:sk-s)) );
    end
end

% Finish renormalizing last basis vector and assign last diagonal entry of
% RR.  Note that this requires just one more sync, no IntraOrtho needed.
tmp = InnerProd(qp([QQ(:,1:n-s) U]), qp(U), musc);  % returned in qp
Y = tmp(1:n-s,:);  % returned in qp
R_diag = chol_switch(tmp(end-s+1:end,:) - Y' * Y, param);  % returned in qp
RR(kk, kk) = double(R_diag);
RR(1:n-s, kk) = RR(1:n-s, kk) + double(Y);
QQ(:,kk) = double(qp(U - QQ(:,1:n-s) * Y) / R_diag);

if param.verbose
    fprintf('%3.0d:', k);
    fprintf('  %2.4e  |', norm( eye(n) - QQ' * QQ ) );
    fprintf('  %2.4e\n', norm( XX - QQ * RR ) / norm(XX) );
end
end

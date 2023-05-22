function [QQ, RR] = bcgs_pio_mp(XX, s, musc, param)
% [QQ, RR] = BCGS_PIO(XX, s, musc, param) performs Block Classical
% Gram-Schmidt with Pythagorean Intra-Orthogonalization modification on the
% m x n matrix XX with p = n/s block partitions each of size s with
% intra-orthogonalization procedure determined by musc.  BCGS_PIO is a
% block generalization of CGS-P/Algorithm 2 from [Smoktunowicz et. al.
% 2006] derived in [Carson et al. 2021].
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
sk = s;

% Pull out block vector (keeps MATLAB from copying full X repeatedly)
W = XX(:,kk);

% First step
[QQ(:,kk), RR(kk,kk)] = IntraOrtho(W, musc, param);

if param.verbose
    fprintf('         LOO      |    RelRes\n');
    fprintf('-----------------------------------\n');
    fprintf('%3.0d:', 1);
    fprintf('  %2.4e  |',...
        norm( eye(s) - InnerProd(QQ(:, 1:s), QQ(:, 1:s), musc) ) );
    fprintf('  %2.4e\n',...
        norm( XX(:,1:s) - QQ(:,1:s) * RR(1:s,1:s) ) / norm(XX(:,1:s)) );
end

for k = 2:p
    % Update block indices
    kk = kk + s;
    sk = sk + s;
    
    % Pull out block vector (keeps MATLAB from copying full X repeatedly)
    W = XX(:,kk);

    % Project W
    R_col = InnerProd(QQ(:,1:sk-s), W, musc);

    % Batched IntraOrtho; cast to qp; note that standard operations are
    % overloaded
    [~, tmp] = qp(IntraOrtho([W zeros(size(W)); zeros(sk-s,s) R_col],...
        musc, param)); % returned in qp
    tmp = tmp' * tmp; % returned in qp

    % Compute Cholesky in qp
    R_diag = chol_switch(tmp(1:s,1:s) - tmp(end-s+1:end, end-s+1:end), param); % returned in qp
    
    % Assign finished entries of RR; cast to double
    RR(kk,kk) = double(R_diag);
    RR(1:sk-s,kk) = R_col;
    
    % Compute next basis vector    
    QQ(:,kk) = double(qp(W - QQ(:,1:sk-s) * R_col) / R_diag);
    
    if param.verbose
        fprintf('%3.0d:', k);
        fprintf('  %2.4e  |',...
            norm( eye(sk) - InnerProd(QQ(:, 1:sk), QQ(:, 1:sk), musc) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk) - QQ(:,1:sk) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
    end
end
end
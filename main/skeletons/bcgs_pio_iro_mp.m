function [QQ, RR] = bcgs_pio_iro_mp(XX, s, musc, param)
% [QQ, RR] = BCGS_PIO_IRO_MP(XX, s, musc, param) performs BCGS_PIO with
% Inner ReOrthonormalization on the m x n matrix XX with p = n/s block
% partitions each of size s and with intra-orthogonalization procedure
% determined by musc.
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
    
    %% First step
    % Project W
    RR1 = InnerProd(QQ(:,1:sk-s), W, musc);

    % Batched IntraOrtho; cast to qp; note that standard operations are
    % overloaded
    [~, tmp] = IntraOrtho([W zeros(size(W)); zeros(sk-s, s) RR1],...
        musc, param); % returned in double
    tmp = qp(tmp)' * qp(tmp); % returned in qp

    % Compute Cholesky in qp
    R1 = chol_switch(tmp(1:s,1:s) - tmp(end-s+1:end, end-s+1:end), param); % returned in qp

    % Compute intermediate basis vector; cast to double
    W = double(qp(W - QQ(:,1:sk-s) * RR1) / R1);
    
    %% Second step
    RR2 = InnerProd(QQ(:,1:sk-s), W, musc);
    [~, tmp] = IntraOrtho([W zeros(size(W)); zeros(sk-s, s) RR2],...
        musc, param);  % returned in double
    tmp = qp(tmp)' * qp(tmp); % returned in qp

    % Compute Cholesky in qp
    R2 = chol_switch(tmp(1:s,1:s) - tmp(end-s+1:end, end-s+1:end), param); % returned in qp

    % Compute next basis vector; cast to double
    QQ(:,kk) = double(qp(W - QQ(:,1:sk-s) * RR2) / R2);
    
    % Assign RR in double and combine both steps
    RR(kk,kk) = double(R2) * double(R1);
    RR(1:sk-s,kk) = RR1 + RR2 * double(R1);
    
    if param.verbose
        fprintf('%3.0d:', k);
        fprintf('  %2.4e  |',...
            norm( eye(sk) - InnerProd(QQ(:, 1:sk), QQ(:, 1:sk), musc) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk) - QQ(:,1:sk) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
    end
end
end
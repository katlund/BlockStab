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
% See BGS for more details about the parameters, and INTRAORTHO for musc
% options.
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
% Defaults
if nargin < 4
    param.verbose = 0;
    param.mp_package = 'advanpix';
end
if ~isfield(param, 'chol')
    param.chol = 'chol_free';
else
    if isempty(param.chol)
        param.chol = 'chol_free';
    end
end

% Set up quad-precision subroutine
qp = @(x) mp_switch(x, param);

% Pre-allocate memory for QQ and RR
[m, n] = size(XX);
RR = zeros(n,n);
QQ = qp(zeros(m,n));
p = n/s;

% Set up block indices
kk = 1:s;
sk = s;

W = qp(XX(:,kk));
[QQ(:,kk), R1] = IntraOrtho(W, musc);
RR(kk,kk) = double(R1);

if param.verbose
    fprintf('         LOO      |    RelRes\n');
    fprintf('-----------------------------------\n');
    fprintf('%3.0d:', 1);
    fprintf('  %2.4e  |',...
        norm( eye(s) - InnerProd(double(QQ(:, 1:s)), double(QQ(:, 1:s)), musc) ) );
    fprintf('  %2.4e\n',...
        norm( XX(:,1:s) - double(QQ(:,1:s)) * RR(1:s,1:s) ) / norm(XX(:,1:s)) );
end

for k = 2:p
    % Update block indices
    kk = kk + s;
    sk = sk + s;
    
    % Extract next block vector in quad storage
    W = qp(XX(:,kk));
    
    %% First step
    RR1 = InnerProd(QQ(:,1:sk-s), W, musc);
    [~, tmp] = IntraOrtho([W zeros(size(W)); zeros(sk-s, s) RR1], musc);
    tmp = tmp' * tmp;

    % Compute Cholesky in quad
    R1 = chol_switch(tmp(1:s,1:s) - tmp(end-s+1:end, end-s+1:end));

    % Compute intermediate basis vector
    W = ( W - QQ(:,1:sk-s) * RR1 ) / R1;
    
    %% Second step
    RR2 = InnerProd(QQ(:,1:sk-s), W, musc);
    [~, tmp] = IntraOrtho([W zeros(size(W)); zeros(sk-s, s) RR2], musc); 
    tmp = tmp' * tmp;

    % Compute Cholesky in quad
    R2 = chol_switch(tmp(1:s,1:s) - tmp(end-s+1:end, end-s+1:end));

    % Compute next basis vector
    QQ(:,kk) = ( W - QQ(:,1:sk-s) * tmp(1:sk-s,:) ) / R2;
    
    % Assign RR in double and combine both steps
    RR(kk,kk) = double(R2) * double(R1);
    RR(1:sk-s,kk) = double(RR1) + double(RR2) * double(RR1);
    
    if param.verbose
        fprintf('%3.0d:', k);
        fprintf('  %2.4e  |',...
            norm( eye(sk) - InnerProd(double(QQ(:, 1:sk)), double(QQ(:, 1:sk)), musc) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk) - double(QQ(:,1:sk)) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
    end
end
end
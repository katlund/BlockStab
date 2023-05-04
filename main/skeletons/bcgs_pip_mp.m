function [QQ, RR] = bcgs_pip_mp(XX, s, musc, verbose, param)
% [QQ, RR] = BCGS_PIP_MP(XX, s, musc, verbose, param) performs Block
% Classical Gram-Schmidt with Pythagorean Inner Product modification on the
% m x n matrix XX with p = n/s block partitions each of size s with
% intra-orthogonalization procedure determined by musc. BCGS_PIP is a block
% generalization of CGS-P/Algorithm 2 from [Smoktunowicz et. al. 2006]
% derived in [Carson et al. 2021].
%
% This mixed precision version computes the inputs to Cholesky and the
% Cholesky factorization itself in simulated quadruple precision.  See
% MP_SWITCH for details on the param struct.
%
% See BGS for more details about the parameters, and INTRAORTHO for musc
% options.
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
[QQ(:,kk), R_diag] = IntraOrtho(W, musc);
RR(kk,kk) = double(R_diag);

if verbose
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
    
    % Only sync point; in quad precision
    tmp = InnerProd([QQ(:,1:sk-s) W], W, musc);

    % Compute Cholesky in quad
    R_diag = chol_free_mp( tmp(kk,:) - tmp(1:sk-s,:)' * tmp(1:sk-s,:), param);
    
    % Assign RR in double
    RR(kk,kk) = double(R_diag); 
    RR(1:sk-s,kk) = double(tmp(1:sk-s,:));

    % Compute next basis vector
    QQ(:,kk) = ( W - QQ(:,1:sk-s) * tmp(1:sk-s,:) ) / R_diag;
    
    if verbose
        fprintf('%3.0d:', k);
        fprintf('  %2.4e  |',...
            norm( eye(sk) - InnerProd(double(QQ(:, 1:sk)), double(QQ(:, 1:sk)), musc) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk) - double(QQ(:,1:sk)) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
    end
end

% Cast QQ back to double
QQ = double(QQ);
end
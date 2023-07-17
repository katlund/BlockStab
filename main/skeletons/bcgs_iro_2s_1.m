function [QQ, RR] = bcgs_iro_2s_1(XX, s, musc, param)
% [QQ, RR] = BCGS_IRO_2S_1(XX, s, musc, param) performs Block Classical
% Gram-Schmidt with Inner ReOrthogonalization on the m x n matrix XX with p
% = n/s block partitions each of size s with intra-orthogonalization
% procedure determined by musc.  Skips the first norm of the first step and
% hard-codes CholQR as the IntraOrtho for columns k > 1 to reduce the sync
% points to 2 (hence '2s').
%
% See BGS for more details about the parameters, and INTRAORTHO for musc
% options.
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%% TODO: Update me with _1 explanation!

%%
% Default: debugging off
if nargin < 4
    param.verbose = 0;
end

% Pre-allocate memory for QQ and RR
[m, n] = size(XX);
RR = zeros(n,n);
QQ = zeros(m,n);
p = n/s;

% Set up block indices
kk = 1:s;
sk = s;

% Initial step
W = XX(:,kk);
[W, RR1] = IntraOrtho(W, musc, param);
[QQ(:,kk), RR(kk,kk)] = IntraOrtho(W, musc, param);   % reorthogonalize first step
RR(kk,kk) = RR(kk,kk) * RR1;

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

    % k.1
    S_col = InnerProd(QQ(:,1:sk-s), XX(:,kk), musc);
    W = XX(:,kk) - QQ(:,1:sk-s) * S_col;

    % k.2
    tmp = InnerProd([QQ(:,1:sk-s) W], W, musc);
    Y_col = tmp(1:sk-s,:);
    Y_diag = chol_switch(tmp(kk,:) - Y_col' * Y_col, param);
    QQ(:,kk) = (W - QQ(:,1:sk-s) * Y_col ) / Y_diag;
    
    % k.3
    RR(1:sk-s,kk) = S_col + Y_col;
    RR(kk,kk) = Y_diag;
    
    if param.verbose
        fprintf('%3.0d:', k);
        fprintf('  %2.4e  |',...
            norm( eye(sk) - InnerProd(QQ(:, 1:sk), QQ(:, 1:sk), musc) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk) - QQ(:,1:sk) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
    end
end
end
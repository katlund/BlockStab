function [QQ, RR] = bcgs_iro_2s(XX, s, musc, verbose)
% [QQ, RR] = BCGS_IRO_2S(XX, s, musc, verbose) performs Block Classical
% Gram-Schmidt with Inner ReOrthogonalization on the m x n matrix XX with p
% = n/s block partitions each of size s with intra-orthogonalization
% procedure determined by musc.  Skips the first norm of the first step and
% hard-codes CholQR as the IntraOrtho for columns k > 1 to reduce the sync
% points to 2 (hence '2s').
%
% See BGS for more details about the parameters, and INTRAORTHO for musc
% options.

%%
addpath(genpath('../'))

% Default: debugging off
if nargin < 4
    verbose = 0;
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
[QQ(:,kk), RR(kk,kk)] = IntraOrtho(XX(:,kk), musc);

if verbose
    fprintf('         LOO      |    RelRes\n');
    fprintf('-----------------------------------\n');
    fprintf('%3.0d:', 1);
    fprintf('  %2.4e  |',...
        norm( eye(s) - QQ(:, 1:s)' * QQ(:, 1:s) ) );
    fprintf('  %2.4e\n',...
        norm( XX(:,1:s) - QQ(:,1:s) * RR(1:s,1:s) ) / norm(XX(:,1:s)) );
end

for k = 2:p
    % Update block indices
    kk = kk + s;
    sk = sk + s;

    % k.1
    S_col = QQ(:,1:sk-s)' * XX(:,kk);
    W = XX(:,kk) - QQ(:,1:sk-s) * S_col;

    % k.2
    tmp = [QQ(:,1:sk-s) W]' * W;
    T_col = tmp(1:sk-s,:);
    diff = tmp(kk,:) - T_col' * T_col;
    T_diag = chol_nan(diff);
    QQ(:,kk) = (W - QQ(:,1:sk-s) * T_col ) / T_diag;
    
    % k.3
    RR(1:sk-s,kk) = S_col + T_col;
    RR(kk,kk) = T_diag;
    
    if verbose
        fprintf('%3.0d:', k);
        fprintf('  %2.4e  |',...
            norm( eye(sk) - QQ(:, 1:sk)' * QQ(:, 1:sk) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk) - QQ(:,1:sk) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
    end
end
end
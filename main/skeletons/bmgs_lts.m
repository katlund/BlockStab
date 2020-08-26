function [QQ, RR, TT] = bmgs_lts(XX, s, IOstr, verbose)
% [QQ, RR, TT] = BMGS_LTS(XX, s, IOstr, verbose) performs Block Modified
% Gram-Schmidt with Lower Triangular Solve on the m x n matrix XX with p =
% n/s block partitions each of size s with inner orthogonalization
% procedure determined by IOstr. BMGS_LTS is the block generalization of
% MGS_LTS.
%
% See BGS for more details about the parameters, and INTRAORTHO for IOstr
% options.

%%
addpath(genpath('../'))

% Default: debugging off
if nargin < 4
    verbose = 0;
end

% Pre-allocate memory for QQ, RR, and TT and auxiliary vectors
[m, n] = size(XX);
QQ = zeros(m,n);
RR = zeros(n,n);
TT = zeros(n,n);
p = n/s;

% Set up block indices
kk = 1:s;
sk = s;

W = XX(:,kk);
[QQ(:,kk), RR(kk,kk), TT(kk,kk)] = IntraOrtho(W, IOstr);

if verbose
    fprintf('         LOO      |    RelRes\n');
    fprintf('-----------------------------------\n');
    fprintf('%3.0d:', 1);
    fprintf('  %2.4e  |',...
        norm( eye(s) - QQ(:, 1:s)' * QQ(:, 1:s) ) );
    fprintf('  %2.4e\n',...
        norm( XX(:,1:s) - QQ(:,1:s) * RR(1:s,1:s) ) / norm(XX(:,1:s)) );
end

for k = 1:p-1
    % Update block indices
    kk = kk + s;
    
    W = XX(:,kk);
    
    RR(1:sk,kk) = TT(1:sk, 1:sk)' \ InnerProd(QQ(:,1:sk), W, IOstr);
    
    W = W - QQ(:,1:sk) * RR(1:sk,kk);
    
    [QQ(:,kk), RR(kk,kk), TT(kk,kk)] = IntraOrtho(W, IOstr);
    
    TT(1:sk,kk) = InnerProd(QQ(:,1:sk), QQ(:,kk), IOstr) * TT(kk,kk);
    
    sk = sk + s;
    if verbose
        fprintf('%3.0d:', k+1);
        fprintf('  %2.4e  |',...
            norm( eye(sk) - QQ(:, 1:sk)' * QQ(:, 1:sk) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk) - QQ(:,1:sk) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
    end
end
end
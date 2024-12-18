function [QQ, RR, TT] = bmgs_svl(XX, s, musc, param)
% [QQ, RR, TT] = BMGS_SVL(XX, s, musc, param) performs Block Modified
% Gram-Schmidt with the Schreiber-Van-Loan reformulation on the m x n
% matrix XX with p = n/s block partitions each of size s with inner
% orthogonalization procedure determined by musc.  When musc = 'HouseQR',
% this algorithm reproduces BMGS_H from [Barlow 2019], and when musc =
% 'MGS_SVL', this algorithm reproduces MGS3 from [Barlow 2019].  BMGS_SVL
% is the block generalization of MGS_SVL.
%
% See BGS for more details about the parameters, and INTRAORTHO for musc
% options.
%
% Part of [BlockStab](https://github.com/katlund/BlockStab) package.  Check README
% for how to properly cite and reuse this file.

%%
% Default: debugging off
if nargin < 4
    param.verbose = 0;
end

% Pre-allocate memory for QQ, RR, and TT
[m, n] = size(XX);
QQ = zeros(m,n);
RR = zeros(n,n);
TT = zeros(n,n);
p = n/s;

% Set up block indices
kk = 1:s;
sk = s;

W = XX(:,kk);
[QQ(:,kk), RR(kk,kk), TT(kk,kk)] = IntraOrtho(W, musc, param);

if param.verbose
    fprintf('         LOO      |    RelRes\n');
    fprintf('-----------------------------------\n');
    fprintf('%3.0d:', 1);
    fprintf('  %2.4e  |',...
        norm( eye(s) - InnerProd(QQ(:, 1:s), QQ(:, 1:s), musc) ) );
    fprintf('  %2.4e\n',...
        norm( XX(:,1:s) - QQ(:,1:s) * RR(1:s,1:s) ) / norm(XX(:,1:s)) );
end

for k = 1:p-1
    % Update block indices
    kk = kk + s;
    
    W = XX(:,kk);
    RR(1:sk,kk) = TT(1:sk, 1:sk)' * InnerProd(QQ(:,1:sk), W, musc);    
    W = W - QQ(:,1:sk) * RR(1:sk,kk);
    [QQ(:,kk), RR(kk,kk), TT(kk,kk)] = IntraOrtho(W, musc, param);
    TT(1:sk,kk) = -TT(1:sk,1:sk)*...
        InnerProd(QQ(:,1:sk), QQ(:,kk), musc) * TT(kk,kk);
    
    sk = sk + s;
    if param.verbose
        fprintf('%3.0d:', k+1);
        fprintf('  %2.4e  |',...
            norm( eye(sk) - InnerProd(QQ(:, 1:sk), QQ(:, 1:sk), musc) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk) - QQ(:,1:sk) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
    end
end
end
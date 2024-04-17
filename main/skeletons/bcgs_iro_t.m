function [QQ, RR, TT] = bcgs_iro_t(XX, s, musc, param)
% [QQ, RR, TT] = BCGS_IRO_T(XX, s, musc, param) performs the T-variant of
% BCGS_IRO on the m x n matrix XX with p = n/s block partitions each of
% size s and with intra-orthogonalization procedure determined by musc.
%
% See BGS for more details about the parameters, and INTRAORTHO for musc
% options.
%
% Part of [BlockStab](https://github.com/katlund) package.  Check README
% for how to properly cite and reuse this file.

%%
% Default: debugging off
if nargin < 4
    param.verbose = 0;
end

% Pre-allocate memory for QQ and RR
[m, n] = size(XX);
RR = zeros(n,n);
QQ = zeros(m,n);
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
    
    % First BCGS step
    RR1 = TT(1:sk,1:sk)' * InnerProd(QQ(:,1:sk), W, musc);
    W = W - QQ(:,1:sk) * TT(1:sk,1:sk) * RR1;
    [W, R1, T1] = IntraOrtho(W, musc, param);
    
    % Second BCGS step
    RR(1:sk,kk) = TT(1:sk,1:sk)' * InnerProd(QQ(:,1:sk), W, musc);
    W = W - QQ(:,1:sk) * TT(1:sk,1:sk) * RR(1:sk,kk);
    [QQ(:,kk), RR(kk,kk), TT(kk,kk)] = IntraOrtho(W, musc, param);
    
    % Combine both steps
    RR(1:sk,kk) = RR1 + RR(1:sk,kk) * R1;
    RR(kk,kk) = RR(kk,kk) * R1;
    TT(kk,kk) = TT(kk,kk) * T1;
    
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
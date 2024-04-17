function [QQ, RR] = bmgs(XX, s, musc, param)
% [QQ, RR] = BMGS(XX, s, musc, param) performs Block Modified
% Gram-Schmidt on the m x n matrix XX with p = n/s block partitions each of
% size s with inner orthogonalization procedure determined by musc.  BMGS
% is the block generalization of MGS.
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
p = n/s;

% Set up block indices
kk = 1:s;

W = XX(:,kk);
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

for k = 1:p-1
    % Update block indices
    kk = kk + s;
    jj = 1:s;
    
    % Pull out vectors (keeps MATLAB from copying full X repeatedly)
    W = XX(:,kk);
    
    for j = 1:k
        RR(jj,kk) = InnerProd(QQ(:,jj), W, musc);
        W = W - QQ(:,jj)*RR(jj,kk);
        jj = jj + s;
    end
    [QQ(:,kk), RR(kk,kk)] = IntraOrtho(W, musc, param);
    
    if param.verbose
        fprintf('%3.0d:', k+1);
        sk = s*(k+1);
        fprintf('  %2.4e  |',...
            norm( eye(sk) - InnerProd(QQ(:, 1:sk), QQ(:, 1:sk), musc) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk) - QQ(:,1:sk) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
    end
end
end
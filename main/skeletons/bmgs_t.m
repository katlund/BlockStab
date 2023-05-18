function [QQ, RR, TT] = bmgs_t(XX, s, musc, param)
% [QQ, RR, TT] = BMGS_T(XX, s, musc, param) performs the T-variant of BMGS
% on the m x n matrix XX with p = n/s block partitions each of size s with
% inner orthogonalization procedure determined by musc.
%
% See BGS for more details about the parameters, and INTRAORTHO for musc
% options.
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

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
    jj = 1:s;
    
    % Pull out vectors (keeps MATLAB from copying full X repeatedly)
    W = XX(:,kk);
    
    for j = 1:k
        RR(jj,kk) = TT(jj,jj)' * InnerProd(QQ(:,jj), W, musc);

        W = W - QQ(:,jj) * TT(jj,jj) * RR(jj,kk);
        jj = jj + s;
    end
    [QQ(:,kk), RR(kk,kk), TT(kk,kk)] = IntraOrtho(W, musc, param);
    
    if param.verbose
        fprintf('%3.0d:', k+1);
        sk = s*(k+1);
        fprintf('  %2.4e  |',...
            norm( eye(sk) - InnerProd(QQ(:, 1:sk),  QQ(:, 1:sk), musc) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk) - QQ(:,1:sk) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
    end
end
end
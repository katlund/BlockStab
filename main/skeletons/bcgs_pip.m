function [QQ, RR] = bcgs_pip(XX, s, musc, verbose)
% [QQ, RR] = BCGS_PIP(XX, s, musc, verbose) performs Block Classical
% Gram-Schmidt with Pythagorean Inner Product modification on the m x n
% matrix XX with p = n/s block partitions each of size s with
% intra-orthogonalization procedure determined by musc. BCGS_PIP is a
% block generalization of CGS-P/Algorithm 2 from [Smoktunowicz et. al.
% 2006] derived in [Carson et al. 2021].
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
end

% Pre-allocate memory for QQ and RR
[m, n] = size(XX);
RR = zeros(n,n);
QQ = zeros(m,n);
p = n/s;

% Set up block indices
kk = 1:s;
sk = s;

W = XX(:,kk);
[QQ(:,kk), RR(kk,kk)] = IntraOrtho(W, musc);

if verbose
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
    
    S = InnerProd([QQ(:,1:sk) W], W, musc);
    Y = S(1:sk,:);
    RR(kk,kk) = chol_nan(S(kk,:) - Y'*Y);
    W = W - QQ(:,1:sk) * Y;
    
    RR(1:sk,kk) = Y;
    QQ(:,kk) = W / RR(kk,kk);
    
    sk = sk + s;
    if verbose
        fprintf('%3.0d:', k+1);
        fprintf('  %2.4e  |',...
            norm( eye(sk) - InnerProd(QQ(:, 1:sk), QQ(:, 1:sk), musc) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk) - QQ(:,1:sk) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
    end
end
end
function [QQ, RR] = bcgs_pio(XX, s, musc, param)
% [QQ, RR] = BCGS_PIO(XX, s, musc, param) performs Block Classical
% Gram-Schmidt with Pythagorean Intra-Orthogonalization modification on the
% m x n matrix XX with p = n/s block partitions each of size s with
% intra-orthogonalization procedure determined by musc.  BCGS_PIO is a
% block generalization of CGS-P/Algorithm 2 from [Smoktunowicz et. al.
% 2006] derived in [Carson et al. 2021].
%
% See BGS for more details about the parameters, and INTRAORTHO for musc
% options.

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

for k = 2:p
    % Update block indices
    kk = kk + s;
    sk = sk + s;
    
    % Set up next vector
    W = XX(:,kk);

    % Sync points
    RR(1:sk-s,kk) = InnerProd(QQ(:,1:sk-s), W, musc);
    [~, tmp] = IntraOrtho([W zeros(size(W)); zeros(sk-s, s) RR(1:sk-s,kk)],...
        musc, param);
    tmp = tmp' * tmp;
    RR(kk,kk) = chol_switch( tmp(1:s,1:s) - tmp(end-s+1:end, end-s+1:end), param);
    QQ(:,kk) = ( W - QQ(:,1:sk-s) * RR(1:sk-s,kk) ) / RR(kk,kk);
    
    if param.verbose
        fprintf('%3.0d:', k);
        fprintf('  %2.4e  |',...
            norm( eye(sk) - InnerProd(QQ(:, 1:sk), QQ(:, 1:sk), musc) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk) - QQ(:,1:sk) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
    end
end
end
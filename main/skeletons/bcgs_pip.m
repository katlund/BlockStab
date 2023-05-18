function [QQ, RR] = bcgs_pip(XX, s, musc, param)
% [QQ, RR] = BCGS_PIP(XX, s, musc, param) performs Block Classical
% Gram-Schmidt with Pythagorean Inner Product modification on the m x n
% matrix XX with p = n/s block partitions each of size s with
% intra-orthogonalization procedure determined by musc. BCGS_PIP is a block
% generalization of CGS-P/Algorithm 2 from [Smoktunowicz et. al. 2006]
% derived in [Carson et al. 2021].
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
   
    % Extract next block vector
    W = XX(:,kk);
    
    % Only sync point
    tmp = InnerProd([QQ(:,1:sk-s) W], W, musc);
    
    % Assign RR
    RR(kk,kk) = chol_switch(tmp(kk,:) - RR(1:sk-s,kk)'* RR(1:sk-s,kk), param);
    RR(1:sk,kk) = tmp(1:sk,:);

    % Compute next basis vector
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
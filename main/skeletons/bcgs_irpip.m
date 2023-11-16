function [QQ, RR] = bcgs_irpip(XX, s, musc, param)
% BCGS+BCGS-PIP in each inner loop.

%%
% Default: debugging off
if nargin < 4
    param.verbose = 0;
end
% param.verbose = 1;
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

for k = 1:p-1
    % Update block indices
    kk = kk + s;

    W = XX(:,kk);

    % BCGS step
    % First BGS step
    RR1 = InnerProd(QQ(:,1:sk), W, musc);    
    W = W - QQ(:,1:sk) * RR1;
    [W, R1] = IntraOrtho(W, musc, param);

    % Second BGS-PIP step
    tmp = InnerProd([QQ(:,1:sk) W], W, musc);
    RR(1:sk,kk) = tmp(1:sk,:);
    RR(kk,kk) = chol_switch(tmp(kk,:) - RR(1:sk,kk)' * RR(1:sk,kk), param);
    QQ(:,kk) = ( W - QQ(:,1:sk) * RR(1:sk,kk) ) / RR(kk,kk);
    RR(1:sk,kk) = RR1 + RR(1:sk,kk) * R1;
    RR(kk,kk) = RR(kk,kk) * R1;

    
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
function [QQ, RR] = bcgs_pio_iro(XX, s, musc, param)
% [QQ, RR] = BCGS_PIO_IRO(XX, s, musc, param) performs BCGS_PIO with
% Inner ReOrthonormalization on the m x n matrix XX with p = n/s block
% partitions each of size s and with intra-orthogonalization procedure
% determined by musc.
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
    
    % First step
    RR1 = InnerProd(QQ(:,1:sk-s), W, musc);
    [~, tmp] = IntraOrtho([W zeros(size(W)); zeros(sk-s, s) RR1],...
        musc, param);
    tmp = tmp' * tmp;
    R1 = chol_switch(tmp(1:s,1:s) - tmp(end-s+1:end, end-s+1:end), param);
    W = ( W - QQ(:,1:sk-s) * RR1 ) / R1;
    
    % Second step
    RR(1:sk-s,kk) = InnerProd(QQ(:,1:sk-s), W, musc);
    [~, tmp] = IntraOrtho([W zeros(size(W)); zeros(sk-s, s) RR(1:sk-s,kk)],...
        musc, param);
    tmp = tmp' * tmp;
    RR(kk,kk) = chol_switch(tmp(1:s,1:s) - tmp(end-s+1:end, end-s+1:end), param);
    QQ(:,kk) = ( W - QQ(:,1:sk-s) * RR(1:sk-s,kk) ) / RR(kk,kk);
    
    % Combine both steps
    RR(1:sk-s,kk) = RR1 + RR(1:sk-s,kk) * R1;
    RR(kk,kk) = RR(kk,kk) * R1;
    
    if param.verbose
        fprintf('%3.0d:', k);
        fprintf('  %2.4e  |',...
            norm( eye(sk) - InnerProd(QQ(:, 1:sk), QQ(:, 1:sk), musc) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk) - QQ(:,1:sk) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
    end
end
end
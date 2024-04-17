function [QQ, RR] = bcgs_iro_1s(XX, s, musc, param)
% [QQ, RR] = BCGS_IRO_1S(XX, s, musc, param) performs Block Classical
% Gram-Schmidt with Inner ReOrthogonalization with 1 sync on the m x n
% matrix XX with p = n/s block partitions each of size s with
% intra-orthogonalization procedure determined by musc.  See BCGS_IRO_LS
% for an alternative formulation based on Algorithm 3 from [Swirydowicz,
% et. al. 2020] or Algorithm 2 from [Bielich, et al. 2022].
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
s1 = kk;
s2 = s1 + s;

% Initial step
[QQ(:,kk), RR(kk,kk)] = IntraOrtho(XX(:,kk), musc, param);

if param.verbose
    fprintf('         LOO      |    RelRes\n');
    fprintf('-----------------------------------\n');
    fprintf('%3.0d:', 1);
    fprintf('  %2.4e  |',...
        norm( eye(s) - InnerProd(QQ(:, 1:s), QQ(:, 1:s), musc) ) );
    fprintf('  %2.4e\n',...
        norm( XX(:,1:s) - QQ(:,1:s) * RR(1:s,1:s) ) / norm(XX(:,1:s)) );
end

% k.1.1
S_col = InnerProd(QQ(:,1:sk), XX(:,kk+s), musc);

% Update block indices
kk = kk + s;
sk = sk + s;

% k.1.2
W = XX(:,kk) - QQ(:,1:sk-s) * S_col;

for k = 2:p-1
    % k.2
    tmp = InnerProd([QQ(:,1:sk-s) W], [W XX(:,kk+s)], musc);
    Y_col = tmp(1:sk-s,s1);
    Y_diag = chol_switch(tmp(kk,s1) - Y_col' * Y_col, param);
    QQ(:,kk) = (W - QQ(:,1:sk-s) * Y_col ) / Y_diag;
    
    % k.3
    RR(1:sk-s,kk) = S_col + Y_col;
    RR(kk,kk) = Y_diag;
    
    if param.verbose
        fprintf('%3.0d:', k);
        fprintf('  %2.4e  |',...
            norm( eye(sk) - InnerProd(QQ(:, 1:sk), QQ(:, 1:sk), musc) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk) - QQ(:,1:sk) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
    end

    % (k+1).1.1
    Z_col = tmp(1:sk-s,s2);
    P = tmp(kk,s2);
    S_col = [Z_col; Y_diag' \ (P - Y_col' * Z_col)];

    % Update block indices
    kk = kk + s;
    sk = sk + s;

    % (k+1).1.2
    W = XX(:,kk) - QQ(:,1:sk-s) * S_col;
end

% p.2
tmp = InnerProd([QQ(:,1:sk-s) W], W, musc);
Y_col = tmp(1:sk-s,s1);
Y_diag = chol_switch(tmp(kk,s1) - Y_col' * Y_col, param);
QQ(:,kk) = (W - QQ(:,1:sk-s) * Y_col ) / Y_diag;

% p.3
RR(1:sk-s,kk) = S_col + Y_col;
RR(kk,kk) = Y_diag;

if param.verbose
    fprintf('%3.0d:', k);
    fprintf('  %2.4e  |',...
        norm( eye(sk) - InnerProd(QQ(:, 1:sk), QQ(:, 1:sk), musc) ) );
    fprintf('  %2.4e\n',...
        norm( XX(:,1:sk) - QQ(:,1:sk) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
end
end
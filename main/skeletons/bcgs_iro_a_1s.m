function [QQ, RR] = bcgs_iro_a_1s(XX, s, musc, param)
% [QQ, RR] = BCGS_IRO_A_1S(XX, s, musc, param) performs a reformulated
% version of BCGS_IRO_A with 1 sync point.  Despite being a _a algorithm,
% only one muscle is allowed, namely IO_A in the very first step.
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

% Extract W
W = XX(:,kk);

% IO_A
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
function [QQ, RR] = bcgs_iro_f_1s(XX, s, musc, param)
% [QQ, RR] = BCGS_IRO_F_1S(XX, s, musc, param) performs BCGS_IRO_1s with
% the first (_f) block vector reorthogonalized.
%
% See BGS for more details about the parameters, and INTRAORTHO for musc
% options.
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%% TODO: Update me with _1 explanation!

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
W = XX(:,kk);
[W, RR1] = IntraOrtho(W, musc, param);
[QQ(:,kk), RR(kk,kk)] = IntraOrtho(W, musc, param);   % reorthogonalize first step
RR(kk,kk) = RR(kk,kk) * RR1;

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
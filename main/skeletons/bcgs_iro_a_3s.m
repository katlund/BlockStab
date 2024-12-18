function [QQ, RR] = bcgs_iro_a_3s(XX, s, musc, param)
% [QQ, RR] = BCGS_IRO_A_3S(XX, s, musc, param) performs a reformulated
% version of BCGS_IRO_A with 3 sync points that skips the first IntraOrtho
% (IO_1).  For more information on how musc can be configured as a multiIO,
% see BCGS_IRO_A.
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

% Set up IO_A, IO_1, and IO_2
if ischar(musc)
    % Defaults
    IO_A = @(W) qr(W,0);
    IO_2 = @(W) IntraOrtho(W, musc, param);
elseif isstruct(musc)
    [musc, musc_param] = unpack_multi_io(musc, param);
    IO_A = @(W) IntraOrtho(W, musc{1}, musc_param{1});
    IO_2 = @(W) IntraOrtho(W, musc{3}, musc_param{3});
end

% Pre-allocate memory for QQ and RR
[m, n] = size(XX);
RR = zeros(n,n);
QQ = zeros(m,n);
p = n/s;

% Set up block indices
kk = 1:s;
sk = s;

% Extract W
W = XX(:,kk);

% IO_A
[QQ(:,kk), RR(kk,kk)] = IO_A(W);

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
    
    % First BCGS step w/o normalization
    S_col = InnerProd(QQ(:,1:sk), XX(:,kk), musc);
    W = XX(:,kk) - QQ(:,1:sk) * S_col;
    
    % Second BCGS step
    Y_col = InnerProd(QQ(:,1:sk), W, musc);
    [QQ(:,kk), Y_diag] = IO_2(W - QQ(:,1:sk) * Y_col);
    
    % Combine both steps
    RR(1:sk,kk) = S_col + Y_col;
    RR(kk,kk) = Y_diag;
    
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
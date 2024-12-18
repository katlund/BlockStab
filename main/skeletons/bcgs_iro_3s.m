function [QQ, RR] = bcgs_iro_3s(XX, s, musc, param)
% [QQ, RR] = BCGS_IRO_3s(XX, s, musc, param) performs Block Classical
% Gram-Schmidt with Inner ReOrthonormalization on the m x n matrix XX with
% p = n/s block partitions each of size s as described in [Barlow &
% Smoktunowicz 2013] but skips the first normalization step in order to
% reduce the total sync count per iteration to 3.
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
    
    % First BCGS step w/o normalization
    S_col = InnerProd(QQ(:,1:sk), XX(:,kk), musc);
    W = XX(:,kk) - QQ(:,1:sk) * S_col;
    
    % Second BCGS step
    Y_col = InnerProd(QQ(:,1:sk), W, musc);
    [QQ(:,kk), Y_diag] = IntraOrtho(W - QQ(:,1:sk) * Y_col, musc, param);
    
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
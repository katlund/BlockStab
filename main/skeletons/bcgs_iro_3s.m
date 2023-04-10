function [QQ, RR] = bcgs_iro_3s(XX, s, musc, verbose)
% [QQ, RR] = BCGS_IRO_3s(XX, s, musc, verbose) performs Block Classical
% Gram-Schmidt on the m x n matrix XX with p = n/s block partitions each of
% size s with Inner ReOrthonormalization as described in [Barlow &
% Smoktunowicz 2013], but skips the first normalization step in order to
% reduce the total sync count per iteration to 3.
%
% See BGS for more details about the parameters, and INTRAORTHO for musc
% options.

%%
addpath(genpath('../'))

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
        norm( eye(s) - QQ(:, 1:s)' * QQ(:, 1:s) ) );
    fprintf('  %2.4e\n',...
        norm( XX(:,1:s) - QQ(:,1:s) * RR(1:s,1:s) ) / norm(XX(:,1:s)) );
end

for k = 1:p-1
    % Update block indices
    kk = kk + s;
    
    % Initialize W
    W = XX(:,kk);
    
    % First BCGS step
    S_col = InnerProd(QQ(:,1:sk), W, musc);
    W = W - QQ(:,1:sk) * S_col;
    
    % Second BCGS step
    T_col = InnerProd(QQ(:,1:sk), W, musc);
    W = W - QQ(:,1:sk) * T_col;
    [QQ(:,kk), T_diag] = IntraOrtho(W, musc);
    
    % Combine both steps
    RR(1:sk,kk) = S_col + T_col;
    RR(kk,kk) = T_diag;
    
    sk = sk + s;
    if verbose
        fprintf('%3.0d:', k+1);
        fprintf('  %2.4e  |',...
            norm( eye(sk) - QQ(:, 1:sk)' * QQ(:, 1:sk) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk) - QQ(:,1:sk) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
    end
end
end
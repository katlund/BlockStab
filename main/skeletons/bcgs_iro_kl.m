function [QQ, RR] = bcgs_pip_iro(XX, s, musc, verbose)
% [QQ, RR] = BCGS_PIP_IRO(XX, s, musc, verbose) performs Block Classical
% Gram-Schmidt with Pythagorean Inner Product modification and Inner
% ReOrthogonalization on the m x n matrix XX with p = n/s block partitions
% each of size s with intra-orthogonalization procedure determined by
% musc.
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

% 1.1
[U, S_diag] = IntraOrtho(XX(:,kk), musc);

% 1.2
[QQ(:,kk), T_diag] = IntraOrtho(U, musc);

% 1.3
RR(kk,kk) = T_diag * S_diag;

if verbose
    fprintf('         LOO      |    RelRes\n');
    fprintf('-----------------------------------\n');
    fprintf('%3.0d:', 1);
    fprintf('  %2.4e  |',...
        norm( eye(s) - QQ(:, 1:s)' * QQ(:, 1:s) ) );
    fprintf('  %2.4e\n',...
        norm( XX(:,1:s) - QQ(:,1:s) * RR(1:s,1:s) ) / norm(XX(:,1:s)) );
end

% Update block indices
kk = kk + s;
sk = sk + s;

% k.1
tmp = [QQ(:,1:sk-s) XX(:,kk)]' * XX(:,kk);
S_col = tmp(1:sk-s,:);
diff = tmp(kk,:) - S_col'*S_col;
S_diag = chol_free(diff);
U = ( XX(:,kk) - QQ(:,1:sk-s) * S_col ) / S_diag;

for k = 2:p
    % k.2
    tmp = [QQ(:,1:sk-s) U]' * U;
    T_col = tmp(1:sk-s,:);
    diff = tmp(kk,:) - T_col' * T_col;
    T_diag = chol_free(diff);
    QQ(:,kk) = ( U - QQ(:,1:sk-s) * T_col ) / T_diag;
    
    % k.3
    RR(1:sk-s,kk) = S_col + T_col * S_diag;
    RR(kk,kk) = T_diag * S_diag;
    
    if verbose
        fprintf('%3.0d:', k);
        fprintf('  %2.4e  |',...
            norm( eye(sk) - QQ(:, 1:sk)' * QQ(:, 1:sk) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk) - QQ(:,1:sk) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
    end

    if k < p
        % Update block indices
        kk = kk + s;
        sk = sk + s;
    
        % (k+1).1
        tmp = [QQ(:,1:sk-s) XX(:,kk)]' * XX(:,kk);
        S_col = tmp(1:sk-s,:);
        diff = tmp(kk,:) - S_col'*S_col;
        S_diag = chol_free(diff);
        U = ( XX(:,kk) - QQ(:,1:sk-s) * S_col ) / S_diag;
    end
end
end
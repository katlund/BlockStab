function [QQ, RR] = bcgs_pio_iro(XX, s, musc, verbose)
% [QQ, RR] = BCGS_PIO_IRO(XX, s, musc, verbose) performs Block Classical Gram-Schmidt
% on the m x n matrix XX with p = n/s block partitions each of size s with
% Inner ReOrthonormalization as described in [Barlow & Smoktunowicz 2013] but with 
% the PIO-variant of BCGS
% Intra-orthonormalization procedure determined by musc.
%
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
    [~, tmp] = IntraOrtho([W zeros(size(W)); zeros(sk-s, s) RR1], musc);
    tmp = tmp' * tmp;
    R1 = chol_nan(tmp(1:s,1:s) - tmp(end-s+1:end, end-s+1:end));
    W = ( W - QQ(:,1:sk-s) * RR1 ) / R1;
    
    % Second step
    RR(1:sk-s,kk) = InnerProd(QQ(:,1:sk-s), W, musc);
    [~, tmp] = IntraOrtho([W zeros(size(W)); zeros(sk-s-s, s) RR(1:sk,kk)], musc);
    tmp = tmp' * tmp;
    RR(kk,kk) = chol_nan(tmp(1:s,1:s) - tmp(end-s+1:end, end-s+1:end));
    QQ(:,kk) = ( W - QQ(:,1:sk-s) * RR(1:sk-s,kk) ) / RR(kk,kk);
    
    % Combine both steps
    RR(1:sk-s,kk) = RR1 + RR(1:sk-s,kk) * R1;
    RR(kk,kk) = RR(kk,kk) * R1;
    
    if verbose
        fprintf('%3.0d:', k);
        fprintf('  %2.4e  |',...
            norm( eye(sk) - InnerProd(QQ(:, 1:sk), QQ(:, 1:sk), musc) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk) - QQ(:,1:sk) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
    end
end
end
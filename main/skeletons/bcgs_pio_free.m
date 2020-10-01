function [QQ, RR] = bcgs_pio_free(XX, s, IOstr, verbose)
% [QQ, RR] = BCGS_PIO_FREE(XX, s, IOstr, verbose) performs Block Classical
% Gram-Schmidt with Pythagorean Intra-Orthogonalization modification on the
% m x n matrix XX with p = n/s block partitions each of size s with
% intra-orthogonalization procedure determined by IOstr.  It uses CHOL_FREE
% instead of Matlab's built-in CHOL.
%
% See BGS for more details about the parameters, and INTRAORTHO for IOstr
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
[QQ(:,kk), RR(kk,kk)] = IntraOrtho(W, IOstr);

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
    
    W = XX(:,kk);     
    S = QQ(:,1:sk)' * W;

    [~, RXS] = IntraOrtho([W zeros(size(W)); zeros(size(S)) S], IOstr);
    RXS = RXS' * RXS;
    diff = RXS(1:s,1:s) - RXS(end-s+1:end, end-s+1:end);
    
    RR(kk,kk) = chol_free(diff); % block version of the Pythagorean theorem
    
    W = W - QQ(:,1:sk) * S;
    
    RR(1:sk,kk) = S;
    QQ(:,kk) = W / RR(kk,kk);
    
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
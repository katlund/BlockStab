function [QQ, RR] = bcgs_pip_iro(XX, s, IOstr, verbose)
% [QQ, RR] = BCGS_PIP_IRO(XX, s, IOstr, verbose) performs Block Classical
% Gram-Schmidt with Pythagorean Inner Product modification and Inner
% ReOrthogonalization on the m x n matrix XX with p = n/s block partitions
% each of size s with intra-orthogonalization procedure determined by
% IOstr.
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
    
    % First step
    W = XX(:,kk);
    
    tmp = [QQ(:,1:sk) W]' * W;
    S1 = tmp(1:sk,:);
    diff = tmp(kk,:) - S1'*S1;    
    
    [~, flag] = chol(diff);
    if ~flag
        R1 = chol(diff); % block version of the Pythagorean theorem
    else
        R1 = NaN;
    end
    
    W = ( W - QQ(:,1:sk) * S1 ) / R1;

    % Second step
    tmp = [QQ(:,1:sk) W]' * W;
    S2 = tmp(1:sk,:);
    diff = tmp(kk,:) - S2'*S2;    
    
    [~, flag] = chol(diff);
    if ~flag
        R2 = chol(diff); % block version of the Pythagorean theorem
    else
        R2 = NaN;
    end
    
    QQ(:,kk) = ( W - QQ(:,1:sk) * S2 ) / R2;
    
    RR(1:sk,kk) = S2 * R1 + S1;
    RR(kk,kk) = R2 * R1;
    
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
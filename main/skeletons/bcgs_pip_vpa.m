function [QQ, RR] = bcgs_pip_vpa(XX, s, IOstr, verbose)
% [QQ, RR] = BCGS_PIP(XX, s, IOstr, verbose) performs Block Classical
% Gram-Schmidt with Pythagorean Inner Product modification on the m x n
% matrix XX with p = n/s block partitions each of size s with
% intra-orthogonalization procedure determined by IOstr. 
% The computation of the input to the calls to Cholesky and the Cholesky 
% factorization itself are both performed in simulated quadruple precision
% Using MATLAB vpa
% BCGS_PIP is a block generalization of CGS-P/Algorithm 2 from 
% [Smoktunowicz et. al. 2006].

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
    
    S = vpa([QQ(:,1:sk) W]',32) * vpa(W,32); % quad precision
    Y = vpa(S(1:sk,:),32); % quad precision  
    diff = vpa(S(kk,:),32)-vpa(Y',32)*vpa(Y,32); % quad precision
    RR(kk,kk) = double(chol_free_vpa(diff)); % block version of the Pythagorean theorem; quad precision

    
    W = W - QQ(:,1:sk) * Y;
    
    RR(1:sk,kk) = Y;
    QQ(:,kk) = W / RR(kk,kk);
    
    sk = sk + s;
    if verbose
        fprintf('%3.0d:', k+1);
        fprintf('  %2.4e  |',...
            norm( eye(sk) - QQ(:, 1:sk)' * QQ(:, 1:sk) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk) - QQ(:,1:sk) * double(RR(1:sk,1:sk)) ) / norm(XX(:,1:sk)) );
    end
end
end
function [QQ, RR] = bcgs_pio_vpa(XX, s, IOstr, verbose)
% [QQ, RR] = BCGS_PIO_VPA(XX, s, IOstr, verbose) performs Block Classical
% Gram-Schmidt with Pythagorean Intra-Orthogonalization modification on the
% m x n matrix XX with p = n/s block partitions each of size s with
% intra-orthogonalization procedure determined by IOstr.  The computation 
% of the input to the calls to Cholesky and the Cholesky factorization itself 
% are both performed in simulated quadruple precision using MATLAB vpa. 
% BCGS_PIO is a
% block generalization of CGS-P/Algorithm 2 from [Smoktunowicz et. al.
% 2006].
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
    RXS = vpa(RXS',32) * vpa(RXS,32); % quad precision
    diff = vpa(RXS(1:s,1:s),32)-vpa(RXS(end-s+1:end, end-s+1:end),32); % quad precision
    RRkkkk = chol_free_vpa(diff); % block version of the Pythagorean theorem; quad precision
    RR(kk,kk) = double(RRkkkk);
    
    W = W - QQ(:,1:sk)*S; 

    
    RR(1:sk,kk) = S;
    
    QQ(:,kk) = double(vpa(W,32) / vpa(RRkkkk,32)); % quad precision
    
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
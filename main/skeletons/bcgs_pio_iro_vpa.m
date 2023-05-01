function [QQ, RR] = bcgs_pio_iro_vpa(XX, s, IOstr, verbose)
% [QQ, RR] = BCGS_PIO_IRO_VPA(XX, s, IOstr, verbose) performs Block Classical Gram-Schmidt
% on the m x n matrix XX with p = n/s block partitions each of size s with
% Inner ReOrthonormalization as described in [Barlow & Smoktunowicz 2013] but with 
% the PIO-variant of BCGS
% Intra-orthonormalization procedure determined by IOstr.
%
% The computation of the input to the calls to Cholesky and the Cholesky 
% factorization itself are both performed in simulated quadruple precision
% using MATLAB vpa. 
%
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
    
    % First BCGSPIO step
    S = QQ(:,1:sk)' * W;

    [~, RXS] = IntraOrtho([W zeros(size(W)); zeros(size(S)) S], IOstr); 
    RXS = vpa(RXS',32) * vpa(RXS,32); % quad precision
    diff = RXS(1:s,1:s) - RXS(end-s+1:end, end-s+1:end); % quad precision
    
    R1 = chol_free_vpa(diff); % block version of the Pythagorean theorem; quad precision
   
    
    W = W - QQ(:,1:sk) * S;
    
    RR1 = S;
    QQ(:,kk) = double(vpa(W,32) / vpa(R1,32)); % quad precision

    
    
    % Second BCGSPIO step
    S = QQ(:,1:sk)' * W;

    [~, RXS] = IntraOrtho([W zeros(size(W)); zeros(size(S)) S]), IOstr); 
    RXS = vpa(RXS',32) * vpa(RXS,32); % quad precision
    diff = vpa(RXS(1:s,1:s),32) - vpa(RXS(end-s+1:end, end-s+1:end),32); % quad precision
    
    RRkkkk = chol_free_vpa(diff); % block version of the Pythagorean theorem; quad precision
    RR(kk,kk) = double(RRkkkk); 
    
    W = W - QQ(:,1:sk) * S;
    
    RR(1:sk,kk) = S;
    QQ(:,kk) = double(vpa(W,32) / vpa(RRkkkk,34)); % quad precision

    
    % Combine both steps
    RR(1:sk,kk) = RR1 + RR(1:sk,kk) * R1;
    RR(kk,kk) = RR(kk,kk) * R1;
    
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
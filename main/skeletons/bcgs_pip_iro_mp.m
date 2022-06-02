function [QQ, RR] = bcgs_pip_iro_mp(XX, s, IOstr, verbose)
% [QQ, RR] = BCGS_PIP_IRO_MP(XX, s, IOstr, verbose) performs Block Classical Gram-Schmidt
% on the m x n matrix XX with p = n/s block partitions each of size s with
% Inner ReOrthonormalization as described in [Barlow & Smoktunowicz 2013] but with 
% the PIP-variant of BCGS.
% Intra-orthonormalization procedure determined by IOstr.
%
% The computation of the input to the calls to Cholesky and the Cholesky 
% factorization itself % are both performed in simulated quadruple precision. 
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
    
    % First BCGSPIP step
    S = mp([QQ(:,1:sk) W]',34) * mp(W,34); % quad precision
    Y = mp(S(1:sk,:),34); % quad precision
    diff = mp(S(kk,:),34) - mp(Y',34)*mp(Y,34); % quad precision
    
    R1 = double(chol_free_mp(diff)); % block version of the Pythagorean theorem; quad precision

    
    W = W - QQ(:,1:sk) * Y;
    
    RR1 = Y;
    QQ(:,kk) = W/R1;
%     S = [QQ(:,1:sk) W]' * W;
%     Y = S(1:sk,:);
%     diff = S(kk,:) - Y'*Y;    
%     
%     R1 = chol_free(diff); % block version of the Pythagorean theorem
%     
%     W = W - QQ(:,1:sk) * Y;
%     
%     RR1 = Y;
%     QQ(:,kk) = W/R1;
    
    
    % Second BCGSPIP step
    S = mp([QQ(:,1:sk) W]',34) * mp(W,34);
    Y = mp(S(1:sk,:),34);
    diff = mp(S(kk,:),34) - mp(Y',34)*mp(Y,34);    
    
    RR(kk,kk) = double(chol_free_mp(diff)); % block version of the Pythagorean theorem; quad precision

    
    W = W - QQ(:,1:sk) * Y;
    
    RR(1:sk,kk) = Y;
    QQ(:,kk) = W/RR(kk,kk);
%     S = [QQ(:,1:sk) W]' * W;
%     Y = S(1:sk,:);
%     diff = S(kk,:) - Y'*Y;    
%     
%     RR(kk,kk) = chol_free(diff); % block version of the Pythagorean theorem
%     
%     W = W - QQ(:,1:sk) * Y;
%     
%     RR(1:sk,kk) = Y;
%     QQ(:,kk) = W/RR(kk,kk);
    
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
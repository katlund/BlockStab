function [QQ, RR] = bcgs_iro_kl(XX, s, ~, verbose)
% [QQ, RR] = BCGS_IRO_KL(XX, s, ~, verbose) performs Block Classical
% Gram-Schmidt on the m x n matrix XX with p = n/s block partitions each of
% size s with Inner ReOrthonormalization (including the first step) and a
% low-sync reformulation derived by KLund.
% 
% See BGS for more details about the parameters, and INTRAORTHO for
% musc options.

%%
addpath(genpath('../'))

% Default: debugging off
if nargin < 4
    verbose = 0;
end

% Pre-allocate memory for QQ and RR and auxiliary matrices
[m, n] = size(XX);
QQ = zeros(m,n);
RR = zeros(n,n);
p = n/s;

% Set up block indices
kk = 1:s;
s1 = 1:s;
s2 = s1 + s;

%% First step + 2.1
% Step 1.1
tmp = XX(:,kk)' * XX(:,kk);
S_diag = chol_free(tmp);
U = XX(:,kk) / S_diag;

% Compute combined inner product (2s x 2s)
tmp = [U, XX(:,kk+s)]' * [U, XX(:,kk+s)];

% Step 1.2
T_diag = chol_free(tmp(s1,s1));
QQ(:,kk) = U / T_diag;

% Step 1.3
RR(kk,kk) = T_diag * S_diag;

% Step 2.1
S_col = (tmp(s2,s1) / T_diag)';
S_diag = chol_free( tmp(s2,s2) - S_col'*S_col );
U = ( XX(:,kk+s) - QQ(:,kk) * S_col ) / S_diag;

if verbose
    fprintf('         LOO      |    RelRes\n');
    fprintf('-----------------------------------\n');
    fprintf('%3.0d:', 1);
    fprintf('  %2.4e  |',...
        norm( eye(s) - QQ(:, 1:s)' * QQ(:, 1:s) ) );
    fprintf('  %2.4e\n',...
        norm( XX(:,1:s) - QQ(:,1:s) * RR(1:s,1:s) ) / norm(XX(:,1:s)) );
end

%% Steps k.2, k.3, and (k+1).1
for k = 2:p-1
    % Update block indices
    kk = kk + s;
    sk = k*s;

    % Compute combined inner product (ks+s x 2s)
    tmp = [QQ(:,1:sk-s), U, XX(:,kk+s)]' * [U, XX(:,kk+s)];
    
    % Step k.2
    T_col = tmp(1:sk-s,s1);
    T_diag = chol_free( tmp(sk-s+1:sk,s1) - T_col'*T_col );
    QQ(:,kk) = U / T_diag;
    
    % Step k.3
    RR(1:sk-s,kk) = S_col + T_col * S_diag;
    RR(kk,kk) = T_diag * S_diag;
    
    % Step (k+1).1
    S_col = [tmp(1:sk-s,s1); 
        T_diag' \ ( tmp(sk-s+1:sk,s2) - T_col'*tmp(1:sk-s,s2) )];
    S_diag = chol_free( tmp(sk+1:sk+s,s2) - S_col'*S_col );
    U = ( XX(:,kk+s) - QQ(:,1:sk) * S_col ) / S_diag;

    if verbose
        fprintf('%3.0d:', k);
        fprintf('  %2.4e  |',...
            norm( eye(sk) - QQ(:, 1:sk)' * QQ(:, 1:sk) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk) - QQ(:,1:sk) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
    end
end

%% Finish up steps p.2 and p.3
% Update block indices
kk = kk + s;
sk = p*s;

% Compute combined inner product (ps x s)
tmp = [QQ(:,1:sk-s), U]' * U;

% Step p.2
T_col = tmp(1:sk-s,s1);
T_diag = chol_free( tmp(sk-s+1:sk,s1) - T_col'*T_col );
QQ(:,kk) = U / T_diag;

% Step k.3
RR(1:sk-s,kk) = S_col + T_col * S_diag;
RR(kk,kk) = T_diag * S_diag;

if verbose
    fprintf('%3.0d:', k);
    fprintf('  %2.4e  |',...
        norm( eye(sk) - QQ(:, 1:sk)' * QQ(:, 1:sk) ) );
    fprintf('  %2.4e\n',...
        norm( XX(:,1:sk) - QQ(:,1:sk) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
end
end
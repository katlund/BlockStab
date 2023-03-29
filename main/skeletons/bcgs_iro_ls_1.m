function [QQ, RR] = bcgs_iro_ls_1(XX, s, musc, verbose)
% [QQ, RR] = BCGS_IRO_LS_1(XX, s, musc, verbose) performs Block Classical
% Gram-Schmidt with Reorthogonalization and a Low-Sync formulation
% (actually, a one-sync) on the n x m matrix XX, along with a tweak for the
% first step (hence the _1).
%
% See BGS for more details about the parameters.

%%
addpath(genpath('../'))

% Default: debugging off
if nargin < 4
    verbose = 0;
end

% Pre-allocate memory for QQ and RR
[m, n] = size(XX);
QQ = zeros(m,n);
RR = zeros(n,n);
p = n/s;

% Set up block indices
kk = 1:s;
s1 = 1:s;
s2 = s+1:2*s;

[Q_tmp, R1] = IntraOrtho(XX(:,kk), musc);

if verbose
    fprintf('         LOO      |    RelRes\n');
    fprintf('-----------------------------------\n');
end

for k = 2:p
    % Increment block index
    kk = kk + s;
    sk = s*k;
    
    % Pull out block vector (keeps MATLAB from copying full X repeatedly)
    Xk = XX(:,kk);
    
    % Compute temporary quantities -- the only sync point!
    if k == 2
        R_tmp = Q_tmp' * [Q_tmp Xk];
    else
        tmp = [QQ(:, 1:sk-2*s) Q_tmp]' * [Q_tmp Xk];
        W = tmp(1:sk-2*s, s1);
        Z = tmp(1:sk-2*s, s2);
        R_tmp = tmp(end-s+1:end,:) - W' * [W Z];
    end
    
    % Pythagorean trick for RR diagonals; assign finished entry
    [~, flag] = chol(R_tmp(:, s1));
    if flag == 0
        R_diag = chol(R_tmp(:, s1));
    else
        R_diag = NaN(s);
    end
    
    if k == 2
    	% Account for help from extra first step
        RR(kk-s, kk-s) = R_diag * R1;
        RR(kk-s, kk) = R_diag' \ R_tmp(:, s2);

        % Finish normalizing QQ(:,k-1)
        QQ(:,kk-s) = Q_tmp / R_diag;
    else
        % Assign finished entries of RR
        RR(kk-s, kk-s) = R_diag;
        RR(kk-s, kk) = R_diag' \ R_tmp(:, s2);
        RR(1:sk-2*s, kk-s) = RR(1:sk-2*s, kk-s) + W;
        RR(1:sk-2*s, kk) = Z;
        
        % Finish normalizing QQ(:,k-1)
        QQ(:,kk-s) = (Q_tmp - QQ(:, 1:sk-2*s) * W) / R_diag;
    end

    % Set up temporary block vector for next iteration
    Q_tmp = Xk - QQ(:, 1:sk-s) * RR(1:sk-s, kk);
    
    if verbose
        fprintf('%3.0d:', k-1);
        fprintf('  %2.4e  |',...
            norm( eye(sk-s) - QQ(:, 1:sk-s)' * QQ(:, 1:sk-s) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk-s) - QQ(:,1:sk-s) * RR(1:sk-s,1:sk-s) ) / norm(XX(:,1:sk-s)) );
    end
end

% Finish renormalizing last basis vector and assign last diagonal entry of
% RR.  Note that this requires just one more sync, no IntraOrtho needed.
tmp = [QQ(:,1:n-s) Q_tmp]' * Q_tmp;
W = tmp(1:n-s,:);
R_tmp = tmp(end-s+1:end,:) - W' * W;
[~, flag] = chol(R_tmp);
if flag == 0
    R_diag = chol(R_tmp);
else
    R_diag = NaN(s);
end
RR(kk, kk) = R_diag;
RR(1:n-s, kk) = RR(1:n-s, kk) + W;
QQ(:,kk) = (Q_tmp - QQ(:,1:n-s) * W) / R_diag;

if verbose
    fprintf('%3.0d:', k+1);
    fprintf('  %2.4e  |', norm( eye(n) - QQ' * QQ ) );
    fprintf('  %2.4e\n', norm( XX - QQ * RR ) / norm(XX) );
end
end

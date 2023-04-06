function [QQ, RR] = bcgs_iro_ls_play(XX, s, ~, verbose)
% [QQ, RR] = BCGS_IRO_LS(XX, s, ~, verbose) performs Block Classical
% Gram-Schmidt with Reorthogonalization and a Low-Sync formulation
% (actually, a one-sync) on the n x m matrix XX. It is the block
% generalization of CGS_IRO_LS, i.e., Algorithm 3 from [Swirydowicz, et.
% al., 2020].  Note that no muscle is explicitly required.
%
% See BGS for more details about the parameters.

% This script has been modified to separate the raw RR entries from the
% final ones, which we hope will simplify the backward error analysis.

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

U = XX(:,kk);

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
        R_tmp = U' * [U Xk];

        % Pythagorean trick for RR diagonals; assign finished entry
        [~, flag] = chol(R_tmp(:, s1));
        if flag == 0
            R_diag = chol(R_tmp(:, s1));
        else
            R_diag = NaN(s);
        end

        % Assign finished entries of RR
        RR(kk-s, kk-s) = R_diag;

        % Finish normalizing QQ(:,k-1)
        QQ(:,kk-s) = U / R_diag;

        % Update S
        S = R_diag' \ R_tmp(:, s2);
        
        % Set up temporary block vector for next iteration
        U = Xk - QQ(:, 1:sk-s) * S;

    elseif k >= 3
        tmp = [QQ(:, 1:sk-2*s) U]' * [U Xk];
        W = tmp(1:sk-2*s, s1);
        Z = tmp(1:sk-2*s, s2);
        R_tmp = tmp(end-s+1:end,:) - W' * [W Z];

        % Pythagorean trick for RR diagonals
        [~, flag] = chol(R_tmp(:, s1));
        if flag == 0
            R_diag = chol(R_tmp(:, s1));
        else
            R_diag = NaN(s);
        end

        % Assign finished entries of RR
        RR(1:sk-2*s, kk-s) = S + W;
        RR(kk-s, kk-s) = R_diag;
        
        % Finish normalizing QQ(:,k-1)
        QQ(:,kk-s) = (U - QQ(:, 1:sk-2*s) * W) / R_diag;

        % Update S (note that it increases by an s x s block)
        S = [Z; R_diag' \ R_tmp(:, s2)];
        
        % Set up temporary block vector for next iteration
        U = Xk - QQ(:, 1:sk-s) * S;
    end
    
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
tmp = [QQ(:,1:n-s) U]' * U;
W = tmp(1:n-s,:);
R_tmp = tmp(end-s+1:end,:) - W' * W;
[~, flag] = chol(R_tmp);
if flag == 0
    R_diag = chol(R_tmp);
else
    R_diag = NaN(s);
end

% Assign finished RR entries
RR(kk, kk) = R_diag;
RR(1:n-s, kk) = S + W;

% Final basis vector
QQ(:,kk) = (U - QQ(:,1:n-s) * W) / R_diag;

if verbose
    fprintf('%3.0d:', k+1);
    fprintf('  %2.4e  |', norm( eye(n) - QQ' * QQ ) );
    fprintf('  %2.4e\n', norm( XX - QQ * RR ) / norm(XX) );
end
end

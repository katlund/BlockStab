function [QQ, RR] = bcgs_iro_ls_vpa(XX, s, ~, verbose)
% [QQ, RR] = BCGS_IRO_LS_vpa(XX, s, ~, verbose) performs Block Classical
% Gram-Schmidt with Reorthogonalization and a Low-Sync formulation
% (actually, a one-sync) on the n x m matrix XX. It is the block
% generalization of CGS_IRO_LS, i.e., Algorithm 3 from [Swirydowicz, et.
% al., 2020].  Note that no muscle is explicitly required.
% Selectively uses quad precision via MATLAB vpa. 
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

Q_tmp = XX(:,kk);

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
        R_tmp = vpa(Q_tmp',32) * vpa([Q_tmp Xk],32); % quad precision
    else
        tmp = vpa([QQ(:, 1:sk-2*s) Q_tmp]',32) * vpa([Q_tmp Xk],32); % quad precision
        W = tmp(1:sk-2*s, s1);
        Z = tmp(1:sk-2*s, s2);
        R_tmp = tmp(end-s+1:end,:) - W' * [W Z]; % quad precision
    end
    
    % Pythagorean trick for RR diagonals; assign finished entry
    R_diag = chol_free_vpa(R_tmp(:, s1)); % quad precision

    
    % Assign finished entries of RR
    RR(kk-s, kk-s) = double(R_diag);
    RR(kk-s, kk) = double(R_diag'\ R_tmp(:, s2));
    
    if k == 2
        % Finish normalizing QQ(:,k-1)
        QQ(:,kk-s) = double(vpa(Q_tmp,32) / vpa(R_diag,32)); % quad precision
    else
        % Assign finished entries of RR
        RR(1:sk-2*s, kk-s) = RR(1:sk-2*s, kk-s) + double(W);
        RR(1:sk-2*s, kk) = double(Z);
        
        % Finish normalizing QQ(:,k-1)
        QQ(:,kk-s) = double(vpa((Q_tmp - QQ(:, 1:sk-2*s) * W),32) / vpa(R_diag,32)); % quad precision
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
tmp = vpa([QQ(:,1:n-s) Q_tmp]',32) * vpa(Q_tmp,32); % quad precision
W = tmp(1:n-s,:);
R_tmp = tmp(end-s+1:end,:) - W' * W; % quad precision
R_diag = chol_free_vpa(R_tmp); % quad precision
RR(kk, kk) = double(R_diag);
RR(1:n-s, kk) = RR(1:n-s, kk) + double(W);
QQ(:,kk) = double(vpa((Q_tmp - QQ(:,1:n-s) * W),32) / vpa(R_diag,32)); % quad precision

if verbose
    fprintf('%3.0d:', k+1);
    fprintf('  %2.4e  |', norm( eye(n) - QQ' * QQ ) );
    fprintf('  %2.4e\n', norm( XX - QQ * RR ) / norm(XX) );
end
end

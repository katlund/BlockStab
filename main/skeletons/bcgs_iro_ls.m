function [QQ, RR] = bcgs_iro_ls(XX, s, musc, param)
% [QQ, RR] = BCGS_IRO_LS(XX, s, musc, param) performs Block Classical
% Gram-Schmidt with Reorthogonalization and a Low-Sync formulation
% (actually, a one-sync) on the n x m matrix XX. It is a block
% generalization of CGS_IRO_LS, i.e., Algorithm 3 from [Swirydowicz, et.
% al. 2020] or Algorithm 2 from [Bielich, et al. 2022].  Note that no
% muscle is explicitly required, because CholQR is hard-coded for all
% intra-orthogonalizations; it can, however, be passed to InnerProd.
%
% This version first appeared in [Carson, et al. 2022].  For an alternative
% 1-sync variation derived from different principles, see BCGS_IRO_1S.
%
% See BGS for more details about the parameters.

%%
% Default: debugging off
if nargin < 4
    param.verbose = 0;
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

% Initialize
U = XX(:,kk);

if param.verbose
    fprintf('         LOO      |    RelRes\n');
    fprintf('-----------------------------------\n');
end

for k = 2:p
    % Increment block index
    kk = kk + s;
    sk = s*k;
    
    % Pull out block vector (keeps MATLAB from copying full X repeatedly)
    W = XX(:,kk);
    
    % Compute temporary quantities -- the only sync point!
    if k == 2
        R_tmp = InnerProd(U, [U W], musc);
    else
        tmp = InnerProd([QQ(:, 1:sk-2*s) U], [U W], musc);
        Y = tmp(1:sk-2*s, s1);
        Z = tmp(1:sk-2*s, s2);
        R_tmp = tmp(end-s+1:end,:) - Y' * [Y Z];
    end
    
    % Pythagorean trick for RR diagonals; assign finished entry
    R_diag = chol_switch(R_tmp(:, s1), param);
    
    % Assign finished entries of RR
    RR(kk-s, kk-s) = R_diag;
    RR(kk-s, kk) = R_diag' \ R_tmp(:, s2);
    
    if k == 2
        % Finish normalizing QQ(:,k-1)
        QQ(:,kk-s) = U / R_diag;
    else
        % Assign finished entries of RR
        RR(1:sk-2*s, kk-s) = RR(1:sk-2*s, kk-s) + Y;
        RR(1:sk-2*s, kk) = Z;
        
        % Finish normalizing QQ(:,k-1)
        QQ(:,kk-s) = (U - QQ(:, 1:sk-2*s) * Y) / R_diag;
    end
    
    % Set up temporary block vector for next iteration
    U = W - QQ(:, 1:sk-s) * RR(1:sk-s, kk);
    
    if param.verbose
        fprintf('%3.0d:', k-1);
        fprintf('  %2.4e  |',...
            norm( eye(sk-s) - InnerProd(QQ(:, 1:sk-s), QQ(:, 1:sk-s), musc) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk-s) - QQ(:,1:sk-s) * RR(1:sk-s,1:sk-s) ) / norm(XX(:,1:sk-s)) );
    end
end

% Finish renormalizing last basis vector and assign last diagonal entry of
% RR.  Note that this requires just one more sync, no IntraOrtho needed.
tmp = InnerProd([QQ(:,1:n-s) U], U, musc);
Y = tmp(1:n-s,:);
R_diag = chol_switch(tmp(end-s+1:end,:) - Y' * Y, param);
RR(kk, kk) = R_diag;
RR(1:n-s, kk) = RR(1:n-s, kk) + Y;
QQ(:,kk) = (U - QQ(:,1:n-s) * Y) / R_diag;

if param.verbose
    fprintf('%3.0d:', k);
    fprintf('  %2.4e  |', norm( eye(n) - InnerProd(QQ, QQ, musc) ) );
    fprintf('  %2.4e\n', norm( XX - QQ * RR ) / norm(XX) );
end
end

function [QQ, RR, TT] = bmgs_icwy(XX, s, IOstr, verbose)
% [QQ, RR, TT] = BMGS_ICWY(XX, s, IOstr, verbose) performs Block Modified
% Gram-Schmidt with Inverse Compact WY form on the m x n matrix XX with p =
% n/s block partitions each of size s with inner orthogonalization
% procedure determined by IOstr. BMGS_ICWY is the block generalization of
% MGS_ICWY.
%
% See BGS for more details about the parameters, and INTRAORTHO for IOstr
% options.

%%
addpath(genpath('../'))

% Default: debugging off
if nargin < 4
    verbose = 0;
end

% Pre-allocate memory for QQ, RR, and TT
[m, n] = size(XX);
QQ = zeros(m,n);
RR = zeros(n,n);
TT = eye(n,n);
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

for k = 1:p-1
    % Update block index
    sk = s*k;
    
    % Pull out vectors (keeps MATLAB from copying full X repeatedly)
    Xk = XX(:,kk+s);
    
    % Compute temporary quantities-- only sync point!
    if k == 1
        tmp = Q_tmp' * [Q_tmp Xk];
        
        [~, flag] = chol(tmp(kk, s1));
        if flag == 0
            R_diag = chol(tmp(kk, s1));
        else
            R_diag = NaN(s);
        end
        
    else
        % I think the problem with LOO stability is here: we are not
        % updating TT(kk,kk), but rather assuming it is identity.  Perhaps
        % PIO version of Pythagorean trick will work?
        
        tmp = [QQ(:,1:sk-s) Q_tmp]' * [Q_tmp Xk];

        
        [~, flag] = chol(tmp(kk, s1));
        if flag == 0
            R_diag = chol(tmp(kk, s1));
        else
            R_diag = NaN(s);
        end
        
        TT(1:sk-s,kk) = tmp(1:sk-s,s1) / R_diag;
    end
    
    RR(kk,kk) = R_diag;
    tmp(kk,s2) = tmp(kk,s2) / R_diag;
    RR(1:sk,kk+s) = TT(1:sk,1:sk)' \ tmp(:,s2);
    
    QQ(:,kk) = Q_tmp / R_diag;
    Q_tmp = Xk - QQ(:,1:sk) * RR(1:sk,kk+s);
    
    % Update block index
    kk = kk + s;
    
    if verbose
        fprintf('%3.0d:', k);
        fprintf('  %2.4e  |',...
            norm( eye(sk) - QQ(:, 1:sk)' * QQ(:, 1:sk) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk) - QQ(:,1:sk) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
    end
end
[QQ(:,kk), RR(kk,kk)] = IntraOrtho(Q_tmp, IOstr);

if verbose
    fprintf('%3.0d:', k+1);
    fprintf('  %2.4e  |', norm( eye(n) - QQ' * QQ ) );
    fprintf('  %2.4e\n', norm( XX - QQ * RR ) / norm(XX) );
end
end

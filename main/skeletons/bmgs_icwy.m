function [QQ, RR, TT] = bmgs_icwy(XX, s, musc, param)
% [QQ, RR, TT] = BMGS_ICWY(XX, s, musc, param) performs Block Modified
% Gram-Schmidt with Inverse Compact WY form on the m x n matrix XX with p =
% n/s block partitions each of size s with inner orthogonalization
% procedure determined by musc. BMGS_ICWY is the block generalization of
% MGS_ICWY.
%
% See BGS for more details about the parameters, and INTRAORTHO for musc
% options.
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
addpath(genpath('../'))

% Default: debugging off
if nargin < 4
    param.verbose = 0;
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

if param.verbose
    fprintf('         LOO      |    RelRes\n');
    fprintf('-----------------------------------\n');
end

for k = 1:p-1
    % Update block index
    sk = s*k;
    
    % Pull out vectors (keeps MATLAB from copying full X repeatedly)
    W = XX(:,kk+s);
    
    % Compute temporary quantities-- only sync point!
    if k == 1
        tmp = InnerProd(Q_tmp, [Q_tmp W], musc);
        R_diag = chol_switch(tmp(kk, s1), param);
        
    else
        tmp = InnerProd([QQ(:,1:sk-s) Q_tmp], [Q_tmp W], musc);
        R_diag = chol_switch(tmp(kk, s1), param);        
        TT(1:sk-s,kk) = tmp(1:sk-s,s1) / R_diag;
    end
    
    RR(kk,kk) = R_diag;
    tmp(kk,s2) =  R_diag' \ tmp(kk,s2);  % tmp(kk,s2) = Q_tmp, W and R_diag is block norm of Q_tmp, so we pre-multiply by the inverse of the transpose
    RR(1:sk,kk+s) = TT(1:sk,1:sk)' \ tmp(:,s2);
    
    QQ(:,kk) = Q_tmp / R_diag;
    Q_tmp = W - QQ(:,1:sk) * RR(1:sk,kk+s);
    
    % Update block index
    kk = kk + s;
    
    if param.verbose
        fprintf('%3.0d:', k);
        fprintf('  %2.4e  |',...
            norm( eye(sk) - InnerProd(QQ(:, 1:sk), QQ(:, 1:sk), musc) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk) - QQ(:,1:sk) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
    end
end
[QQ(:,kk), RR(kk,kk)] = IntraOrtho(Q_tmp, musc);

if param.verbose
    fprintf('%3.0d:', k+1);
    fprintf('  %2.4e  |', norm( eye(n) - InnerProd(QQ, QQ, musc) ) );
    fprintf('  %2.4e\n', norm( XX - QQ * RR ) / norm(XX) );
end
end

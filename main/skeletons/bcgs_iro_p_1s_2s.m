function [QQ, RR] = bcgs_iro_p_1s_2s(XX, s, musc, param)
% [QQ, RR] = BCGS_IRO_P_1S_2S(XX, s, musc, param) performs Block Classical
% Gram--Schmidt on the m x n matrix XX with p = n/s block partitions each
% of size s with intra-orthogonalization procedure determined by musc.
% This subroutine combines BCGS_IRO_P_1S and BCGS_IRO_P_2S to reduce sync 
% points as much as possible.
% HouseQR is fixed for the first vector (_a).
%
% See Algorithm 4 from [Carson & Ma, 2024] for details.
%
% See BGS for more details about the parameters, and INTRAORTHO for musc
% options.
%
% Part of [BlockStab](https://github.com/katlund) package.  Check README
% for how to properly cite and reuse this file.


[QQ, RR, flag, sk] = bcgs_iro_p_1s_sub(XX, s, musc, param);

if flag ~= 0
    [QQ, RR] = bcgs_iro_p_2s_sub(XX, s, musc, param, QQ, RR, sk);
end

end

function [QQ, RR, flag, sk] = bcgs_iro_p_1s_sub(XX, s, musc, param)
% [QQ, RR] = BCGS_IRO_P_1S_sub(XX, s, musc, param) performs Block Classical
% Gram--Schmidt with one sync on the m x n matrix XX with p = n/s block
% partitions each of size s.
% This subroutine is similar to BCGS_IRO_P_1S, but employs additional
% condition to determine when to switch to BCGS_IRO_P_2S iterations
% adaptively.

%%
% Default: debugging off
if nargin < 4
    param.verbose = 0;
end

% Pre-allocate memory for QQ and RR
[m, n] = size(XX);
RR = zeros(n,n);
QQ = zeros(m,n);
p = n/s;
flag = 0;

% Set up block indices
kk = 1:s;
sk = s;
s1 = kk;
s2 = s1 + s;

% Strong first step
W = XX(:,kk);
normX = norm(W, 'fro');
[QQ(:,kk), RR(kk,kk)] = qr(W, 0);

if param.verbose
    fprintf('         LOO      |    RelRes\n');
    fprintf('-----------------------------------\n');
    fprintf('%3.0d:', 1);
    fprintf('  %2.4e  |',...
        norm( eye(s) - InnerProd(QQ(:, 1:s), QQ(:, 1:s), musc) ) );
    fprintf('  %2.4e\n',...
        norm( XX(:,1:s) - QQ(:,1:s) * RR(1:s,1:s) ) / norm(XX(:,1:s)) );
end

% k.1.1
tmp = InnerProd([QQ(:,1:sk), XX(:,kk+s)], XX(:,kk+s), musc);
S_col = tmp(1:sk,:);
omega = tmp(kk+s,:);
[Skk, nanflag] = chol_switch(omega - S_col'*S_col, param);
if nanflag ~= 0
    flag = 4;
    sk = sk-s;
    return
end

% Update block indices
kk = kk + s;
sk = sk + s;

% k.1.2
W = XX(:,kk) - QQ(:,1:sk-s) * S_col;
U = W / Skk;

for k = 2:p-1
    % k.2
    tmp = InnerProd([QQ(:,1:sk-s) U XX(:,kk+s)], [U XX(:,kk+s)], musc);
    normX = sqrt(normX^2 + norm(XX(:,kk+s), 'fro')^2);
    Y_col = tmp(1:sk-s,s1);
    
    eigvU2 = eig((tmp(kk,s1)+tmp(kk,s1)')/2.0);
    if 3*min(abs(eigvU2)) <= max(abs(eigvU2)) 
        flag = 1;
        sk = sk - s;
        return
    end
    [Y_diag, nanflag] = chol_switch(tmp(kk,s1) - Y_col' * Y_col, param);
    if nanflag ~= 0
        flag = 3;
        sk = sk-s;
        return
    end
    QQ(:,kk) = (U - QQ(:,1:sk-s) * Y_col ) / Y_diag;
    
    % k.3
    RR(1:sk-s,kk) = S_col + Y_col*Skk;
    RR(kk,kk) = Y_diag*Skk;
    
    if param.verbose
        fprintf('%3.0d:', k);
        fprintf('  %2.4e  |',...
            norm( eye(sk) - InnerProd(QQ(:, 1:sk), QQ(:, 1:sk), musc) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk) - QQ(:,1:sk) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
    end

    % (k+1).1.1
    Z_col = tmp(1:sk-s,s2);
    P = tmp(kk,s2);
    S_col = [Z_col; Y_diag' \ (P - Y_col' * Z_col)];
    

    % Update block indices
    kk = kk + s;
    sk = sk + s;

    % (k+1).1.2
    [Skk, nanflag] = chol_switch(tmp(kk,s2) - S_col'*S_col, param);
    if nanflag ~= 0
        flag = 2;
        sk = sk-s;
        return
    end

    W = XX(:,kk) - QQ(:,1:sk-s) * S_col;
    U = W/Skk;
    
end

% p.2
tmp = InnerProd([QQ(:,1:sk-s) U], U, musc);
Y_col = tmp(1:sk-s,s1);
Y_diag = chol_switch(tmp(kk,s1) - Y_col' * Y_col, param);
QQ(:,kk) = (U - QQ(:,1:sk-s) * Y_col ) / Y_diag;

% p.3
RR(1:sk-s,kk) = S_col + Y_col*Skk;
RR(kk,kk) = Y_diag*Skk;

if param.verbose
    fprintf('%3.0d:', k);
    fprintf('  %2.4e  |',...
        norm( eye(sk) - InnerProd(QQ(:, 1:sk), QQ(:, 1:sk), musc) ) );
    fprintf('  %2.4e\n',...
        norm( XX(:,1:sk) - QQ(:,1:sk) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
end
end


function [QQ, RR] = bcgs_iro_p_2s_sub(XX, s, musc, param, QQ, RR, sk)
% [QQ, RR] = BCGS_IRO_P_2S_sub(XX, s, musc, param) performs Block Classical
% Gram--Schmidt with one sync on the matrix XX(:, sk-s+1:n) with
% p = (n-sk+1)/s block partitions each of size s.
% This subroutine is similar to BCGS_IRO_P_2S, but only orthogonalizes
% XX(:, sk-s+1:n). The first sk-s columns has been handled by
% BCGS_IRO_P_1S_sub.

%%
% Default: debugging off
if nargin < 4
    param.verbose = 0;
end

% Pre-allocate memory for QQ and RR
[m, n] = size(XX);
p = (n-sk+1)/s;

% Set up block indices
kk = sk-s+1:sk;
s1 = 1:s;
s2 = s1 + s;


if param.verbose
    fprintf('         LOO      |    RelRes\n');
    fprintf('-----------------------------------\n');
    fprintf('%3.0d:', 1);
    fprintf('  %2.4e  |',...
        norm( eye(s) - InnerProd(QQ(:, 1:s), QQ(:, 1:s), musc) ) );
    fprintf('  %2.4e\n',...
        norm( XX(:,1:s) - QQ(:,1:s) * RR(1:s,1:s) ) / norm(XX(:,1:s)) );
end

% k.1.1
S_col = InnerProd(QQ(:,1:sk), XX(:,kk+s), musc);

% Update block indices
kk = kk + s;
sk = sk + s;

% % k.1.2
W = XX(:,kk) - QQ(:,1:sk-s) * S_col;
[U, Skk] = IntraOrtho(W, musc, param);

for k = 2:p
    % k.2
    tmp = InnerProd([QQ(:,1:sk-s) U], [U XX(:,kk+s)], musc);
    Y_col = tmp(1:sk-s,s1);
    Y_diag = chol_switch(tmp(kk,s1) - Y_col' * Y_col, param);
    QQ(:,kk) = (U - QQ(:,1:sk-s) * Y_col ) / Y_diag;
    
    % k.3
    RR(1:sk-s,kk) = S_col + Y_col*Skk;
    RR(kk,kk) = Y_diag*Skk;
    
    if param.verbose
        fprintf('%3.0d:', k);
        fprintf('  %2.4e  |',...
            norm( eye(sk) - InnerProd(QQ(:, 1:sk), QQ(:, 1:sk), musc) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk) - QQ(:,1:sk) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
    end

    % (k+1).1.1
    Z_col = tmp(1:sk-s,s2);
    P = tmp(kk,s2);
    S_col = [Z_col; Y_diag' \ (P - Y_col' * Z_col)];
    

    % Update block indices
    kk = kk + s;
    sk = sk + s;

    % (k+1).1.2
    W = XX(:,kk) - QQ(:,1:sk-s) * S_col;
    [U, Skk] = IntraOrtho(W, musc, param);
end

% p.2
tmp = InnerProd([QQ(:,1:sk-s) U], U, musc);
Y_col = tmp(1:sk-s,s1);
Y_diag = chol_switch(tmp(kk,s1) - Y_col' * Y_col, param);
QQ(:,kk) = (U - QQ(:,1:sk-s) * Y_col ) / Y_diag;

% p.3
RR(1:sk-s,kk) = S_col + Y_col*Skk;
RR(kk,kk) = Y_diag*Skk;

if param.verbose
    fprintf('%3.0d:', k);
    fprintf('  %2.4e  |',...
        norm( eye(sk) - InnerProd(QQ(:, 1:sk), QQ(:, 1:sk), musc) ) );
    fprintf('  %2.4e\n',...
        norm( XX(:,1:sk) - QQ(:,1:sk) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
end
end
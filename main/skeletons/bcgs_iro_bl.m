function [QQ, RR] = bcgs_iro_bl(XX, s, ~, verbose)
% [QQ, RR] = BCGS_IRO_BL(XX, s, ~, verbose) performs Block Classical
% Gram-Schmidt with Reorthogonalization a la Bielich and Langou
% (https://github.com/dbielich/DCGS2) on the n x m matrix XX. It is the
% block generalization of CGS_IRO_BL, i.e., Algorithm 2 from [Bielich
% et. al., 2022].  Note that no muscle is explicitly required.
%
% See BGS for more details about the parameters.

%%
addpath(genpath('../'))

% Default: debugging off
if nargin < 4
    verbose = 0;
end

n = size(XX,2);
p = n/s;
QQ = XX;
RR = zeros(n,n);

if verbose
    fprintf('         LOO      |    RelRes\n');
    fprintf('-----------------------------------\n');
end

kk = 1:s;
for k = 2:p
    sk = s*k;
    kk = kk + s;
    [QQ(:,sk-2*s+1:sk), RR(1:sk-s,sk-2*s+1:sk)] =...
        orth_dcgs2_qr(QQ(:, 1:sk), RR(1:sk-2*s, kk-s), s);

    if verbose
        fprintf('%3.0d:', k-1);
        fprintf('  %2.4e  |',...
            norm( eye(sk-s) - QQ(:,1:sk-s)' * QQ(:,1:sk-s) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk-s) - QQ(:,1:sk-s) * RR(1:sk-s,1:sk-s) ) / norm(XX(:,1:sk-s)) );
    end
end
[QQ(:,kk), RR(:,kk)] = orth_dcgs2_qr_cleanup(QQ, RR(1:end-s,kk), s);
if verbose
    fprintf('%3.0d:', k);
    fprintf('  %2.4e  |', norm( eye(n) - QQ' * QQ ) );
    fprintf('  %2.4e\n', norm( XX - QQ * RR ) / norm(XX) );
end
end

%% Subroutines that add one vector at a time
function [Q, R] = orth_dcgs2_qr(QQ, RR, s)
[m, sk] = size(QQ);
k = sk/s;

R = zeros(sk-s, 2*s);
Q = zeros(m, 2*s);

s1 = 1:s;
s2 = s+1:2*s;
if k == 2
    R_tmp = QQ(:, s1)' * QQ;

    [~, flag] = chol( R_tmp(:, s1) );
    if flag == 0
        R_diag = chol( R_tmp(:, s1) );
    else
        R_diag = NaN(s);
    end

    R(s1, s1) = R_diag;
    R(s1, s2) = R_diag' \ R_tmp(:, s2);

    Q(:,s1) = QQ(:,s1) / R_diag;
    Q(:,s2) = QQ(:,s2) - Q(:,s1) * R(s1,s2); 
end

kk = sk-s+1:sk;
if k >= 3
    tmp = QQ(:, 1:sk-s)' * QQ(:, sk-2*s+1:sk);
    W = tmp(1:sk-2*s, s1);
    Z = tmp(1:sk-2*s, s2);
    R_tmp = tmp(kk-s, :) - W' * [W Z];
    
    [~, flag] = chol( R_tmp(:, s1) );
    if flag == 0
        R_diag = chol( R_tmp(:, s1) );
    else
        R_diag = NaN(s);
    end

    R(kk-s, s1) = R_diag;
    R(kk-s, s2) = R_diag' \ R_tmp(:, s2);

    R(1:sk-2*s, s1) = RR(1:sk-2*s, s1) + W;
    R(1:sk-2*s, s2) = Z;
        
    Q(:, s1) = ( QQ(:, kk-s) - QQ(:, 1:sk-2*s) * W ) / R_diag;
    Q(:, s2) = QQ(:, kk) - [QQ(:, 1:sk-2*s) Q(:, s1)] * R(1:sk-s, s2);
end
end

function [Q, R] = orth_dcgs2_qr_cleanup(QQ, RR, s)
[m, n] = size(QQ);

R = zeros(n,1);
Q = zeros(m,1);

s1 = 1:s;
ps = n;
pp = ps-s+1:ps;

R(:, s1) = QQ' * QQ(:, pp);
tmp = QQ' * QQ(:, pp);
W = tmp(1:ps-s, s1);
R_tmp = tmp(pp, s1) - W' * W;

[~, flag] = chol( R_tmp );
if flag == 0
    R_diag = chol( R_tmp );
else
    R_diag = NaN(s);
end

R(pp, s1) = R_diag;
R(1:ps-s, s1) = RR(1:ps-s, s1) + W;

Q(:,s1) = ( QQ(:, pp) - QQ(:, 1:ps-s) * W ) / R_diag;
end
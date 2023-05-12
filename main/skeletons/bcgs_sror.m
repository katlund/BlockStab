function [QQ, RR] = bcgs_sror(XX, s, rpltol, verbose)
% [QQ, RR] = BCGS_SROR(XX, s, rpltol) performs Block Classical Gram-Schmidt
% with Selective ReOrthogonalization and Replacement on the m x n matrix XX
% with p = n/s block partitions each of size s and with CGS_SROR as the
% fixed IntraOrtho procedure.  BCGS_SROR is the block generalization of
% CGS_SROR as proposed in [Stewart 2008].  rpltol is nonnegative scalar
% that determines when to replace "bad" vectors; the larger the number, the
% more replacement allowed, which can have an adverse effect on the
% residual.  The default is set to 1.
%
% See BGS for more details about the parameters.
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
addpath(genpath('../'))

% Default: debugging off
if nargin < 4
    verbose = 0;
end

% Pre-allocate memory for QQ and RR
[m, n] = size(XX);
p = n/s;
QQ = zeros(m, n);
RR = zeros(p);

% Default for rpltol
if isempty(rpltol)
    rpltol = 1;
end

% Set up block indices
kk = 1:s;
sk = s;

[QQ(:,kk), ~, RR(kk,kk)] = bcgs_step_sror(zeros(m,0), XX(:,kk), rpltol);

if verbose
    fprintf('         LOO      |    RelRes\n');
    fprintf('-----------------------------------\n');
    fprintf('%3.0d:', 1);
    fprintf('  %2.4e  |',...
        norm( eye(s) - InnerProd(QQ(:, 1:s), QQ(:, 1:s)) ) );
    fprintf('  %2.4e\n',...
        norm( XX(:,1:s) - QQ(:,1:s) * RR(1:s,1:s) ) / norm(XX(:,1:s)) );
end

for k = 1:p-1
    % Update block indices
    kk = kk + s;
    
    [QQ(:,kk), RR(1:sk,kk), RR(kk,kk)] =...
        bcgs_step_sror(QQ(:,1:sk), XX(:,kk), rpltol);
    
    sk = sk + s;
    if verbose
        fprintf('%3.0d:', k+1);
        fprintf('  %2.4e  |',...
            norm( eye(sk) - InnerProd(QQ(:, 1:sk), QQ(:, 1:sk)) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk) - QQ(:,1:sk) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
    end
end
end


function [Y, R12, R22] = bcgs_step_sror(QQ, X, rpltol)
% [Y, R12, R22] = BCGS_STEP_SROR(QQ, X, rpltol) returns an orthonormal
% matrix Y whose columns lie in the orthogonal complement of the column
% space of QQ, a matrix R12, and an upper triangular matrix R22 satisfying
% 
% (*)   X = QQ*R12 + Y*R22.
% 
% The optional argument rpltol has a default value of 1.  Increasing it,
% say to 100, may improve performance, but will degrade the accuracy of the
% relation (*).
%
% Originally written by G. W. Stewart, 2008.  Modified by Kathryn Lund,
% 2020.

%%
%  Initialization.
reorth = false;
if nargin < 3
  rpltol = 1;
else
    if rpltol < 1
        rpltol = 1;
    end
end

s = size(X,2);
sk = size(QQ,2);

nu = zeros(1,s);
for k = 1:s
  nu(k) = norm(X(:,k));
end

%  Beginning of the first orthogonalization step.  Project Y
%  onto the orthogonal complement of Q.
R12 = InnerProd(QQ, X);
Y = X - QQ*R12;
R22 = zeros(s);

%  Orthogonalize the columns of Y.
for k = 1:s
    [yk, r, rho] = cgs_step_sror(Y(:,1:k-1), Y(:,k), nu(k), rpltol);
    Y(:,k) = yk;
    R22(1:k-1,k) = r;
    R22(k,k) = rho;
    if rho <= 0.5*nu(k)
        reorth = true;
    end
end

%  Return if Q has zero columns.
if sk == 0 || ~reorth
   return
end

% Beginning of the reorthogonalization.  Project Y onto the orthogonal
% complement of Q.
S12 = InnerProd(QQ, Y);
Y = Y - QQ*S12;

% Orthogonalize the columns of Y.
S22 = zeros(s);
for k = 1:s
    [yk, r, sig] = cgs_step_sror(Y(:,1:k-1), Y(:,k), 1, rpltol);

    if sig < 0.5
    % The result, yk, is not satisfactorily orthogonal.  Perform an
    % unblocked orthogonalization against Q and the previous columns of Y.
        [yk, r, sig] = cgs_step_sror([QQ, Y(:,1:k-1)], Y(:,k), rpltol);
        S12(:,k) = S12(:,k) + r(1:sk);
        r = r(sk+1:sk+k-1);
    end
    Y(:,k) = yk;
    S22(1:k-1,k) = r;
    S22(k,k) = sig;
end

R12 = R12 + S12*R22;
R22 = S22*R22;
end
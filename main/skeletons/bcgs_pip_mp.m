function [QQ, RR] = bcgs_pip_mp(XX, s, musc, param)
% [QQ, RR] = BCGS_PIP_MP(XX, s, musc, param) performs Block Classical
% Gram-Schmidt with Pythagorean Inner Product modification on the m x n
% matrix XX with p = n/s block partitions each of size s with
% intra-orthogonalization procedure determined by musc. BCGS_PIP is a block
% generalization of CGS-P/Algorithm 2 from [Smoktunowicz et. al. 2006]
% derived in [Carson et al. 2021].
%
% This multiprecision version computes the inputs to Cholesky and the
% Cholesky factorization itself in simulated quadruple (or other,
% user-specified precision) precision.  See MP_SWITCH for details on the
% param struct.
%
% See BGS and MP_SWITCH for more details about the parameters, and
% INTRAORTHO for musc options.
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
% Defaults
if nargin < 4
    param = mp_param_init;
elseif nargin == 4
    param = mp_param_init(param);
end

% Set up primary and secondary precision routines (p1 & p2, respectively)
p1 = @(x) mp_switch(x, param.mp_package, param.mp_pair{1});
p2 = @(x) mp_switch(x, param.mp_package, param.mp_pair{2});

% Pre-allocate memory for QQ and RR
[m, n] = size(XX);
RR = p1(zeros(n,n));
QQ = p1(zeros(m,n));
p = n/s;

% Set up block indices
kk = 1:s;
sk = s;

% Pull out block vector (keeps MATLAB from copying full X repeatedly)
W = p1(XX(:,kk)); % cast in case p1 is single

% First step
[QQ(:,kk), RR(kk,kk)] = IntraOrtho(W, musc, param);

if param.verbose
    fprintf('         LOO      |    RelRes\n');
    fprintf('-----------------------------------\n');
    fprintf('%3.0d:', 1);
    fprintf('  %2.4e  |',...
        norm( eye(s) - InnerProd(QQ(:, 1:s), QQ(:, 1:s), musc) ) );
    fprintf('  %2.4e\n',...
        norm( XX(:,1:s) - QQ(:,1:s) * RR(1:s,1:s) ) / norm(XX(:,1:s)) );
end

for k = 2:p
    % Update block indices
    kk = kk + s;
    sk = sk + s;
    
    % Pull out block vector (keeps MATLAB from copying full X repeatedly)
    W = p1(XX(:,kk));  % cast in case p1 is single
    
    % Compute only W^T * W in p2.  Note that standard operations are
    % overloaded.  Note also that although we split InnerProd across two
    % steps to handle different precisions, we still regard these steps as
    % a single sync point.
    RRk = InnerProd(QQ(:,1:sk-s), W);
    tmp = InnerProd(p2(W), p2(W), musc); % returned in p2

    % Compute Cholesky in p2
    R_diag = chol_switch(tmp - p2(RRk)' * p2(RRk), param); % returned in p2
    
    % Assign RR in p1
    RR(1:sk-s,kk) = RRk;
    RR(kk,kk) = p1(R_diag); 

    % Compute next basis vector in p2 and cast to p1
    QQ(:,kk) = p1(p2(W - QQ(:,1:sk-s) * RRk) / R_diag);
    
    if param.verbose
        fprintf('%3.0d:', k-1);
        fprintf('  %2.4e  |',...
            norm( eye(sk-s) - InnerProd(QQ(:, 1:sk-s), QQ(:, 1:sk-s), musc) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk-s) - QQ(:,1:sk-s) * RR(1:sk-s,1:sk-s) ) / norm(XX(:,1:sk-s)) );
    end
end
end
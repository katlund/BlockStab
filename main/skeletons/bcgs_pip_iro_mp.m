function [QQ, RR] = bcgs_pip_iro_mp(XX, s, musc, param)
% [QQ, RR] = BCGS_PIP_IRO_MP(XX, s, musc, param) performs BCGS_PIP with
% Inner ReOrthogonalization on the m x n matrix XX with p = n/s block
% partitions each of size s and with intra-orthogonalization procedure
% determined by musc.
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

% Set up quad-precision subroutine
qp = @(x) mp_switch(x, param);

% Pre-allocate memory for QQ and RR
[m, n] = size(XX);
RR = zeros(n,n);
QQ = zeros(m,n);
p = n/s;

% Set up block indices
kk = 1:s;
sk = s;

% Pull out block vector (keeps MATLAB from copying full X repeatedly)
W = XX(:,kk);

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
    W = XX(:,kk);
    
    %% First sync point
    % Compute temporary quantities (inner product) in qp.  Note that
    % standard operations are overloaded.
    tmp = InnerProd(qp([QQ(:,1:sk-s) W]), qp(W), musc); % returned in qp
    RR1 = tmp(1:sk-s,:); % returned in qp

    % Compute Cholesky in qp
    R1 = chol_switch( tmp(kk,:) - RR1' * RR1, param ); % returned in qp

    % Compute intermediate basis vector and cast to double
    W = double(qp(W - QQ(:,1:sk-s) * RR1) / R1);
    
    %% Second sync point
    % Again compute inner product in qp
    tmp = InnerProd(qp([QQ(:,1:sk-s) W]), qp(W), musc); % returned in qp

    % Compute Cholesky in qp
    R2 = chol_switch( tmp(kk,:) - tmp(1:sk-s,:)' * tmp(1:sk-s,:), param); % returned in qp
    
    % Compute next basis vector and cast to double
    QQ(:,kk) = double(qp(W - QQ(:,1:sk-s) * tmp(1:sk-s,:)) / R2);
    
    % Assign RR in double and combine both steps
    RR(1:sk-s,kk) = double(RR1) + double(tmp(1:sk-s,:)) * double(R1);
    RR(kk,kk) = double(R2) * double(R1);
    
    if param.verbose
        fprintf('%3.0d:', k-1);
        fprintf('  %2.4e  |',...
            norm( eye(sk-s) - InnerProd(QQ(:, 1:sk-s), QQ(:, 1:sk-s), musc) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk-s) - QQ(:,1:sk-s) * RR(1:sk-s,1:sk-s) ) / norm(XX(:,1:sk-s)) );
    end
end
end
function [QQ, RR] = bcgs_iro_a(XX, s, musc, param)
% [QQ, RR] = BCGS_IRO_A(XX, s, musc, param) performs BCGS_IRO with the
% option of musc being a multiIO.
%
% _A methods can accept a multiIO struct for musc.  When musc is provided
% as a single string, then IO_A is HouseQR and IO_1 (i.e., the IntraOrtho
% within the loop) is musc.  Otherwise, musc must be a struct with .io_a
% and .io_1 fields.  In particular,
%       struct('io_a', 'houseqr', 'io_1', 'cholqr')
% reproduces the default for musc = 'cholqr'.
%
% Note that a multiIO struct may have parameters for the musc encoded in
% the struct as a subfield param.
%
% See BGS for more details about the parameters, and INTRAORTHO for musc
% options.
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
% Default: debugging off
if nargin < 4
    param.verbose = 0;
end

% Set up IO_A, IO_1, and IO_2
if ischar(musc)
    % Defaults
    IO_A = @(W) qr(W,0);
    IO_1 = @(W) IntraOrtho(W, musc, param);
    IO_2 = @(W);
elseif isstruct(musc)
    [musc, musc_param] = unpack_multi_io(musc, param);
    IO_A = @(W) IntraOrtho(W, musc{1}, musc_param{1});
    IO_1 = @(W) IntraOrtho(W, musc{2}, musc_param{2});
    IO_2 = @(W) IntraOrtho(W, musc{3}, musc_param{3});
end

% Pre-allocate memory for QQ and RR
[m, n] = size(XX);
RR = zeros(n,n);
QQ = zeros(m,n);
p = n/s;

% Set up block indices
kk = 1:s;
sk = s;

% Extract W
W = XX(:,kk);

% IO_A
[QQ(:,kk), RR(kk,kk)] = IO_A(W);

if param.verbose
    fprintf('         LOO      |    RelRes\n');
    fprintf('-----------------------------------\n');
    fprintf('%3.0d:', 1);
    fprintf('  %2.4e  |',...
        norm( eye(s) - InnerProd(QQ(:, 1:s), QQ(:, 1:s), musc{1}) ) );
    fprintf('  %2.4e\n',...
        norm( XX(:,1:s) - QQ(:,1:s) * RR(1:s,1:s) ) / norm(XX(:,1:s)) );
end

for k = 1:p-1
    % Update block indices
    kk = kk + s;

    W = XX(:,kk);
    
    % First BCGS step
    RR1 = InnerProd(QQ(:,1:sk), W, musc{2});
    W = W - QQ(:,1:sk) * RR1;
    [W, R1] = IO_1(W);
    
    % Second BCGS step
    RR(1:sk,kk) = InnerProd(QQ(:,1:sk), W, musc{3});
    W = W - QQ(:,1:sk) * RR(1:sk,kk);
    [QQ(:,kk), RR(kk,kk)] = IO_2(W);
    
    % Combine both steps
    RR(1:sk,kk) = RR1 + RR(1:sk,kk) * R1;
    RR(kk,kk) = RR(kk,kk) * R1;
    
    sk = sk + s;
    if param.verbose
        fprintf('%3.0d:', k+1);
        fprintf('  %2.4e  |',...
            norm( eye(sk) - InnerProd(QQ(:, 1:sk), QQ(:, 1:sk), musc{3}) ) );
        fprintf('  %2.4e\n',...
            norm( XX(:,1:sk) - QQ(:,1:sk) * RR(1:sk,1:sk) ) / norm(XX(:,1:sk)) );
    end
end
end
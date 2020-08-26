 function RunTestSingValLowSync(XXdim, singval, musc)
% RUNTESTSINGVALLOWSYNC(XXdim, singval, musc) is a wrapper function that
% compares stability for different muscles for a set of matrices of size
% XXdim = [m s] with varying singular values specified by the vector array
% singval.
%
% musc should be given as either a char array or a cell of char arrays
% (i.e., text strings with single quotes).
%
% Options for XXdim:
%   XXdim = [m s], where m is the number of rows and s is the number of
%   columns.
%
%   Default: XXdim = [1000 200]
%
% Options for singval:
%   singval should be a vector of negative powers.  Each entry of singval
%   is treated as the power of the smallest singular value for a matrix
%   whose singular values range from 10^0 to 10^t, where t is an entry of
%   singval.  For each t, a different matrix is generated, in order to see
%   how BGS behaves for matrices with different condition numbers.
%   
%   Default: singval = -(1:16)
%
% Options for musc: see INTRAORTHO
%
%   Default: musc = {'CGS', 'MGS', 'MGS_SVL', 'MGS_LTS',...
%               'MGS_ICWY', 'MGS_CWY'};
%
% When run without arguments, RUNTESTSINGVAL returns loss of orthogonality
% and residual plots for default settings.
%
% (c) Kathryn Lund, Charles University, 2020

%%
addpath(genpath('../main/'))                                                % path to main routines
fstr = 'sing_val_low_sync';

% Defaults for inputs
if nargin == 0
    XXdim = [1000, 200];
    singval = -(1:16);
    musc = {'CGS', 'MGS', 'MGS_SVL', 'MGS_LTS',... 
        'MGS_ICWY', 'MGS_CWY'};
    musc_str = {'CGS', 'MGS', 'MGS\_SVL (MGS2)', 'MGS\_LTS (Alg 4)',...
        'MGS\_ICWY (Alg 5)', 'MGS\_CWY (Alg 6)'};
elseif nargin == 1
    singval = -(1:16);
    musc = {'CGS', 'MGS', 'MGS_SVL', 'MGS_LTS',... 
        'MGS_ICWY', 'MGS_CWY'};
    musc_str = {'CGS', 'MGS', 'MGS\_SVL (MGS2)', 'MGS\_LTS (Alg 4)',...
        'MGS\_ICWY (Alg 5)', 'MGS\_CWY (Alg 6)'};
elseif nargin == 2
    musc = {'CGS', 'MGS', 'MGS_SVL', 'MGS_LTS',... 
        'MGS_ICWY', 'MGS_CWY'};
    musc_str = {'CGS', 'MGS', 'MGS\_SVL (MGS2)', 'MGS\_LTS (Alg 4)',...
        'MGS\_ICWY (Alg 5)', 'MGS\_CWY (Alg 6)'};
end

% Defaults for empty arguments
if isempty(XXdim)
    XXdim = [1000, 200];
end
if isempty(singval)
    singval = -(1:16);
end
if isempty(musc)
    musc = {'CGS', 'MGS', 'MGS_SVL', 'MGS_LTS',... 
        'MGS_ICWY', 'MGS_CWY'};
    musc_str = {'CGS', 'MGS', 'MGS\_SVL (MGS2)', 'MGS\_LTS (Alg 4)',...
        'MGS\_ICWY (Alg 5)', 'MGS\_CWY (Alg 6)'};
end    

% Defaults for processing a single char array
if ischar(musc)
    musc = {musc};
end

% Default strings and replace underscore with tex underscore
if ~exist('musc_str','var')
    musc_str = musc;
    musc_str = strrep(musc_str, '_RO', '+');
    musc_str = strrep(musc_str, 'RO', '+');
    musc_str = strrep(musc_str, '_', '\_');
end

% Pre-allocate memory for measures
nmat = length(singval);
nmusc = length(musc);
loss_ortho = zeros(nmat, nmusc);
XXcond = zeros(1,nmat);

% Extract dimensions
m = XXdim(1); s = XXdim(2);
I = eye(s);

% Plot settings
musc_cmap = lines(nmusc);
musc_lbl = {'s-', 'o-', '*-', '^-', 'p-', '.-', 'h-', 'd-'};

[U,~] = qr(randn(m,s),0);
[V,~] = qr(randn(s,s),0);
for i = 1:nmat
    Sigma = diag(logspace(0,singval(i),s)');
    XX = U*Sigma*V';
    XXcond(i) = cond(XX);
    
        for k = 1:nmusc
            % Call BGS skeleton-muscle configuration
            QQ = IntraOrtho(XX, musc{k});

            % Compute loss of orthonormality
            loss_ortho(i, k) = norm(I - QQ'*QQ, 2);
            
            % Clear computed variables before next run
            clear QQ RR
        end
end

% Plot
x = XXcond;
lgd_str = musc_str;
fig_loss_ortho = clf; ax_loss_ortho = gca; hold on;
for k = 1:nmusc
    plot(ax_loss_ortho, x, loss_ortho(:,k),...
        musc_lbl{k}, 'Color', musc_cmap(k,:));
end
plot(ax_loss_ortho, x, eps*x, 'k--', x, eps*(x.^2), 'k-')
set(ax_loss_ortho, 'Yscale', 'log', 'Xscale', 'log');
title(ax_loss_ortho, 'Loss of Orthogonality ');
xlabel(ax_loss_ortho, '\kappa(X)')

lgd_str{end+1} = 'O(\epsilon) \kappa(X)';
lgd_str{end+1} = 'O(\epsilon) \kappa(X)^2';
legend(ax_loss_ortho, lgd_str, 'Location', 'BestOutside');

% Save plots
folderstr = sprintf('results/%s_m%d_s%d', fstr, m, s);
mkdir(folderstr)

savestr = sprintf('%s/out', folderstr);
save(savestr,'loss_ortho');

savestr = sprintf('%s/loss_ortho', folderstr);
savefig(fig_loss_ortho, savestr, 'compact');
saveas(fig_loss_ortho, savestr, 'epsc')

% close all;
end
 function RunTestSingValLowSyncBlock(XXdim, singval, skel, musc)
% RUNTESTSINGVALLOWSYNCBLOCK(XXdim, singval, skel, musc) is a wrapper
% function that compares stability for different skeleton-muscle
% combinations for a set of matrices of size XXdim = [m p s] with varying
% singular values specified by the vector array singval.
%
% skel and musc should be given as either a char array or a cell of char
% arrays (i.e., text strings with single quotes).
%
% Options for XXdim:
%   XXdim = [m p s], where m is the number of rows, p is the number of
%   block vectors, and s is the number of columns per block vector.
%
%   Default: XXdim = [2000 3 5]
%
% Options for singval:
%   singval should be a vector of negative powers.  Each entry of singval
%   is treated as the power of the smallest singular value for a matrix
%   whose singular values range from 10^0 to 10^t, where t is an entry of
%   singval.  For each t, a different matrix is generated, in order to see
%   how BGS behaves for matrices with different condition numbers.
%   
%   Default: singval = 1:16
%
% Options for skel: see BGS
%
%   Default: skel = BMGS_SVL
%
% Options for musc: see INTRAORTHO
%
%   Default: musc = {'CGS', 'MGS', 'MGS_SVL', 'MGS_LTS',...
%               'MGS_ICWY', 'MGS_CWY', 'HouseQR'};
%
% When run without arguments, RUNTESTSINGVAL returns loss of orthogonality
% and residual plots for default settings.
%
% (c) Kathryn Lund, Charles University, 2020

%%
addpath(genpath('../main/'))                                                % path to main routines
fstr = 'sing_val_low_sync_block';

% Defaults for inputs
if nargin == 0
    XXdim = [2000 3 5];
    singval = 1:16;
    skel = 'BMGS_SVL';
    musc = {'CGS', 'MGS', 'MGS_SVL', 'MGS_LTS',... 
        'MGS_ICWY', 'MGS_CWY', 'HouseQR'};
elseif nargin == 1
    singval = 1:16;
    skel = 'BMGS_SVL';
    musc = {'CGS', 'MGS', 'MGS_SVL', 'MGS_LTS',... 
        'MGS_ICWY', 'MGS_CWY', 'HouseQR'};
elseif nargin == 2
    skel = 'BMGS_SVL';
    musc = {'CGS', 'MGS', 'MGS_SVL', 'MGS_LTS',... 
        'MGS_ICWY', 'MGS_CWY', 'HouseQR'};
elseif nargin == 3
    musc = {'CGS', 'MGS', 'MGS_SVL', 'MGS_LTS',... 
        'MGS_ICWY', 'MGS_CWY', 'HouseQR'};
end

% Defaults for empty arguments
if isempty(XXdim)
    XXdim = [2000 3 5];
end
if isempty(singval)
    singval = 1:16;
end
if isempty(skel)
    skel = 'BMGS_SVL';
end
if isempty(musc)
    musc = {'CGS', 'MGS', 'MGS_SVL', 'MGS_LTS',... 
        'MGS_ICWY', 'MGS_CWY', 'HouseQR'};
end    

% Defaults for processing a single char array
if ischar(skel)
    skel = {skel};
end
if ischar(musc)
    musc = {musc};
end

% Default strings and replace underscore with tex underscore
if ~exist('skel_str','var')
    skel_str = skel;
    skel_str = strrep(skel_str, '_RO', '+');
    skel_str = strrep(skel_str, 'RO', '+');
    skel_str = strrep(skel_str, '_', '\_');
end
if ~exist('musc_str','var')
    musc_str = musc;
    musc_str = strrep(musc_str, '_RO', '+');
    musc_str = strrep(musc_str, 'RO', '+');
    musc_str = strrep(musc_str, '_', '\_');
end

% Pre-allocate memory for measures
nmat = length(singval);
nskel = length(skel);
nmusc = length(musc);
loss_ortho = zeros(nmat, nskel, nmusc);
res = zeros(nmat, nskel, nmusc);
XXcond = zeros(1,nmat);

% Extract dimensions
m = XXdim(1); p = XXdim(2); s = XXdim(3);
n = p*s;
I = eye(n);

% Run through matrices for each skel-musc combo
[U,~] = qr(randn(m,n),0);
[V,~] = qr(randn(n,n),0);
for i = 1:nmat
    Sigma = diag(logspace(0,singval(i),n)');
    XX = U*Sigma*V';
    XXcond(i) = cond(XX);
    XXnorm = norm(XX);
    
    for j = 1:nskel
        for k = 1:nmusc
            if strcmpi(skel{j}, 'bcgs_sror')
                if strcmpi(musc{k}, 'cgs_sror')
                    [QQ, RR] = BGS(XX, s, 'bcgs_sror', 'cgs_sror', rpltol);
                elseif strcmpi(musc{k}, 'cgs_sro')
                    [QQ, RR] = BGS(XX, s, 'bcgs_sror', 'cgs_sro', 0);
                else
                    QQ = NaN;
                    RR = NaN;
                end
            else
                % Call BGS skeleton-muscle configuration
                [QQ, RR] = BGS(XX, s, skel{j}, musc{k});
            end

            % Compute loss of orthonormality
            loss_ortho(i, j, k) = norm(I - QQ'*QQ, 2);
            
            % Compute relative residual
            res(i, j, k) = norm(XX - QQ*RR, 2)/XXnorm;

            % Clear computed variables before next run
            clear QQ RR            
        end
    end
end

% Plot
skel_cmap = lines(nskel);
musc_lbl = {'s-', 'o-', '*-', '^-', 'p-', '.-', 'h-', 'd-'};

x = XXcond;
lgd_str = {};
fig_loss_ortho = clf; ax_loss_ortho = gca; hold on;
% fig_res = figure; ax_res = gca; hold on;
for j = 1:nskel
    for k = 1:nmusc
        plot(ax_loss_ortho, x, loss_ortho(:,j,k),...
            musc_lbl{k}, 'Color', skel_cmap(j,:));
%         plot(ax_res, x, res(:,j,k),...
%             musc_lbl{k}, 'Color', skel_cmap(j,:));
        lgd_str{end+1} = sprintf('%s \\circ %s, s = %d',...
            skel_str{j}, musc_str{k}, s);
    end
end
plot(ax_loss_ortho, x, eps*x, 'k--', x, eps*(x.^2), 'k-')
set(ax_loss_ortho, 'Yscale', 'log', 'Xscale', 'log');
title(ax_loss_ortho, 'Loss of Orthogonality');
xlabel(ax_loss_ortho, '\kappa(X)')
set(ax_loss_ortho, 'XGrid', 'on', 'YGrid', 'on',...
    'XMinorGrid', 'off', 'YMinorGrid', 'off');

% set(ax_res, 'Yscale', 'log', 'Xscale', 'log');
% title(ax_res, 'Relative Residual');
% xlabel(ax_res, '\kappa(X)')
% set(ax_res, 'XGrid', 'on', 'YGrid', 'on',...
%     'XMinorGrid', 'off', 'YMinorGrid', 'off');
% legend(ax_res, lgd_str, 'Location', 'BestOutside');

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
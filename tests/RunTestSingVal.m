function RunTestSingVal(XXdim, singval, skel, musc)
% RUNTESTSINGVAL(XXdim, singval, musc, verbose) is a wrapper function that
% executes a series of Block Gram-Schmidt (BGS) variants for a series of
% matrices specified by singval via the skeleton-muscle paradigm and
% returns loss of orthogonality and residual plots.  In more detail,
% RUNTESTSINGVAL
% 1. executes a loop that loads/generates the matrix or matrices specified
%    by singval for the dimensions [m p s] = XXdim;
% 2. executes a second loop running through skeleton options specfied by
%    skel;
% 3. executes a third loop running through muscle options specified by
%    musc;
% 4. plots and saves results.
%
% skel and musc should be given as either a char array or a cell of char
% arrays (i.e., text strings with single quotes).
%
% Options for XXdim:
%   XXdim = [m p s], where m is the number of rows, p is the number of
%   block vectors, and s is the number of columns per block vector.  The
%   resulting matrix XX has dimensions m x p*s.
%
%   Default: XXdim = [1000 50 5]
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
% Options for skel (see also BGS):
%   'BCGS' - Block Classical Gram-Schmidt
%   'BCGS_IRO' - Block Classical Gram-Schmidt with Inner
%       ReOrthonormalization
%   'BMGS' - Block Modified Gram-Schmidt
%
%   Default: skel = {'BCGS', 'BMGS', 'BCGS_IRO'}
%
% Options for musc (see also IntraOrtho):
%   'CGS' - Classical Gram-Schmidt.
%   'CGS_RO' - Classical Gram-Schmidt with ReOrthonormalization. Calls CGS
%       twice.
%   'CGS_IRO' - Classical Gram-Schmidt with Inner ReOrthonormalization.
%   'MGS' - Modified Gram-Schmidt.
%   'MGS_RO' - Modified Gram-Schmidt with ReOrthonormalization. Cals MGS
%      twice.
%   'MGS_IRO' - Modified Gram-Schmidt with Inner ReOrthonormalization.
%   'HouseQR' - Householder-based QR. Calls built-in QR with economic
%       setting.
%
%   Default: musc = {'CGS', 'MGS', 'HouseQR'}
%
% When run without arguments, RUNTESTSINGVAL returns loss of orthogonality
% and residual plots for default settings.
%
% When specifying only a subset of arguments, set non-specified arguments
% to [] to ensure the default value.  For example,
%
%     RunTestSingVal([2000 4 5], [], [], 'HouseQR')
%
% runs tests for matrices with dimensions [2000 4 5], default singval and
% skel settings, and only the IntraOrtho HouseQR.
%
% (c) Kathryn Lund, Charles University, 2020

%%
addpath(genpath('../main/'))                                                % path to main routines
fstr = 'sing_val';

% Defaults for inputs
if nargin == 0
    XXdim = [1000 50 5];
    singval = -(1:16);
    skel = {'BCGS', 'BCGS_PIP', 'BCGS_PIO', 'BMGS', 'BCGS_IRO'};
    skel_str = {'BCGS', 'BCGS-PIP', 'BCGS-PIO', 'BMGS', 'BCGS2'};
    musc = {'CGS', 'MGS', 'HouseQR'};
    musc_str = {'CGS', 'MGS', 'HouseQR'};
elseif nargin == 1
    singval = -(1:16);
    skel = {'BCGS', 'BCGS_PIP', 'BCGS_PIO', 'BMGS', 'BCGS_IRO'};
    skel_str = {'BCGS', 'BCGS-PIP', 'BCGS-PIO', 'BMGS', 'BCGS2'};
    musc = {'CGS', 'MGS', 'HouseQR'};
    musc_str = {'CGS', 'MGS', 'HouseQR'};
elseif nargin == 2
    skel = {'BCGS', 'BCGS_PIP', 'BCGS_PIO', 'BMGS', 'BCGS_IRO'};
    skel_str = {'BCGS', 'BCGS-PIP', 'BCGS-PIO', 'BMGS', 'BCGS2'};
    musc = {'CGS', 'MGS', 'HouseQR'};
    musc_str = {'CGS', 'MGS', 'HouseQR'};
elseif nargin == 3
    musc = {'CGS', 'MGS', 'HouseQR'};
    musc_str = {'CGS', 'MGS', 'HouseQR'};
end

% Defaults for empty arguments
if isempty(XXdim)
    XXdim = [1000 50 5];
end
if isempty(singval)
    singval = -(1:16);
end
if isempty(skel)
    skel = {'BCGS', 'BCGS_PIP', 'BCGS_PIO', 'BMGS', 'BCGS_IRO'};
    skel_str = {'BCGS', 'BCGS-PIP', 'BCGS-PIO', 'BMGS', 'BCGS2'};
end
if isempty(musc)
    musc = {'CGS', 'MGS', 'HouseQR'};
    musc_str = {'CGS', 'MGS', 'HouseQR'};
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

% Extract dimensions
m = XXdim(1); p = XXdim(2); s = XXdim(3);
n = p*s;
I = eye(n);
XXnorm = zeros(1,nmat);
XXcond = zeros(1,nmat);

U = orth(randn(m,n));
V = orth(randn(n,m));
for i = 1:nmat
    Sigma = diag(logspace(0,singval(i),n)');
    XX = U*Sigma*V';
    XXnorm(i) = norm(XX);
    XXcond(i) = cond(XX);
    
    for j = 1:nskel
        for k = 1:nmusc
            % Call BGS skeleton-muscle configuration
            [QQ, RR] = BGS(XX, s, skel{j}, musc{k});

            % Compute loss of orthonormality
            loss_ortho(i, j, k) = norm(I - QQ'*QQ, 2);

%             % Compute relative residual
%             res(i, j, k) = norm(XX - QQ*RR, 2)/XXnorm(i);
            
            % Compute relative residual for Cholesky relation
            res(i, j, k) = norm(XX'*XX - RR'*RR, 2)/XXnorm(i)^2;
            
            % Clear computed variables before next run
            clear QQ RR            
        end
    end
end

%% Plot
fig_loss_ortho = figure; ax_loss_ortho = gca; hold on;
fig_res = figure; ax_res = gca; hold on;

% Position options for nice viewing and paper image size; comment if plots
% do not render correctly on your machine
set(fig_loss_ortho, 'Position', [11 158 645 420])
set(fig_res, 'Position', [613 158 645 420])

skel_cmap = lines(nskel);
musc_lbl = {'s-', 'o-', '*-', 'p-', 'h-', '.-', '^-'};

x = XXcond;
lgd_str = {};
for j = 1:nskel
    for k = 1:nmusc
        plot(ax_loss_ortho, x, loss_ortho(:,j,k),...
            musc_lbl{k}, 'Color', skel_cmap(j,:));
        plot(ax_res, x, res(:,j,k), ... 
            musc_lbl{k}, 'Color', skel_cmap(j,:));
        lgd_str{end+1} = sprintf('%s \\circ %s', skel_str{j}, musc_str{k});
    end
end
plot(ax_loss_ortho, x, eps*x, 'k--', x, eps*(x.^2), 'k-')
plot(ax_res, x, eps*XXnorm, 'k--', x, eps*XXnorm.^2, 'k-')

set(ax_loss_ortho, 'Yscale', 'log', 'Xscale', 'log','XGrid', 'on', 'YGrid', 'on',...
    'XMinorGrid', 'off', 'YMinorGrid', 'off');
ylabel(ax_loss_ortho, 'loss of orthogonality ');
xlabel(ax_loss_ortho, '\kappa(X)')
lgd_str{end+1} = '\epsilon_M \kappa(X)';
lgd_str{end+1} = '\epsilon_M \kappa(X)^2';
legend(ax_loss_ortho, lgd_str, 'Location', 'BestOutside');

set(ax_res, 'Yscale', 'log', 'Xscale', 'log','XGrid', 'on', 'YGrid', 'on',...
    'XMinorGrid', 'off', 'YMinorGrid', 'off');
ylabel(ax_res, 'relative residual');
xlabel(ax_res, '\kappa(X)')
lgd_str{end-1} = '\epsilon_M ||X||';
lgd_str{end} = '\epsilon_M ||X||^2';
legend(ax_res, lgd_str, 'Location', 'BestOutside');

% Save plots
folderstr = sprintf('results/%s_m%d_p%d_s%d', fstr, m, p, s);
mkdir(folderstr)

savestr = sprintf('%s/out', folderstr);
save(savestr, 'x', 'loss_ortho','res');

savestr = sprintf('%s/loss_ortho', folderstr);
savefig(fig_loss_ortho, savestr, 'compact');
saveas(fig_loss_ortho, savestr, 'epsc')

savestr = sprintf('%s/res', folderstr);
savefig(fig_res, savestr, 'compact');
saveas(fig_res, savestr, 'epsc')

% close all;
end
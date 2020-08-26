function RunTestMonomial(XXdim, svec, skel, musc)
% RUNTESTMONOMIAL(XXdim, svec, musc, verbose) is a wrapper function that
% executes a series of Block Gram-Schmidt (BGS) variants for a series of
% matrices specified by svec via the skeleton-muscle paradigm and
% returns loss of orthogonality and residual plots.  In more detail,
% RUNTESTMONOMIAL
% 1. executes a loop that loads/generates the matrix or matrices specified
%    by svec for the dimensions [m n] = XXdim;
% 2. executes a second loop running through skeleton options specfied by
%    skel;
% 3. executes a third loop running through muscle options specified by
%    musc;
% 4. plots and saves results.
%
% skel and musc should be given as either a char array or a cell of char
% arrays (i.e., text strings with single quotes).
%
% Options for XXdim and svec:
%   XXdim = [m n], where m is the number of rows and n the number of
%   columns. svec is a vector of scalars that divide n.  With s denoting
%   one such scalar, a matrix XX is produced with m rows and p = n/s block
%   vectors each with s columns.  Each block vector X_k of XX is built from
%   a random vector v_k like
%
%       X_k = [v_k A*v_k ... A^(s-1)*v_k].
%
%   Defaults:
%       XXdim = [5000 840];
%       svec = 2:2:14;
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
% When run without arguments, RUNTESTMONOMIAL returns loss of orthogonality
% and residual plots for default settings.
%
% When specifying only a subset of arguments, set non-specified arguments
% to [] to ensure the default value.  For example,
%
%     RunTestSingVal([1000 200], [], [], 'HouseQR')
%
% runs tests for matrices with dimensions [1000 200], default svec and
% skel settings, and only the IntraOrtho HouseQR.
%
% (c) Kathryn Lund, Charles University, 2020

%%
addpath(genpath('../main/'))                                                % path to main routines
fstr = 'monomial';

% Defaults for inputs
if nargin == 0
    XXdim = [5000 840];
    svec = 2:2:14;
    skel = {'BCGS', 'BCGS_PIP', 'BCGS_PIO', 'BMGS', 'BCGS_IRO'};
    skel_str = {'BCGS', 'BCGS-PIP', 'BCGS-PIO', 'BMGS', 'BCGS2'};
    musc = {'CGS', 'MGS', 'HouseQR'};
    musc_str = {'CGS', 'MGS', 'HouseQR'};
elseif nargin == 1
    svec = 2:2:14;
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
    XXdim = [5000 840];
end
if isempty(svec)
    svec = 2:2:14;
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
nmat = length(svec);
nskel = length(skel);
nmusc = length(musc);
loss_ortho = zeros(nmat, nskel, nmusc);
res = zeros(nmat, nskel, nmusc);
XXnorm = zeros(1,nmat);
XXcond = zeros(1,nmat);

% Extract dimensions
m = XXdim(1); n = XXdim(2);
I = eye(n);

A = spdiags(linspace(.1,1,m)',0,m,m);
for i = 1:nmat
    % Create or load XX
    s = svec(i); p = n/s;
    matstr = sprintf('%s_m%d_p%d_s%d.mat', fstr, m, p, s);
    
    cd matrices
    if exist(matstr, 'file')
        load(matstr, 'XX');
    else
        XXhat = rand(m,p);
        pp = 1:p;
        XXhat(:,pp) = XXhat(:,pp)/norm(XXhat(:,pp));
        for k = 2:s
            pp = pp + p;
            XXhat(:,pp) = A*XXhat(:,pp - p);
        end
        % Reshape XX
        XX = zeros(m,n);
        ind = 1:s:n;
        kk = 1:p;
        for k = 1:s
            XX(:,ind) = XXhat(:,kk);
            ind = ind + 1;
            kk = kk + p;
        end
        
        save(matstr, 'XX');
    end
    cd ..
    XXnorm(i) = norm(XX);
    XXcond(i) = cond(XX);
    
    for j = 1:nskel
        for k = 1:nmusc
            % Call BGS skeleton-muscle configuration
            [QQ, RR] = BGS(XX, s, skel{j}, musc{k});

            % Compute loss of orthonormality
            loss_ortho(i, j, k) = norm(I - QQ'*QQ, 2);

%             % Compute relative residual
%             res(i, j, k) = norm(XX - QQ*RR, 2)/XXnorm;

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
folderstr = sprintf('results/%s_m%d_n%d', fstr, m, n);
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
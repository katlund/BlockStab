function GluedBlockKappaPlot(XXdim, logcondXX, skel, musc)
% RUNTESTGLUED(XXdim, logcondXX, skel, s_bgs) is a wrapper function
% that executes a series of Block Gram-Schmidt (BGS) variants for a series
% of matrices specified by singval via the skeleton-muscle paradigm and
% returns loss of orthogonality and residual plots.  In more detail,
% RUNTESTGLUED
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
%   singval should be a vector of positive powers, where the power
%   denotes the log of the desired condition number of XX.
%   
%   Default: singval = 1:8
%
% Options for skel (see also BGS):
%   'BCGS' - Block Classical Gram-Schmidt
%   'BCGS_IRO' - Block Classical Gram-Schmidt with Inner
%       ReOrthonormalization
%   'BMGS' - Block Modified Gram-Schmidt
%   'BCGS_PIP' - BCGS with Pythagorean Inner Product modification
%   'BCGS_PIO' - BCGS with Pythagorean Intra-Orthogonalization modification
%
%   Default: skel = {'BCGS', 'BCGS_PIP', 'BCGS_PIO'}
%
% Options for musc:
%   See INTRAORTHO.
%
%   Default: musc = {'CGS', 'MGS', 'HouseQR'};

%%
addpath(genpath('../main/'))                                                % path to main routines
fstr = 'glued';

% Defaults for inputs
if nargin == 0
    XXdim = [1000 50 5];
    logcondXX = 1:8;
    skel = {'BCGS', 'BCGS_PIP', 'BCGS_PIO', 'BMGS', 'BCGS_IRO'};
    skel_str = {'BCGS', 'BCGS-PIP', 'BCGS-PIO', 'BMGS', 'BCGS2'};
    musc = {'CGS', 'MGS', 'HouseQR'};
elseif nargin == 1
    logcondXX = 1:8;
    skel = {'BCGS', 'BCGS_PIP', 'BCGS_PIO', 'BMGS', 'BCGS_IRO'};
    skel_str = {'BCGS', 'BCGS-PIP', 'BCGS-PIO', 'BMGS', 'BCGS2'};
    musc = {'CGS', 'MGS', 'HouseQR'};
elseif nargin == 2
    skel = {'BCGS', 'BCGS_PIP', 'BCGS_PIO', 'BMGS', 'BCGS_IRO'};
    skel_str = {'BCGS', 'BCGS-PIP', 'BCGS-PIO', 'BMGS', 'BCGS2'};
    musc = {'CGS', 'MGS', 'HouseQR'};
elseif nargin == 3
    musc = {'CGS', 'MGS', 'HouseQR'};
end

% Defaults for empty arguments
if isempty(XXdim)
    XXdim = [1000 50 5];
end
if isempty(logcondXX)
    logcondXX = 1:8;
end
if isempty(skel)
    skel = {'BCGS', 'BCGS_PIP', 'BCGS_PIO', 'BMGS', 'BCGS_IRO'};
    skel_str = {'BCGS', 'BCGS-PIP', 'BCGS-PIO', 'BMGS', 'BCGS2'};
end
if isempty(musc)
    musc = {'CGS', 'MGS', 'HouseQR'};
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
% Default strings and replace underscore with tex underscore
if ~exist('musc_str','var')
    musc_str = musc;
    musc_str = strrep(musc_str, '_RO', '+');
    musc_str = strrep(musc_str, 'RO', '+');
    musc_str = strrep(musc_str, '_', '\_');
end

% Pre-allocate memory for measures
nmat = length(logcondXX);
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

for i = 1:nmat
    % Create glued matrix; if s matches glued block size, results are
    % skewed; see RUNTESTGLUEDVARYS
    sglued = p;
    pglued = s;
    matstr = sprintf('%s_cond%d_m%d_p%d_s%d.mat', fstr, logcondXX(i), m, pglued, sglued);
    cd matrices
    if exist(matstr, 'file')
        load(matstr, 'XX')
    else
        XX = MatGen([m, pglued, sglued], .5*logcondXX(i), logcondXX(i));
        
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

            % Compute relative residual
            res(i, j, k) = norm(XX - QQ*RR, 2)/XXnorm(i);
            
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

% Position options for nice viewing and paper image size; comment out if
% plots do not render correctly
set(fig_loss_ortho, 'Position', [11 158 645 420])
set(fig_res, 'Position', [613 158 645 420])

skel_cmap = lines(nskel);
musc_lbl = {'s-', 'o-', '*-', 'p-', 'h-', '.-', '^-'};
lgd_str = {};

x = XXcond; % condition number 
for j = 1:nskel
    for k = 1:nmusc
        plot(ax_loss_ortho, x, loss_ortho(:,j,k),...
            musc_lbl{k}, 'Color', skel_cmap(j,:));
        plot(ax_res, x, res(:,j,k), ... 
            musc_lbl{k}, 'Color', skel_cmap(j,:));
        if s == 1
            lgd_str{end+1} = sprintf('%s',...
                skel_str{j}(2:end));
        else
            lgd_str{end+1} = sprintf('%s \\circ %s',...
                skel_str{j}, musc_str{k});
        end
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
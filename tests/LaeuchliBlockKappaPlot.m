function LaeuchliBlockKappaPlot(XXdim, eps_vec, skel, musc)
% BLOCKKAPPAPLOT(XXdim, eps_vs, skel, musc) compares loss of orthogonality
% and relative residual for different skeleton-muscle combinationss for a
% set of matrices of size XXdim = [m p s] with varying condition numbers
% determined by the vector array eps_vec.
%
% skel and musc should be given as either a char array or a cell of char
% arrays (i.e., text strings with single quotes).
%
% Options for XXdim:
%   XXdim = [m p s], where m is the number of rows, p is the number of
%   block vectors, and s is the number of columns per block vector.  The
%   resulting matrix XX has dimensions m x p*s.
%
%   Default: XXdim = [100 10 2]
%
% Options for eps_vec:
%   eps_vec should be a vector of values strictly between 0 and 1
%   
%   Default: eps_vec = logspace(-1, -16, 10);
%
% Options for skel: see BGS.
%
% Options for musc: see INTRAORTHO.
%
% When specifying only a subset of arguments, set non-specified arguments
% to [] to ensure the default value.  For example,
%
%     LaeuchliBlockKappaPlot([], [], 'BCGS', 'HouseQR')
%
% runs tests for matrices with dimensions [100 20 2], default logcondXX
% settings, and the skeleton-muscle combination BCGS \circ HouseQR.
%
% Note that the routine may take a while for large XXdim and many skel/musc
% options, since all possible combinations of skel/musc are performed.  For
% example,
%
%     LaeuchliBlockKappaPlot([], [], {'BCGS', 'BCGS_PIP', 'BCGS_PIO'}, {'CGS', 'MGS', 'HouseQR'})
%
% runs 9 BGS combinations in total, given the 3 skeletons and 3 muscles.

%%
addpath(genpath('../main/'))                                                % path to main routines
addpath(genpath('auxiliary/'))
addpath(genpath('matrices/'))
fstr = 'laeuchli_block_kappa_plot';

% Defaults for empty arguments
if isempty(XXdim)
    XXdim = [100 20 2];
end
if isempty(eps_vec)
    eps_vec = logspace(-1, -16, 10);
end  

% Defaults for processing a single char array
if ischar(skel)
    skel = {skel};
end
if ischar(musc)
    musc = {musc};
end

% Default strings and format strings for plot legends
skel_str = AlgString(skel);
musc_str = AlgString(musc);

% Pre-allocate memory for measures
nmat = length(eps_vec);
nskel = length(skel);
nmusc = length(musc);
loss_ortho = zeros(nmat, nskel, nmusc);
res = zeros(nmat, nskel, nmusc);
res_chol = zeros(nmat, nskel, nmusc);
XXcond = zeros(1,nmat);
XXnorm = zeros(1,nmat);

% Extract dimensions
m = XXdim(1); p = XXdim(2); s = XXdim(3);
n = p*s;
I = eye(n);

for i = 1:nmat
    XX = laeuchli(m, n, eps_vec(i));
    XXcond(i) = cond(XX);
    XXnorm(i) = norm(XX, 2);
    
    for j = 1:nskel
        for k = 1:nmusc
            % Call BGS skeleton-muscle configuration
            [QQ, RR] = BGS(XX, s, skel{j}, musc{k});

            % Compute loss of orthonormality
            loss_ortho(i, j, k) = norm(I - QQ' * QQ, 2);

            % Compute relative residual
            res(i, j, k) = norm(XX - QQ * RR, 2) / XXnorm(i);
            
            % Compute relative residual for Cholesky relation
            res_chol(i, j, k) = norm(XX' * XX - RR' * RR, 2) / XXnorm(i)^2;
            
            % Clear computed variables before next run
            clear QQ RR            
        end
    end
end

%% Plots
skel_cmap = lines(nskel);
musc_lbl = {'s-', 'o-', '*-', '^-', 'p-', '.-', 'h-', 'd-'};

x = XXcond;
lgd_str = {};

% Initialize figures and axes
fg = cell(1,3); ax = cell(1,3);
for i = 1:3
    fg{i} = figure;
    ax{i} = gca;
    hold on;
end

% Plot data
for j = 1:nskel
    for k = 1:nmusc
        plot(ax{1}, x, loss_ortho(:,j,k),...
            musc_lbl{k}, 'Color', skel_cmap(j,:),'MarkerSize',10,'LineWidth',1);
        plot(ax{2}, x, res(:,j,k),... 
            musc_lbl{k}, 'Color', skel_cmap(j,:),'MarkerSize',10,'LineWidth',1);
        plot(ax{3}, x, res_chol(:,j,k),...
            musc_lbl{k}, 'Color', skel_cmap(j,:),'MarkerSize',10,'LineWidth',1);
        lgd_str{end+1} = sprintf('%s $\\circ$ %s', skel_str{j}, musc_str{k});
    end
end
plot(ax{1}, x, eps*x, 'k--', x, eps*(x.^2), 'k-')

% Make plots pretty and save figures
folder_str = sprintf('results/%s_m%d_p%d_s%d', fstr, m, p, s);
mkdir(folder_str)
PrettyKappaPlot(fg, ax, lgd_str, folder_str);

% Save data
savestr = sprintf('%s/out', folder_str);
save(savestr, 'x', 'loss_ortho','res');

% close all;
end
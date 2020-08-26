function BlockKappaPlot(XXdim, logcondXX, skel, musc)
% BLOCKKAPPAPLOT(XXdim, logcondXX, skel, musc) compares loss of
% orthogonality and relative residual for different skeleton-muscle
% combinationss for a set of matrices of size XXdim = [m p s] with varying
% singular values specified by the vector array logcondXX.
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
% Options for logcondXX:
%   logcondXX should be a vector of negative powers.  Each entry of
%   logcondXX is treated as the power of the smallest singular value for a
%   matrix whose singular values range from 10^0 to 10^t, where t is an
%   entry of logcondXX.  For each t, a different matrix is generated, in
%   order to see how BGS behaves for matrices with different condition
%   numbers.
%   
%   Default: logcondXX = -(1:16)
%
% Options for skel: see BGS
%
% Options for musc: see INTRAORTHO
%
% When specifying only a subset of arguments, set non-specified arguments
% to [] to ensure the default value.  For example,
%
%     BlockKappaPlot([100 20 2], [], 'BCGS', 'HouseQR')
%
% runs tests for matrices with dimensions [100 20 2], default logcondXX
% settings, and the skeleton-muscle combination BCGS \circ HouseQR.
%
% Note that the routine may take a while for large XXdim and many skel/musc
% options, since all possible combinations of skel/musc are performed.  For
% example,
%
%     BlockKappaPlot([100 20 2], [], {'BCGS', 'BCGS_PIP', 'BCGS_PIO'}, {'CGS', 'MGS', 'HouseQR'})
%
% runs 9 BGS combinations in total, given the 3 skeletons and 3 muscles.

%%
addpath(genpath('../main/'))                                                % path to main routines
fstr = 'block_kappa_plot';

% Defaults for empty arguments
if isempty(XXdim)
    XXdim = [1000 50 5];
end
if isempty(logcondXX)
    logcondXX = -(1:16);
end  

% Defaults for processing a single char array
if ischar(skel)
    skel = {skel};
end
if ischar(musc)
    musc = {musc};
end

% Default strings and replace underscore with tex underscore
skel_str = AlgString(skel);
musc_str = AlgString(musc);

% Pre-allocate memory for measures
nmat = length(logcondXX);
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

U = orth(randn(m,n));
V = orth(randn(n,n));
for i = 1:nmat
    Sigma = diag(logspace(0, logcondXX(i), n)');
    XX = U * Sigma * V';
    XXnorm(i) = norm(XX);
    XXcond(i) = cond(XX);
    
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
musc_lbl = {'s-', 'o-', '*-', 'p-', 'h-', '.-', '^-'};

x = XXcond;
lgd_str = {};

fg = cell(1,3); ax = cell(1,3);
for i = 1:3
    fg{i} = figure;
    ax{i} = gca;
    hold on;
end

for j = 1:nskel
    for k = 1:nmusc
        plot(ax{1}, x, loss_ortho(:,j,k),...
            musc_lbl{k}, 'Color', skel_cmap(j,:));
        plot(ax{2}, x, res(:,j,k),... 
            musc_lbl{k}, 'Color', skel_cmap(j,:));
        plot(ax{3}, x, res_chol(:,j,k),...
            musc_lbl{k}, 'Color', skel_cmap(j,:));
        lgd_str{end+1} = sprintf('%s \\circ %s', skel_str{j}, musc_str{k});
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
function BlockKappaPlot(XXdim, logcondXX, skel, musc)
% BLOCKKAPPAPLOT(XXdim, logcondXX, musc, verbose) compares loss of
% orthogonality and relative residual for different skeleton-muscle
% combinations for a set of matrices of size XXdim = [m p s] with varying
% singular values specified by the vector array logcondX.
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
%     BlockKappaPlot([2000 4 5], [], 'BCGS', 'HouseQR')
%
% runs tests for matrices with dimensions [2000 4 5], default logcondXX
% settings, and the skeleton-muscle combination BCGS \circ HouseQR.

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

% Extract dimensions
m = XXdim(1); p = XXdim(2); s = XXdim(3);
n = p*s;
I = eye(n);
XXnorm = zeros(1,nmat);
XXcond = zeros(1,nmat);

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
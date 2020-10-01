 function LaeuchliKappaPlot(Xdim, eps_vec, musc)
% LAEUCHLIKAPPAPLOT(Xdim, eps_vec, musc) compares loss of orthogonality and
% relative residual for different muscles for a set of matrices of size
% Xdim = [m s] with varying condition numbers determined by the vector
% array eps_vec.
%
% musc should be given as either a char array or a cell of char arrays
% (i.e., text strings with single quotes).
%
% Options for Xdim:
%   Xdim = [m s], where m is the number of rows and s is the number of
%   columns.
%
%   Default: Xdim = [100 20]
%
% Options for eps_vec:
%   eps_vec should be a vector of values strictly between 0 and 1
%   
%   Default: eps_vec = logspace(-1, -16, 10);
%
% Options for musc: see INTRAORTHO.
%
% When specifying only a subset of arguments, set non-specified arguments
% to [] to ensure the default value.  For example,
%
%     LaeuchliKappaPlot([], [], 'CholQR')
%
% runs tests for matrices with dimensions [100 20], default eps_vec
% settings, and the muscle CholQR.

%%
addpath(genpath('../main/'))                                                % path to main routines
addpath(genpath('auxiliary/'))
addpath(genpath('matrices/'))
fstr = 'laeuchli_kappa_plot';

% Defaults for empty arguments
if isempty(Xdim)
    Xdim = [100, 20];
end
if isempty(eps_vec)
    eps_vec = logspace(-1, -16, 10);
end 

% Defaults for processing a single char array
if ischar(musc)
    musc = {musc};
end

% Default strings and format strings for plot legends
musc_str = AlgString(musc);

% Pre-allocate memory for measures
nmat = length(eps_vec);
nmusc = length(musc);
loss_ortho = zeros(nmat, nmusc);
res = zeros(nmat, nmusc);
res_chol = zeros(nmat, nmusc);
Xcond = zeros(1,nmat);
Xnorm = zeros(1,nmat);

% Extract dimensions
m = Xdim(1); s = Xdim(2);
I = eye(s);

for i = 1:nmat
    X = laeuchli(m, s, eps_vec(i));
    Xcond(i) = cond(X);
    Xnorm(i) = norm(X, 2);
    
    for k = 1:nmusc
        % Call IntraOrtho
        [Q, R] = IntraOrtho(X, musc{k});

        % Compute loss of orthonormality
        loss_ortho(i, k) = norm(I - Q' * Q, 2);

        % Compute relative residual
        res(i, k) = norm(X - Q * R, 2) / Xnorm(i);

        % Compute relative residual for Cholesky relation
        res_chol(i, k) = norm(X' * X - R' * R, 2) / Xnorm(i)^2;

        % Clear computed variables before next run
        clear Q R
    end
end

%% Plots
musc_cmap = lines(nmusc);
musc_lbl = {'s-', 'o-', '*-', '^-', 'p-', '.-', 'h-', 'd-'};

x = Xcond;
lgd_str = musc_str;

% Initialize figures and axes
fg = cell(1,3); ax = cell(1,3);
for i = 1:3
    fg{i} = figure;
    ax{i} = gca;
    hold on;
end

% Plot data
for k = 1:nmusc
    plot(ax{1}, x, loss_ortho(:,k),...
        musc_lbl{k}, 'Color', musc_cmap(k,:));
    plot(ax{2}, x, res(:,k),...
        musc_lbl{k}, 'Color', musc_cmap(k,:));
    plot(ax{3}, x, res_chol(:,k),...
        musc_lbl{k}, 'Color', musc_cmap(k,:));
end
plot(ax{1}, x, eps*x, 'k--', x, eps*(x.^2), 'k-')

% Make plots pretty and save figures
folder_str = sprintf('results/%s_m%d_s%d', fstr, m, s);
mkdir(folder_str)
PrettyKappaPlot(fg, ax, lgd_str, folder_str);

% Save data
save_str = sprintf('%s/out', folder_str);
save(save_str, 'x', 'loss_ortho', 'res', 'res_chol');

% close all;
end
function GluedBlockKappaPlot(XXdim, logcondXX, skel, musc)
% GLUEDBLOCKKAPPAPLOT(XXdim, logcondXX, skel, s_bgs) compares loss of
% orthogonality and relative residual for different skeleton-muscle
% combinationss for a set of glued matrices of size XXdim = [m p s] with
% varying condition numbers loosely specified by the exponents in
% logcondXX.
%
% skel and musc should be given as either a char array or a cell of char
% arrays (i.e., text strings with single quotes).
%
% Options for XXdim:
%   XXdim = [m p s], where m is the number of rows, p is the number of
%   block vectors, and s is the number of columns per block vector.  The
%   resulting matrix XX has dimensions m x p*s.
%
%   Default: XXdim = [1000 50 4]
%
% Options for logcondXX:
%   logcondXX should be a vector of positive powers, where the power
%   denotes the log of the desired condition number of XX.
%   
%   Default: logcondXX = 1:8
%
% Options for skel: see BGS
%
% Options for musc: see INTRAORTHO.
%
% See also BLOCKKAPPAPLOT for more details about basic functionalities.

%%
addpath(genpath('../main/'))                                                % path to main routines
addpath(genpath('auxiliary/'))
fstr = 'glued_block_kappa_plot';

% Defaults for inputs
if nargin == 0
    XXdim = [1000 50 4];
    logcondXX = 1:8;
elseif nargin == 1
    logcondXX = 1:8;
end

% Defaults for empty arguments
if isempty(XXdim)
    XXdim = [1000 50 4];
end
if isempty(logcondXX)
    logcondXX = 1:8;
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

% Fix glued dimensions
factors = factor(n);
mid_ind = round(length(n)/2);
r = prod(factors(1:mid_ind));
t = prod(factors(mid_ind+1:end));

for i = 1:nmat
    % Create glued matrix; if s matches glued block size, results are
    % skewed; see GLUEDBLOCKKAPPAPLOTVARYS
    matstr = sprintf('glued_cond%d_m%d_p%d_s%d.mat', logcondXX(i), m, r, t);
    cd matrices
    if exist(matstr, 'file')
        load(matstr, 'XX')
    else
        XX = CreateGluedMatrix(m, r, t,...
            .5 * logcondXX(i), logcondXX(i));
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
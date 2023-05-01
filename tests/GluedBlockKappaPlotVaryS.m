function GluedBlockKappaPlotVaryS(XXdim, logcondXX, s_bgs, skel)
% GLUEDBLOCKKAPPAPLOTVARYS(XXdim, logcondXX, skel, s_bgs, skel) compares
% loss of orthogonality and relative residual for different skeleton-muscle
% combinationss for a set of glued matrices of size XXdim = [m n] with
% varying condition numbers loosely specified by the exponents in
% logcondXX.
%
% skel should be given as either a char array or a cell of char arrays
% (i.e., text strings with single quotes).
%
% Options for XXdim:
%   XXdim = [m n], where m is the number of rows, and n is the number of
%   columns.
%
%   Default: XXdim = [1000 200]
%
% Options for logcondXX:
%   logcondXX should be a vector of positive powers, where the power
%   denotes the log of the desired condition number of XX.
%   
%   Default: logcondXX = 1:8
%
% Options for s_bgs:
%   s_bgs is a vector of block sizes for the call to BGS. This differs
%   from the block size used to generate the matrix XX.
%
%   Default: s_bgs = [1 2 5 10 20]
%
% Options for skel: see BGS.

%%
addpath(genpath('../main/'))                                                % path to main routines
addpath(genpath('auxiliary/'))
fstr = 'glued_block_kappa_plot_vary_s';

% Defaults for inputs
if nargin == 0
    XXdim = [1000 200];
    logcondXX = 1:8;
    s_bgs = [1 2 5 10 20];
elseif nargin == 1
    logcondXX = 1:8;
    s_bgs = [1 2 5 10 20];
elseif nargin == 2
    s_bgs = [1 2 5 10 20];
end

% Defaults for empty arguments
if isempty(XXdim)
    XXdim = [10000 200];
end
if isempty(logcondXX)
    logcondXX = 1:8;
end
if isempty(s_bgs)
    s_bgs = [1 2 5 10 20];
end    

% Default strings and format strings for plot legends
skel_str = AlgString(skel);

% Pre-allocate memory for measures
nmat = length(logcondXX);
nskel = length(skel);
ns_bgs = length(s_bgs);
loss_ortho = zeros(nmat, nskel, ns_bgs);
res = zeros(nmat, nskel, ns_bgs);
res_chol = zeros(nmat, nskel, ns_bgs);
XXcond = zeros(1,nmat);
XXnorm = zeros(1,nmat);

% Extract dimensions
m = XXdim(1); n = XXdim(2);
I = eye(n);

% Fix glued dimensions
factors = factor(n);
mid_ind = round(length(n)/2);
r = prod(factors(1:mid_ind));
t = prod(factors(mid_ind+1:end));

for i = 1:nmat
    % Create glued matrix
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
        for k = 1:ns_bgs
            % Call BGS skeleton-muscle configuration
            [QQ, RR] = BGS(XX, s_bgs(k), skel{j}, 'HouseQR');

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

%% Plot
skel_cmap = lines(nskel);
s_bgs_lbl = {'s-', 'o-', '*-', '^-', 'p-', '.-', 'h-', 'd-'};

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
    for k = 1:ns_bgs
        plot(ax{1}, x, loss_ortho(:,j,k),...
            s_bgs_lbl{k}, 'Color', skel_cmap(j,:));
        plot(ax{2}, x, res(:,j,k), ... 
            s_bgs_lbl{k}, 'Color', skel_cmap(j,:));
        plot(ax{3}, x, res_chol(:,j,k),...
            s_bgs_lbl{k}, 'Color', skel_cmap(j,:));
        if s_bgs(k) == 1
            lgd_str{end+1} = sprintf('%s',...
                skel_str{j}(2:end));
        else
            lgd_str{end+1} = sprintf('%s $\\circ$ HouseQR, s = %d',...
                skel_str{j}, s_bgs(k));
        end
    end
end
plot(ax{1}, x, eps*x, 'k--', x, eps*(x.^2), 'k-')

% Make plots pretty and save figures
folder_str = sprintf('results/%s_m%d_n%d', fstr, m, n);
mkdir(folder_str)
PrettyKappaPlot(fg, ax, lgd_str, folder_str);

% Save data
savestr = sprintf('%s/out', folder_str);
save(savestr, 'x', 'loss_ortho','res');

% close all;
end
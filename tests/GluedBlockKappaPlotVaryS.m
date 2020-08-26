function GluedBlockKappaPlotVaryS(XXdim, singval, skel, s_bgs)
% RUNTESTGLUEDVARYS(XXdim, singval, skel, s_bgs) is a wrapper function
% that executes a series of Block Gram-Schmidt (BGS) variants for a series
% of matrices specified by singval via the skeleton-muscle paradigm and
% returns loss of orthogonality and residual plots.  In more detail,
% RUNTESTGLUED 1. executes a loop that loads/generates the matrix or
% matrices specified
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
%   'BCGS' - BCGS with Pythagorean modification
%   'BCGS_CHOL' - BCGS with Cholesky modification
%
%   Default: skel = {'BCGS', 'BCGS_P', 'BCGS_CHOL'}
%
% Options for s_bgs:
%   s_bgs is a vector of block sizes for the call to BGS. This differs
%   from the block size used to generate the matrix XX.
%
%   Default: s_bgs = [1 2 5 10 20]
%
% (c) Kathryn Lund, Charles University, 2020

%%
addpath(genpath('../main/'))                                                % path to main routines
fstr = 'glued_vary_s';

% Defaults for inputs
if nargin == 0
    XXdim = [10000, 50, 10];
    singval = 1:8;
    skel = {'BCGS', 'BCGS_P', 'BCGS_CHOL'};
    skel_str = {'BCGS', 'BCGS-P', 'BCGS-chol'};
    s_bgs = [1 2 5 10 20];
elseif nargin == 1
    singval = 1:8;
    skel = {'BCGS', 'BCGS_P', 'BCGS_CHOL'};
    skel_str = {'BCGS', 'BCGS-P', 'BCGS-chol'};
    s_bgs = [1 2 5 10 20];
elseif nargin == 2
    skel = {'BCGS', 'BCGS_P', 'BCGS_CHOL'};
    skel_str = {'BCGS', 'BCGS-P', 'BCGS-chol'};
    s_bgs = [1 2 5 10 20];
elseif nargin == 3
    s_bgs = [1 2 5 10 20];
end

% Defaults for empty arguments
if isempty(XXdim)
    XXdim = [10000, 50, 10];
end
if isempty(singval)
    singval = 1:8;
end
if isempty(skel)
    skel = {'BCGS', 'BCGS_P', 'BCGS_CHOL'};
    skel_str = {'BCGS', 'BCGS-P', 'BCGS-chol'};
end
if isempty(s_bgs)
    s_bgs = [1 2 5 10 20];
end    

% Defaults for processing a single char array
if ischar(skel)
    skel = {skel};
end

% Default strings and replace underscore with tex underscore
if ~exist('skel_str','var')
    skel_str = skel;
    skel_str = strrep(skel_str, '_RO', '+');
    skel_str = strrep(skel_str, 'RO', '+');
    skel_str = strrep(skel_str, '_', '\_');
end

% Pre-allocate memory for measures
nmat = length(singval);
nskel = length(skel);
ns_bgs = length(s_bgs);
loss_ortho = zeros(nmat, nskel, ns_bgs);
res = zeros(nmat, nskel, ns_bgs);

% Extract dimensions
m = XXdim(1); p = XXdim(2); s = XXdim(3);
n = p*s;
I = eye(n);
XXnorm = zeros(1,nmat);
XXcond = zeros(1,nmat);

for i = 1:nmat
    % Create glued matrix
    matstr = sprintf('%s_cond%d_m%d_p%d_s%d.mat', fstr, singval(i), m, p, s);
    cd matrices
    if exist(matstr, 'file')
        load(matstr, 'XX')
    else
        XX = create_gluedmatrix(.5*singval(i), singval(i), m, p, s);
        
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
s_bgs_lbl = {'s-', 'o-', '*-', 'p-', 'h-', '.-', '^-'};
lgd_str = {};

x = XXcond; % condition number 
for j = 1:nskel
    for k = 1:ns_bgs
        plot(ax_loss_ortho, x, loss_ortho(:,j,k),...
            s_bgs_lbl{k}, 'Color', skel_cmap(j,:));
        plot(ax_res, x, res(:,j,k), ... 
            s_bgs_lbl{k}, 'Color', skel_cmap(j,:));
        if s_bgs(k) == 1
            lgd_str{end+1} = sprintf('%s',...
                skel_str{j}(2:end));
        else
            lgd_str{end+1} = sprintf('%s \\circ HouseQR, s = %d',...
                skel_str{j}, s_bgs(k));
        end
    end
end
plot(ax_loss_ortho, x, eps*x, 'k--', x, eps*(x.^2), 'k-')
plot(ax_res, x, eps*XXnorm, 'k--', x, eps*XXnorm.^2, 'k-')

titlestr = sprintf('s_{X} = %d', s);

set(ax_loss_ortho, 'Yscale', 'log', 'Xscale', 'log');
title(ax_loss_ortho, titlestr);
ylabel(ax_loss_ortho, 'loss of orthogonality ');
xlabel(ax_loss_ortho, '\kappa(X)')
lgd_str{end+1} = '\epsilon_M \kappa(X)';
lgd_str{end+1} = '\epsilon_M \kappa(X)^2';
legend(ax_loss_ortho, lgd_str, 'Location', 'BestOutside');

set(ax_res, 'Yscale', 'log', 'Xscale', 'log');
title(ax_res, titlestr);
ylabel(ax_res, 'relative residual');
xlabel(ax_res, '\kappa(X)')
lgd_str{end-1} = '\epsilon_M ||X||';
lgd_str{end} = '\epsilon_M ||X||^2';
legend(ax_res, lgd_str, 'Location', 'BestOutside');

% Save plots
folderstr = sprintf('results/%s_m%d_p%d_s%d', fstr, m, p, s);
mkdir(folderstr)

savestr = sprintf('%s/out', folderstr);
save(savestr,'loss_ortho','res');

savestr = sprintf('%s/loss_ortho', folderstr);
savefig(fig_loss_ortho, savestr, 'compact');
saveas(fig_loss_ortho, savestr, 'epsc')

savestr = sprintf('%s/res', folderstr);
savefig(fig_res, savestr, 'compact');
saveas(fig_res, savestr, 'epsc')

% close all;
end
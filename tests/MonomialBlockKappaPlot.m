function MonomialBlockKappaPlot(XXdim, svec, skel, musc)
% MONOMIALBLOCKKAPPAPLOT(XXdim, svec, skel, musc) compares loss of
% orthogonality and relative residual for different skeleton-muscle
% combinationss for a set of monomial matrices of size XXdim = [m p s] with
% varying block widths specified by the vector array svec.
%
% skel and musc should be given as either a char array or a cell of char
% arrays (i.e., text strings with single quotes).
%
% Options for XXdim and svec:
%   XXdim = [m p s], where m is the number of rows, p the number of block
%   vectors, and s the number of columns per block vector. svec is a vector
%   of scalars that divide n = p*s.  With s0 denoting one such scalar, a
%   matrix XX is produced with m rows and p = n/s0 block vectors each with
%   s0 columns. Each block vector X_k of XX is built from a random vector
%   v_k like
%
%       X_k = [v_k A*v_k ... A^(s0-1)*v_k].
%
%   Defaults:
%       XXdim = [1000 120 2];
%       svec = 2:2:12;
%
% Options for skel: see BGS.
%
% Options for musc: see INTRAORTHO.
%
% See also BLOCKKAPPAPLOT for more details about basic functionalities.
% Note that these tests are particularly slow, especially since large
% matrix sizes are needed to reveal interesting behavior.

%%
addpath(genpath('../main/'))                                                % path to main routines
addpath(genpath('auxiliary/'))
fstr = 'monomial_block_kappa_plot';

% Defaults for empty arguments
if isempty(XXdim)
    XXdim = [1000 120 2];
end
if isempty(svec)
    svec = 2:2:12;
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
nmat = length(svec);
nskel = length(skel);
nmusc = length(musc);
loss_ortho = zeros(nmat, nskel, nmusc);
res = zeros(nmat, nskel, nmusc);
res_chol = zeros(nmat, nskel, nmusc);
XXnorm = zeros(1,nmat);
XXcond = zeros(1,nmat);

% Extract dimensions
m = XXdim(1); p = XXdim(2); s = XXdim(3); n = p*s;
I = eye(n);

A = spdiags(linspace(.1,1,m)',0,m,m);
for i = 1:nmat
    % Create or load XX
    mat_s = svec(i); mat_p = n/mat_s;
    matstr = sprintf('monomial_m%d_p%d_s%d.mat', m, mat_p, mat_s);
    
    cd matrices
    if exist(matstr, 'file')
        load(matstr, 'XX');
    else
        XXhat = rand(m,mat_p);
        pp = 1:mat_p;
        XXhat(:,pp) = XXhat(:,pp)/norm(XXhat(:,pp));
        for k = 2:mat_s
            pp = pp + mat_p;
            XXhat(:,pp) = A*XXhat(:,pp - mat_p);
        end
        % Reshape XX
        XX = zeros(m,n);
        ind = 1:mat_s:n;
        kk = 1:mat_p;
        for k = 1:mat_s
            XX(:,ind) = XXhat(:,kk);
            ind = ind + 1;
            kk = kk + mat_p;
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
folder_str = sprintf('results/%s_m%d_p%d_s%d', fstr, m, mat_p, s);
mkdir(folder_str)
PrettyKappaPlot(fg, ax, lgd_str, folder_str);

% Save data
savestr = sprintf('%s/out', folder_str);
save(savestr, 'x', 'loss_ortho','res');

% close all;
end
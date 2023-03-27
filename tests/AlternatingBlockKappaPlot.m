function AlternatingBlockKappaPlot(XXdim, svec, skel, musc)
% ALTERNATINGBLOCKKAPPAPLOT(XXdim, svec, skel, musc) compares loss of
% orthogonality and relative residual for different skeleton-muscle
% combinationss for a set of so-called "alternating" matrices of size XXdim
% = [m p s] with varying block widths specified by the vector array svec.
%
% skel and musc should be given as either a char array or a cell of char
% arrays (i.e., text strings with single quotes).
%
% Options for XXdim and svec:
%   XXdim = [m p s], where m is the number of rows, p the number of block
%   vectors, and s the number of columns per block vector. svec is a vector
%   of scalars such that each entry of 2*vvec-1 divides n = p*s.  With s0
%   denoting one such scalar, a matrix XX is produced with m rows and p =
%   n/(2*s0-1) block vectors each with 2*s0-1 columns. Each block vector
%   X_k of XX is built from a random vector v_k like
%
%       X_k = [v_k [A^{-1}*v_k A*v_k] ... [A^{-s0+1}*v_k A^(s0-1)*v_k] ].
%
%   Defaults:
%       XXdim = [1000 315 2];
%       svec = 1:5;
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
fstr = 'alternating_block_kappa_plot';

% Defaults for empty arguments
if isempty(XXdim)
    XXdim = [1000 315 2];
end
if isempty(svec)
    svec = 1:5;
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

A = spdiags(linspace(1e-3,1e3,m)',0,m,m);
for i = 1:nmat
    % Create or load XX
    mat_s = 2*svec(i)-1; mat_p = n/mat_s;
    matstr = sprintf('alternating_m%d_p%d_s%d.mat', m, mat_p, mat_s);
    
    cd matrices
    if exist(matstr, 'file')
        load(matstr, 'XX');
    else
        % Simulate building extended Krylov via channels
        XXpower = rand(m, mat_p*svec(i));
        XXinv = rand(m, mat_p*(svec(i)-1));
        pp = 1:mat_p;
        XXpower(:,pp) = XXpower(:,pp)/norm(XXpower(:,pp));
        if svec(i) > 1
            XXinv(:,pp) = A\XXpower(:,pp);
            pp = pp + mat_p;
            XXpower(:,pp) = A*XXpower(:,pp);
            for k = 3:svec(i)
                XXinv(:,pp) = A\XXinv(:,pp - mat_p);
                pp = pp + mat_p;
                XXpower(:,pp) = A*XXpower(:,pp - mat_p);
            end
        end
        % Reshape XX
        XX = zeros(m,n);
        ind = 1:mat_s:n;
        kk = 1:mat_p;
        XX(:,ind) = XXpower(:,kk);
        for k = 2:svec(i)
            ind = ind + 1;
            XX(:,ind) = XXinv(:,kk);
            ind = ind + 1;
            kk = kk + mat_p;
            XX(:,ind) = XXpower(:,kk);
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
folder_str = sprintf('results/%s_m%d_p%d_s%d', fstr, m, p, s);
mkdir(folder_str)
PrettyKappaPlot(fg, ax, lgd_str, folder_str);

% Save data
savestr = sprintf('%s/out', folder_str);
save(savestr, 'x', 'loss_ortho','res');

% close all;
end
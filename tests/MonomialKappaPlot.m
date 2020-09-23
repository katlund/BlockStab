function MonomialKappaPlot(Xdim, svec, musc)
% MONOMIALKAPPAPLOT(Xdim, svec, musc) compares loss of
% orthogonality and relative residual for different skeleton-muscle
% combinationss for a set of monomial matrices of size Xdim = [m p s] with
% varying block widths specified by the vector array svec.
%
% skel and musc should be given as either a char array or a cell of char
% arrays (i.e., text strings with single quotes).
%
% Options for Xdim and svec:
%   Xdim = [m s], where m is the number of rows and s the number of
%   columns. svec is a vector of scalars that divide s.  With s0 denoting
%   one such scalar, a matrix X is produced with m rows and p = s/s0 block
%   vectors each with s0 columns.  Each block vector X_k of X is built from
%   a random vector v_k like
%
%       X_k = [v_k A*v_k ... A^(s0-1)*v_k].
%
%   Defaults:
%       Xdim = [1000 120];
%       svec = 2:2:12;
%
% Options for skel: see BGS.
%
% Options for musc: see INTRAORTHO.
%
% See also KAPPAPLOT for more details about basic functionalities. Note
% that these tests are particularly slow, especially since large matrix
% sizes are needed to reveal interesting behavior.

%%
addpath(genpath('../main/'))                                                % path to main routines
addpath(genpath('auxiliary/'))
fstr = 'monomial_kappa_plot';

% Defaults for empty arguments
if isempty(Xdim)
    Xdim = [1000 120];
end
if isempty(svec)
    svec = 2:2:12;
end

% Defaults for processing a single char array
if ischar(musc)
    musc = {musc};
end

% Default strings and format strings for plot legends
musc_str = AlgString(musc);

% Pre-allocate memory for measures
nmat = length(svec);
nmusc = length(musc);
loss_ortho = zeros(nmat, nmusc);
res = zeros(nmat, nmusc);
res_chol = zeros(nmat, nmusc);
Xnorm = zeros(1,nmat);
Xcond = zeros(1,nmat);

% Extract dimensions
m = Xdim(1); s = Xdim(2);
I = eye(s);

A = spdiags(linspace(.1,1,m)',0,m,m);
for i = 1:nmat
    % Create or load XX
    mat_s = svec(i); mat_p = s/mat_s;
    matstr = sprintf('monomial_m%d_p%d_s%d.mat', m, mat_p, mat_s);
    
    cd matrices
    if exist(matstr, 'file')
        load(matstr, 'XX');
    else
        Xhat = rand(m,mat_p);
        pp = 1:mat_p;
        Xhat(:,pp) = Xhat(:,pp)/norm(Xhat(:,pp));
        for k = 2:mat_s
            pp = pp + mat_p;
            Xhat(:,pp) = A*Xhat(:,pp - mat_p);
        end
        % Reshape X
        XX = zeros(m,s);
        ind = 1:mat_s:s;
        kk = 1:mat_p;
        for k = 1:mat_s
            XX(:,ind) = Xhat(:,kk);
            ind = ind + 1;
            kk = kk + mat_p;
        end
        
        save(matstr, 'XX');
    end
    cd ..
    Xnorm(i) = norm(XX);
    Xcond(i) = cond(XX);
    
    for k = 1:nmusc
        % Call BGS skeleton-muscle configuration
        [Q, R] = IntraOrtho(XX, musc{k});

        % Compute loss of orthonormality
        loss_ortho(i, k) = norm(I - Q' * Q, 2);

        % Compute relative residual
        res(i, k) = norm(XX - Q * R, 2) / Xnorm(i);

        % Compute relative residual for Cholesky relation
        res_chol(i, k) = norm(XX' * XX - R' * R, 2) / Xnorm(i)^2;

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
savestr = sprintf('%s/out', folder_str);
save(savestr, 'x', 'loss_ortho','res');

% close all;
end
 function KappaPlot(XXdim, logcondX, musc)
% KAPPAPLOT(XXdim, logcondX, musc) compares loss of orthogonality and
% relative residual for different muscles for a set of matrices of size
% XXdim = [m s] with varying singular values specified by the vector array
% logcondX.
%
% musc should be given as either a char array or a cell of char arrays
% (i.e., text strings with single quotes).
%
% Options for XXdim:
%   XXdim = [m s], where m is the number of rows and s is the number of
%   columns.
%
%   Default: XXdim = [1000 200]
%
% Options for logcondX:
%   logcondX should be a vector of negative powers.  Each entry of logcondX
%   is treated as the power of the smallest singular value for a matrix
%   whose singular values range from 10^0 to 10^t, where t is an entry of
%   logcondX.  For each t, a different matrix is generated, in order to see
%   how QR algorithms behave across a range of condition numbers.
%   
%   Default: logcondX = -(1:16)
%
% Options for musc: see INTRAORTHO

%%
addpath(genpath('../main/'))                                                % path to main routines
fstr = 'kappa_plot';

% Defaults for empty arguments
if isempty(XXdim)
    XXdim = [1000, 200];
end
if isempty(logcondX)
    logcondX = -(1:16);
end 

% Defaults for processing a single char array
if ischar(musc)
    musc = {musc};
end

% Default strings and replace underscore with tex underscore
musc_str = AlgString(musc);

% Pre-allocate memory for measures
nmat = length(logcondX);
nmusc = length(musc);
loss_ortho = zeros(nmat, nmusc);
res = zeros(nmat, nmusc);
res_chol = zeros(nmat, nmusc);
Xcond = zeros(1,nmat);
Xnorm = zeros(1,nmat);

% Extract dimensions
m = XXdim(1); s = XXdim(2);
I = eye(s);

% Plot settings
musc_cmap = lines(nmusc);
musc_lbl = {'s-', 'o-', '*-', '^-', 'p-', '.-', 'h-', 'd-'};

U = orth(randn(m,s));
V = orth(randn(s,s));
for i = 1:nmat
    Sigma = diag(logspace(0, logcondX(i), s)');
    X = U * Sigma * V';
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
        res_chol(i, k) = norm(X'*X - R'*R, 2) / Xnorm(i)^2;

        % Clear computed variables before next run
        clear QQ RR
    end
end

% Plot
x = Xcond;
lgd_str = musc_str;
fg = cell(1,3); ax = cell(1,3);
for i = 1:3
    fg{i} = figure;
    ax{i} = gca;
    hold on;
end
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
save(save_str,'loss_ortho', 'res', 'res_chol');

close all;
end
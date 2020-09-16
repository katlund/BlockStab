function GluedKappaPlot(Xdim, logcondX, musc)
% GLUEDKAPPAPLOT(Xdim, logcondX, musc) compares stability for different
% muscles for a set of glued matrices of size Xdim = [m s] with varying
% condition numbers loosely specified by the exponents in logcondX.
%
% musc should be given as either a char array or a cell of char arrays
% (i.e., text strings with single quotes).
%
% Options for Xdim:
%   Xdim = [m s], where m is the number of rows and s is the number of
%   columns.
%
%   Default: Xdim = [1000 200]
%
% Options for logcondXX:
%   logcondXX should be a vector of positive powers, where the power
%   denotes the log of the desired condition number of XX.
%   
%   Default: logcondXX = 1:8
%
% Options for musc: see INTRAORTHO.
%
% See also KAPPAPLOT for more details about basic functionalities.

%%
addpath(genpath('../main/'))                                                % path to main routines
addpath(genpath('auxiliary/'))
fstr = 'glued_kappa_plot';

% Defaults for inputs
if nargin == 0
    Xdim = [1000, 200];
    logcondX = 1:8;
elseif nargin == 1
    logcondX = 1:8;
end

% Defaults for empty arguments
if isempty(Xdim)
    Xdim = [1000, 200];
end
if isempty(logcondX)
    logcondX = 1:8;
end  

% Defaults for processing a single char array
if ischar(musc)
    musc = {musc};
end

% Default strings and format strings for plot legends
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
m = Xdim(1); s = Xdim(2);
I = eye(s);

% Fix glued dimensions
factors = factor(s);
mid_ind = round(length(s)/2);
r = prod(factors(1:mid_ind));
t = prod(factors(mid_ind+1:end));

for i = 1:nmat
    % Create glued matrix
    matstr = sprintf('glued_cond%d_m%d_p%d_s%d.mat', logcondX(i), m, r, t);
    cd matrices
    if exist(matstr, 'file')
        load(matstr, 'XX')
    else
        XX = CreateGluedMatrix(m, r, t,...
            .5*logcondX(i), logcondX(i));
        save(matstr, 'XX');
    end
    cd ..
    Xnorm(i) = norm(XX);
    Xcond(i) = cond(XX);
    
    for k = 1:nmusc
        % Call IntraOrtho
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
save_str = sprintf('%s/out', folder_str);
save(save_str, 'x', 'loss_ortho', 'res', 'res_chol');

% close all;
end
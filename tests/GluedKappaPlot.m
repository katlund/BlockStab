function GluedKappaPlot(XXdim, logcondXX, musc)
% RUNTESTGLUEDLOWSYNC(XXdim, glued, musc) is a wrapper function that
% compares stability for different muscles for a set of matrices of size
% XXdim = [m s] with varying singular values specified by the vector array
% glued.
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
% Options for logcondXX:
%   logcondXX should be a vector of positive powers, where the power
%   denotes the log of the desired condition number of XX.
%   
%   Default: logcondXX = 1:8
%
% Options for musc: see INTRAORTHO
%
%   Default: musc = {'CGS', 'MGS', 'MGS_SVL', 'MGS_LTS',...
%               'MGS_ICWY', 'MGS_CWY'};
%
% When run without arguments, RUNTESTGLUED returns loss of orthogonality
% and residual plots for default settings.
%
% (c) Kathryn Lund, Charles University, 2020

%%
addpath(genpath('../main/'))                                                % path to main routines
fstr = 'glued_low_sync';

% Defaults for inputs
if nargin == 0
    XXdim = [1000, 200];
    logcondXX = 1:8;
    musc = {'CGS', 'MGS', 'MGS_SVL', 'MGS_LTS',... 
        'MGS_ICWY', 'MGS_CWY'};
    musc_str = {'CGS', 'MGS', 'MGS\_SVL (MGS2)', 'MGS\_LTS (Alg 4)',...
        'MGS\_ICWY (Alg 5)', 'MGS\_CWY (Alg 6)'};
elseif nargin == 1
    logcondXX = 1:8;
    musc = {'CGS', 'MGS', 'MGS_SVL', 'MGS_LTS',... 
        'MGS_ICWY', 'MGS_CWY'};
    musc_str = {'CGS', 'MGS', 'MGS\_SVL (MGS2)', 'MGS\_LTS (Alg 4)',...
        'MGS\_ICWY (Alg 5)', 'MGS\_CWY (Alg 6)'};
elseif nargin == 2
    musc = {'CGS', 'MGS', 'MGS_SVL', 'MGS_LTS',... 
        'MGS_ICWY', 'MGS_CWY'};
    musc_str = {'CGS', 'MGS', 'MGS\_SVL (MGS2)', 'MGS\_LTS (Alg 4)',...
        'MGS\_ICWY (Alg 5)', 'MGS\_CWY (Alg 6)'};
end

% Defaults for empty arguments
if isempty(XXdim)
    XXdim = [1000, 200];
end
if isempty(logcondXX)
    logcondXX = 1:8;
end
if isempty(musc)
    musc = {'CGS', 'MGS', 'MGS_SVL', 'MGS_LTS',... 
        'MGS_ICWY', 'MGS_CWY'};
    musc_str = {'CGS', 'MGS', 'MGS\_SVL (MGS2)', 'MGS\_LTS (Alg 4)',...
        'MGS\_ICWY (Alg 5)', 'MGS\_CWY (Alg 6)'};
end    

% Defaults for processing a single char array
if ischar(musc)
    musc = {musc};
end

% Default strings and replace underscore with tex underscore
if ~exist('musc_str','var')
    musc_str = musc;
    musc_str = strrep(musc_str, '_RO', '+');
    musc_str = strrep(musc_str, 'RO', '+');
    musc_str = strrep(musc_str, '_', '\_');
end

% Pre-allocate memory for measures
nmat = length(logcondXX);
nmusc = length(musc);
loss_ortho = zeros(nmat, nmusc);

% Extract dimensions
m = XXdim(1); s = XXdim(2);
I = eye(s);
XXnorm = zeros(1,nmat);
XXcond = zeros(1,nmat);

% Plot settings
musc_cmap = lines(nmusc);
musc_lbl = {'s-', 'o-', '*-', '^-', 'p-', '.-', 'h-', 'd-'};

for i = 1:nmat
    % Create glued matrix
    pp = 20; ss = s/pp;
    matstr = sprintf('%s_cond%d_m%d_p%d_s%d.mat', fstr, logcondXX(i), m, pp, ss);
    cd matrices
    if exist(matstr, 'file')
        load(matstr, 'XX')
    else
        XX = create_gluedmatrix(.5*logcondXX(i), logcondXX(i), m, pp, ss);
        
        save(matstr, 'XX');
    end
    cd ..
    XXnorm(i) = norm(XX);
    XXcond(i) = cond(XX);
    
        for k = 1:nmusc
            % Call BGS skeleton-muscle configuration
            [QQ, ~] = IntraOrtho(XX, musc{k});

            % Compute loss of orthonormality
            loss_ortho(i, k) = norm(I - QQ'*QQ, 2);
            
            % Clear computed variables before next run
            clear QQ RR            
        end
end

% Plot
x = XXcond; % condition number 
lgd_str = musc_str;
fig_loss_ortho = clf; ax_loss_ortho = gca; hold on;
for k = 1:nmusc
    plot(ax_loss_ortho, x, loss_ortho(:,k),...
        musc_lbl{k}, 'Color', musc_cmap(k,:));
end
plot(ax_loss_ortho, x, eps*x, 'k--', x, eps*(x.^2), 'k-')
set(ax_loss_ortho, 'Yscale', 'log', 'Xscale', 'log');
title(ax_loss_ortho, 'Loss of Orthogonality ');
xlabel(ax_loss_ortho, '\kappa(X)')

lgd_str{end+1} = 'O(\epsilon) \kappa(X)';
lgd_str{end+1} = 'O(\epsilon) \kappa(X)^2';
legend(ax_loss_ortho, lgd_str, 'Location', 'BestOutside');

% Save plots
folderstr = sprintf('results/%s_m%d_s%d', fstr, m, s);
mkdir(folderstr)

savestr = sprintf('%s/out', folderstr);
save(savestr,'loss_ortho');

savestr = sprintf('%s/loss_ortho', folderstr);
savefig(fig_loss_ortho, savestr, 'compact');
saveas(fig_loss_ortho, savestr, 'epsc')

% close all;
end

%% Auxiliary functions
function XX = create_gluedmatrix(logcondXX_glob, logcondXX_block, m, p, s)
% Example 2 matrix from [Smoktunowicz et. al., 2006]
n = p*s;
XX = orth(rand(m,n));
XX = XX*diag(10.^(0:logcondXX_glob/(n-1):logcondXX_glob)) * orth(randn(n,n));
ibeg = 1;
iend = s;
for i = 1:p
    XX(:,ibeg:iend) = XX(:,ibeg:iend)*...
        diag(10.^(0:logcondXX_block/(s-1):logcondXX_block))*...
        orth(randn(s,s));
    ibeg = ibeg + s;
    iend = iend + s;
end
end
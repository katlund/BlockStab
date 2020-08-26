function RunTest(XXdim, mat, skel, musc, rpltol, verbose)
% RUNTEST(XXdim, mat, skel, musc, rpltol, verbose) is a wrapper function
% that executes panels of block Gram-Schmidt tests by doing the following:
%
% 1. Execute a loop that loads/generates the matrix or matrices specified
%    by mat for the dimensions [m p s] = XXdim.
% 2. Execute a second loop running through skeleton options specfied by
%    skel
% 3. Execute a third loop running through muscle options specified by musc
% 4. For each matrix, generate a heat-map chart for each measure-- loss of
%    orthogonality, residual, and timings
%
% All input variables should be given as either a char array or a cell of
% char arrays (i.e., text strings with single quotes).
%
% Options for mat:
%   'rand_uniform' - random entries drawn from uniform distribution
%   'rand_normal' - random entries drawn from normal distribution
%   'rank_def' - like rand_uniform but with a block vector set to 100 times
%       another, in order to force rank deficiency
%   'laeuchli' - the classic Läuchli matrix
%   'monomial' - a matrix resembling a sequence of block vectors
%       encountered in s-step Arnoldi with a monomial
%   'stewart' - a matrix with a geometric sequence of singular values,
%       ranging from 1 to 10^-t, t = 20
%   'stewart_extreme' - a matrix with a geometric sequence of singular
%       values, ranging from 1 to 10^-t, t = 10, for the first half and
%       exactly zeros for the second half
%   'hilbert' - HILB(m,p*s)
%
% Options for skel (see also BGS):
%   'BCGS' - Block Classical Gram-Schmidt
%   'BCGS_SROR' - Block Classical Gram-Schmidt with Selective
%       ReOrthonormalization and Replacement. Only compatible with
%       CGS_SROR.
%   'BCGS_IRO' - Block Classical Gram-Schmidt with Inner
%       ReOrthonormalization
%   'BCGS_IRO_T' - T-variant of BCGS_IRO
%   'BCGS_IRO_1' - reorthogonalize first step of BCGS_IRO
%   'BMGS' - Block Modified Gram-Schmidt
%   'BMGS_T' - T-variant of BMGS
%   'BMGS_SVL' - Block Modified Gram-Schmidt with Schreiber-Van-Loan
%       reformulation
%
% Options for musc (see also IntraOrtho):
%   'CGS' - Classical Gram-Schmidt.
%   'CGS_RO' - Classical Gram-Schmidt with ReOrthonormalization. Calls CGS
%       twice.
%   'CGS_IRO' - Classical Gram-Schmidt with Inner ReOrthonormalization.
%   'CGS_SRO' - Classical Gram-Schmidt with Selective ReOrthonormalization.
%       Calls CGS_SROR with rpltol = 0.
%   'CGS_SROR' - Classical Gram-Schmidt with Selective ReOrthonormalization
%       and Replacement
%   'MGS' - Modified Gram-Schmidt.
%   'MGS_RO' - Modified Gram-Schmidt with ReOrthonormalization. Cals MGS
%      twice.
%   'MGS_IRO' - Modified Gram-Schmidt with Inner ReOrthonormalization.
%   'MGS_SVL' - Modified Gram-Schmidt with Schreiber-Van-Loan reformulation
%   'HouseQR' - Householder-based QR. Calls built-in QR with economic
%       setting.
%   'CholQR' - Cholesky QR. Shifted CholQR is implemented when X'*X is not
%       positive definite.
%   'CholQR_RO' - Cholesky QR with ReOrthonormalization. Calls CHOLQR
%      twice.
%   'Sh_CholQR_RORO' - Shifted Cholesky QR with double ReOrthonormalization.
%
% When verbose is set to 1, a table for every matrix will print to screen.
% The default is 0.
%
% RUNTEST(XXdim, mat, skel, musc, rpltol) with any of the variables
% set to empty [] will invoke default settings for that variable.
%
% RUNTEST([], 'laeuchli', 'T', {'mgs', 'mgs_svl'}, []) will run
% the specified tests with T-variants of BCGS_IRO and BMGS, as well as
% BMGS_SVL.  Note that BMGS_SVL does not have a T-variant, because it
% already incorporates the "error sponge" matrix T.
%
% RUNTEST([], {'monomial', 'stewart_extreme'}, 'BCGS_IRO_1') will run the
% specified tests with the variant of BCGS_IRO that reorthogonalizes the
% first step.
%
% (c) Kathryn Lund, Charles University, 2020

%%
addpath(genpath('../main/'))                                                % path to main routines
addpath(genpath('matrices'))                                                % path to matrix files

% Defaults for inputs
if nargin == 0
    XXdim = [10000, 50, 10];
    mat = {'rand_uniform', 'rand_normal', 'rank_def',...
        'laeuchli', 'monomial', 'stewart', 'stewart_extreme', 'hilbert'};
    skel = {'BCGS', 'BCGS_IRO', 'BCGS_SROR', 'BMGS', 'BMGS_SVL'};
    skel_str = {'BCGS', 'BCGSI+', 'BCGSS+R', 'BMGS', 'BMGS\_SVL'};
    musc = {'CGS','CGS_RO', 'CGS_IRO' 'CGS_SRO', 'CGS_SROR', 'MGS', 'MGS_RO',...
        'MGS_IRO', 'MGS_SVL', 'HouseQR', 'CholQR', 'CholQR_RO', 'Sh_CholQR_RORO'};
    musc_str = {'CGS','CGS+', 'CGSI+', 'CGSS+', 'CGSS+R', 'MGS', 'MGS+',...
        'MGSI+','MGS\_SVL', 'HouseQR', 'CholQR', 'CholQR+', 'ShCholQR++'};
    rpltol = 100;
    verbose = 0;
elseif nargin == 1
    mat = {'rand_uniform', 'rand_normal', 'rank_def',...
        'laeuchli', 'monomial', 'stewart', 'stewart_extreme', 'hilbert'};
    skel = {'BCGS', 'BCGS_IRO', 'BCGS_SROR', 'BMGS', 'BMGS_SVL'};
    skel_str = {'BCGS', 'BCGSI+', 'BCGSS+R', 'BMGS', 'BMGS\_SVL'};
    musc = {'CGS','CGS_RO', 'CGS_IRO' 'CGS_SRO', 'CGS_SROR', 'MGS', 'MGS_RO',...
        'MGS_IRO', 'MGS_SVL', 'HouseQR', 'CholQR', 'CholQR_RO', 'Sh_CholQR_RORO'};
    musc_str = {'CGS','CGS+', 'CGSI+', 'CGSS+', 'CGSS+R', 'MGS', 'MGS+',...
        'MGSI+','MGS\_SVL', 'HouseQR', 'CholQR', 'CholQR+', 'ShCholQR++'};
    rpltol = 100;
    verbose = 0;
elseif nargin == 2
    skel = {'BCGS', 'BCGS_IRO', 'BCGS_SROR', 'BMGS', 'BMGS_SVL'};
    skel_str = {'BCGS', 'BCGSI+', 'BCGSS+R', 'BMGS', 'BMGS\_SVL'};
    musc = {'CGS','CGS_RO', 'CGS_IRO' 'CGS_SRO', 'CGS_SROR', 'MGS', 'MGS_RO',...
        'MGS_IRO', 'MGS_SVL', 'HouseQR', 'CholQR', 'CholQR_RO', 'Sh_CholQR_RORO'};
    musc_str = {'CGS','CGS+', 'CGSI+', 'CGSS+', 'CGSS+R', 'MGS', 'MGS+',...
        'MGSI+','MGS\_SVL', 'HouseQR', 'CholQR', 'CholQR+', 'ShCholQR++'};
    rpltol = 100;
    verbose = 0;
elseif nargin == 3
    musc = {'CGS','CGS_RO', 'CGS_IRO' 'CGS_SRO', 'CGS_SROR', 'MGS', 'MGS_RO',...
        'MGS_IRO', 'MGS_SVL', 'HouseQR', 'CholQR', 'CholQR_RO', 'Sh_CholQR_RORO'};
    musc_str = {'CGS','CGS+', 'CGSI+', 'CGSS+', 'CGSS+R', 'MGS', 'MGS+',...
        'MGSI+','MGS\_SVL', 'HouseQR', 'CholQR', 'CholQR+', 'ShCholQR++'};
    rpltol = 100;
    verbose = 0;
elseif nargin == 4
    rpltol = 100;
    verbose = 0;
elseif nargin == 5
    verbose = 0;
end

if isempty(XXdim)
    XXdim = [10000, 50, 10];
end
if isempty(mat)
    mat = {'rand_uniform', 'rand_normal', 'rank_def',...
        'laeuchli', 'monomial', 'stewart', 'stewart_extreme', 'hilbert'};
end
if isempty(skel)
    skel = {'BCGS', 'BCGS_IRO', 'BCGS_SROR', 'BMGS', 'BMGS_SVL'};
    skel_str = {'BCGS', 'BCGSI+', 'BCGSS+R', 'BMGS', 'BMGS\_SVL'};
end
if strcmp(skel, 'T')
    skel = {'BCGS_IRO_T', 'BMGS_T', 'BMGS_SVL'};
    skel_str = {'BCGSI+T', 'BMGST', 'BMGS\_SVL'};
end
if strcmp(skel, 'BCGS_IRO_1')
    skel_str = {'BCGSI+1'};
end
if isempty(musc)
    musc = {'CGS','CGS_RO', 'CGS_IRO' 'CGS_SRO', 'CGS_SROR', 'MGS', 'MGS_RO',...
        'MGS_IRO', 'MGS_SVL', 'HouseQR', 'CholQR', 'CholQR_RO', 'Sh_CholQR_RORO'};
    musc_str = {'CGS','CGS+', 'CGSI+', 'CGSS+', 'CGSS+R', 'MGS', 'MGS+',...
        'MGSI+','MGS\_SVL', 'HouseQR', 'CholQR', 'CholQR+', 'ShCholQR++'};
end
if isempty(rpltol)
    rpltol = 100;
end    

% Defaults for processing a single char array
if ischar(mat)
    mat = {mat};
end
if ischar(skel)
    skel = {skel};
end
if ischar(musc)
    musc = {musc};
end

% Default strings and replace underscore with tex underscore
if ~exist('skel_str','var')
    skel_str = skel;
    skel_str = strrep(skel_str, '_RO', '+');
    skel_str = strrep(skel_str, 'RO', '+');
    skel_str = strrep(skel_str, '_', '\_');
end
if ~exist('musc_str','var')
    musc_str = musc;
    musc_str = strrep(musc_str, '_RO', '+');
    musc_str = strrep(musc_str, 'RO', '+');
    musc_str = strrep(musc_str, '_', '\_');
end

nmat = length(mat);
nskel = length(skel);
nmusc = length(musc);

% Extract dimensions
m = XXdim(1); p = XXdim(2); s = XXdim(3);
n = p*s;
I = eye(n);

for i = 1:nmat
    matstr = sprintf('%s_m%d_p%d_s%d', mat{i}, m, p, s);
    checkstr = sprintf('%s.mat',matstr);
    cd matrices
    if exist(checkstr, 'file') == 2
        load(checkstr, 'XX', 'XXstr', 'XXprops');
    else
        [XX, XXstr, XXprops] = MatGen(mat{i}, XXdim);
    end
    cd ..
    fprintf('%s\n',XXstr)
    XXnorm = norm(XX, 2);
    display(XXprops.cond)
    
    % Pre-allocate memory for measures
    loss_ortho = zeros(nmusc, nskel);
    res = zeros(nmusc, nskel);
    run_time = zeros(nmusc, nskel);
    
    for j = 1:nskel
        for k = 1:nmusc
            if strcmpi(skel{j}, 'bcgs_sror')
                if strcmpi(musc{k}, 'cgs_sror')
                    [QQ, RR, ~, tt] = BGS(XX, s, 'bcgs_sror', 'cgs_sror', rpltol);
                elseif strcmpi(musc{k}, 'cgs_sro')
                    [QQ, RR, ~, tt] = BGS(XX, s, 'bcgs_sror', 'cgs_sro', 0);
                else
                    QQ = NaN;
                    RR = NaN;
                    tt = NaN;
                end
            else
                % Call BGS skeleton-muscle configuration
                [QQ, RR, ~, tt] = BGS(XX, s, skel{j}, musc{k}, rpltol);
            end
            
            % Compute loss of orthonormality
            loss_ortho(k, j) = norm(I - QQ'*QQ, 2);

            % Compute relative residual
            res(k, j) = norm(XX - QQ*RR, 2)/XXnorm;

            % Store runtime
            run_time(k, j) = tt;

            % Clear computed variables before next run
            clear QQ RR tt            
        end
    end
    
    % Print tables to screen
    if verbose
        fprintf('||I - Q''*Q||\n')
        display(array2table(loss_ortho, 'VariableNames', skel, 'RowNames', musc))

        fprintf('||X - Q*R|| / ||X||\n')
        display(array2table(res, 'VariableNames', skel, 'RowNames', musc))

%         fprintf('Runtimes\n')
%         display(array2table(run_time, 'VariableNames', skel, 'RowNames', musc))
    end
    
    %% Generate and save heatmaps as .eps files (for TeX use)
    mkdir results
    cd results
    mkdir(matstr);
    cd ..
    
    hax = cell(3,1);
    hfg = cell(3,1);
    hfg{1} = figure('visible', 'off');
    hax{1} = heatmap(loss_ortho);
    hax{1}.Colormap = parula*.9;
    hax{1}.Title = 'Loss of Orthogonality, log10-scale';
    
    hfg{2} = figure('visible', 'off');
    hax{2} = heatmap(res);
    hax{2}.Colormap = spring*.9;
    hax{2}.Title = 'Relative Residual, log10-scale';
    
    hfg{3} = figure('visible', 'off');
    hax{3} = heatmap(run_time);
    hax{3}.Colormap = autumn*.9;
    hax{3}.Title = 'Runtime, seconds';
    
    % Add labels and change norm data to log scale
    measstr = {'loss_ortho', 'res', 'run_time'};
    for j = 1:3
        set(hax{j}, 'XDisplayLabels', skel_str, 'YDisplayLabels', musc_str)
        if j ~= 3
            cdj = get(hax{j}, 'ColorData');
            set(hax{j},'ColorData', log10(cdj), 'ColorLimits', [-16 0]);
        end
        set(hax{j},'CellLabelFormat', '%0.2f',...
            'MissingDataColor', .75*ones(1,3),...
            'FontSize', 14);
        
        savestr = sprintf('results/%s/out',matstr);
        save(savestr,'loss_ortho','res','run_time', 'XXprops');        
        savestr2 = sprintf('results/%s/%s', matstr, measstr{j});
        savefig(hfg{j}, savestr2, 'compact');
        saveas(hfg{j}, savestr2, 'epsc')
    end
    close all;  % frees up memory; otherwise Matlab keeps all figures open in the background
end
fprintf('To open figures, navigate to the appropriate folder and use OPENFIG with ''visible'' option\n');

% matProps(XXdim);    % Generate table for paper
end
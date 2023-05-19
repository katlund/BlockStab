function MakeHeatmap(XXdim, mat, skel, musc, rpltol, verbose)
% MAKEHEATMAP(XXdim, mat, skel, musc, rpltol, verbose) is a wrapper
% function that executes panels of block Gram-Schmidt tests by doing the
% following:
%
% 1. Execute a loop that loads/generates the matrix or matrices specified
%    by mat for the dimensions [m p s] = XXdim.
% 2. Execute a second loop running through skeleton options specfied by
%    skel
% 3. Execute a third loop running through muscle options specified by musc
% 4. For each matrix, generate a heatmap for each measure-- loss of
%    orthogonality, relative residual, and timings
%
% All input variables should be given as either a char array or a cell of
% char arrays (i.e., text strings with single quotes).
%
% XXdim = [m p s]:
%   m - number of rows
%   p - number of block vectors
%   s - number of columns per block vector
%
%   Default is [1000 100 5].
%
% Options for mat (see also MATGEN):
%   'rand_uniform' - random entries drawn from uniform distribution
%   'rand_normal' - random entries drawn from normal distribution
%   'rank_def' - like rand_uniform but with a block vector set to 100 times
%       another, in order to force rank deficiency
%   'laeuchli' - the classic L�uchli matrix
%   'monomial' - a matrix resembling a sequence of block vectors
%       encountered in s-step Arnoldi with a monomial
%   'stewart' - a matrix with a geometric sequence of singular values,
%       ranging from 1 to 10^-t, t = 20
%   'stewart_extreme' - a matrix with a geometric sequence of singular
%       values, ranging from 1 to 10^-t, t = 10, for the first half and
%       exactly zeros for the second half
%   'hilbert' - the Hilbert matrix, generated by the built-in function
%      HILB(m,p*s)
%   's-step' - a matrix with columns spanning a Krylov subspace of size sp
%   'newton' - like 's-step', but better conditioned
%
%   Default is 'laeuchli'.
%
% Options for skel (see also BGS):
%   'BCGS' - Block Classical Gram-Schmidt (1-sync/block)
%   'BCGS_PIP' - BCGS with Pythagorean Inner Product reformulation
%       (1-sync/block)
%   'BCGS_PIO' - BCGS with Pythagorean IntraOrtho reformulation
%       (1-sync/block)
%   'BCGS_SROR' - BCGS with Selective ReOrthonormalization and Replacement.
%       Only compatible with CGS_SROR.
%   'BCGS_IRO' - BCGS with Inner ReOrthonormalization (2-sync/block)
%   'BCGS_IRO_T' - T-variant of BCGS_IRO (2-sync/block)
%   'BCGS_IRO_1' - reorthogonalize first step of BCGS_IRO (2-sync/block)
%   'BCGS_IRO_LS' - Low-sync version of BCGS_IRO (1-sync/block)
%   'BMGS' - Block Modified Gram-Schmidt
%   'BMGS_T' - T-variant of BMGS
%   'BMGS_SVL' - BMGS with Schreiber-Van-Loan reformulation (2-sync/block)
%   'BMGS_LTS' - BMGS with Lower Triangular Solve (2-sync/block)
%   'BMGS_CWY' - BMGS with Compact WY reformulation (1-sync/block)
%   'BMGS_ICWY' - BMGS with Inverse Compact WY reformulation (1-sync/block)
%
%   Default is {'BCGS', 'BCGS_IRO', 'BCGS_SROR', 'BCGS_IRO_LS', 'BMGS',
%   'BMGS_SVL', 'BMGS_CWY'}.
%
% Options for musc (see also IntraOrtho):
%   'CGS' - Classical Gram-Schmidt (1-sync/col)
%   'CGS_RO' - CGS with ReOrthonormalization. Calls CGS twice. (2-sync/col)
%   'CGS_IRO' - CGS with Inner ReOrthonormalization (2-sync/col)
%   'CGS_SRO' - CGS with Selective ReOrthonormalization. Calls CGS_SROR
%       with rpltol = 0.
%   'CGS_SROR' - CGS with Selective ReOrthonormalization and Replacement
%   'CGS_IRO_LS' - Low-sync version of CGS_IRO (1-sync/col)
%   'MGS' - Modified Gram-Schmidt.
%   'MGS_RO' - MGS with ReOrthonormalization. Calls MGS twice.
%   'MGS_IRO' - MGS with Inner ReOrthonormalization.
%   'MGS_SVL' - MGS with Schreiber-Van-Loan reformulation (2-sync/col)
%   'MGS_LTS' - MGS with Lower Triangular Solve (2-sync/col)
%   'MGS_CWY' - MGS with Compact WY reformulation (1-sync/col)
%   'MGS_SVL' - MGS with Inverse Compact WY reformulation (1-sync/col)
%   'HouseQR' - Householder-based QR. Calls built-in QR with economic
%       setting.
%   'CholQR' - Cholesky QR. Shifted CholQR is implemented when X'*X is not
%       positive definite. (1-sync)
%   'CholQR_RO' - Cholesky QR with ReOrthonormalization. Calls CHOLQR
%      twice. (2-sync)
%   'Sh_CholQR_RORO' - Shifted Cholesky QR with double ReOrthonormalization
%      (3-sync)
%   
%   Default is {'CGS', 'CGS_IRO', 'CGS_SRO', 'CGS_SROR', 'CGS_IRO_LS',
%   'MGS', 'MGS_SVL', 'MGS_CWY', 'HouseQR', 'CholQR', 'CholQR_RO',
%   'Sh_CholQR_RORO'};
%
% When verbose is set to 1, a table for every matrix will print to screen.
% The default is 0.
%
% MAKEHEATMAP(XXdim, mat, skel, musc, rpltol) with any of the variables
% set to empty [] will invoke default settings for that variable.
%
% -------------------------------------------------------------------------
% Examples:
%
% MAKEHEATMAP([], [], [], []) will generate a heatmap for the laeuchli
% matrix and a wide sampling of skeletons and muscles.
%
% MAKEHEATMAP([], 'monomial', {'BCGS', 'BCGS_IRO'}, {'HouseQR', 'CGS_RO'})
% will generate a heatmap for the monomial matrix and the four specified
% skeleton-muscle combinations.
%
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
addpath(genpath('../main/'))                                                % path to main routines
addpath(genpath('matrices/'))                                               % path to matrix files
addpath(genpath('auxiliary/'))                                              % path to auxiliary files

% Defaults for inputs
if nargin == 4
    rpltol = 100;
    verbose = 0;
elseif nargin == 5
    verbose = 0;
end

if isempty(XXdim)
    XXdim = [1000 100 5];
end
if isempty(mat)
    mat = {'laeuchli'};
end
if isempty(skel)
    skel = {'BCGS', 'BCGS_IRO', 'BCGS_SROR', 'BCGS_IRO_LS', 'BMGS', 'BMGS_SVL', 'BMGS_CWY'};
end
if isempty(musc)
    musc = {'CGS','CGS_RO', 'CGS_IRO', 'CGS_SRO', 'CGS_SROR', 'CGS_IRO_LS',...
        'MGS', 'MGS_SVL', 'MGS_CWY', 'HouseQR', 'CholQR', 'CholQR_RO', 'Sh_CholQR_RORO'};
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

% Default strings and format strings for plot legends
skel_str = basic_strrep(skel);
musc_str = basic_strrep(musc);

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
        [XX, XXstr, XXprops] = mat_gen(mat{i}, XXdim);
    end
    cd ..
    fprintf('%s, cond(XX) = %2.2e\n', XXstr, XXprops.cond)
    XXnorm = norm(XX, 2);
    
    % Pre-allocate memory for measures
    loss_ortho = zeros(nmusc, nskel);
    res = zeros(nmusc, nskel);
    run_time = zeros(nmusc, nskel);
    
    for j = 1:nskel
        for k = 1:nmusc
            if strcmpi(skel{j}, 'bcgs_sror')
                param = [];
                if strcmpi(musc{k}, 'cgs_sror')
                    param.rpltol = rpltol;
                    [QQ, RR, ~, tt] = BGS(XX, s, 'bcgs_sror', 'cgs_sror', param);
                elseif strcmpi(musc{k}, 'cgs_sro')
                    param.rpltol = 0;
                    [QQ, RR, ~, tt] = BGS(XX, s, 'bcgs_sror', 'cgs_sro', param);
                else
                    QQ = NaN(m, n);
                    RR = NaN(n, n);
                    tt = NaN;
                end
            else
                % Call BGS skeleton-muscle configuration
                [QQ, RR, ~, tt] = BGS(XX, s, skel{j}, musc{k});
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

        fprintf('Runtimes\n')
        display(array2table(run_time, 'VariableNames', skel, 'RowNames', musc))
    end
    
    %% Generate and save heatmaps as .eps files (for TeX use)
    mkdir results\heatmap
    cd results\heatmap
    mkdir(matstr);
    cd ..\..
    
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
        
        savestr = sprintf('results/heatmap/%s/out', matstr);
        save(savestr,'loss_ortho','res','run_time', 'XXprops');        
        savestr2 = sprintf('results/heatmap/%s/%s', matstr, measstr{j});
        savefig(hfg{j}, savestr2, 'compact');
        saveas(hfg{j}, savestr2, 'epsc')
    end
    close all;  % frees up memory; otherwise Matlab keeps all figures open in the background
end
fprintf('To open figures, navigate to ''results/heatmap'' and the desired\n');
fprintf('subfolder therein and use OPENFIG with ''visible'' option.\n');

% matProps(XXdim);    % Generate table for paper
end
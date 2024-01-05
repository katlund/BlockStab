function tex_report(run_data)
% TEX_REPORT(run_data) generates a LaTeX report from RUNKAPPAPLOT requiring
% only basic packages.
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
% Open new file
save_str = sprintf('%s/report.tex', run_data.dir_str);
fID = fopen(save_str,'w');

% Preamble
fprintf(fID, '\\documentclass[10pt]{article}\n');
fprintf(fID, '\\usepackage{booktabs, makecell} %% for column headers\n');
fprintf(fID, '\\usepackage{datetime2}\n');
fprintf(fID, '\\usepackage{enumerate}\n');
if run_data.options.save_eps
    fprintf(fID, '\\usepackage{epsfig, epstopdf}\n');
end
fprintf(fID, '\\usepackage{float}\n');
fprintf(fID, '\\usepackage{geometry, graphicx}\n');
fprintf(fID, '\\geometry{a4paper}\n');
fprintf(fID, '\n');
fprintf(fID, '\\title{Report}\n');
fprintf(fID, '\\author{}\n');
fprintf(fID, '\\date{%s}\n',char(run_data.dtnow,'dd-MMM-yyyy HH:mm:ss'));
fprintf(fID, '\n');

% Begin document
fprintf(fID, '\\begin{document}\n');
fprintf(fID, '\\maketitle\n');
fprintf(fID, '\n');

% Run details
fprintf(fID, '\\section*{Test Details}\n');
fprintf(fID, '\\begin{itemize}\n');
fprintf(fID, '\t\\item Matrix type: \\texttt{%s}\n', ...
    run_data.mat_type);
fprintf(fID, '\t\\item Dimensions of $X$: $m = %d$, $p = %d$, $s = %d$\n', ...
    run_data.options.num_rows, ...
    run_data.options.num_partitions, ...
    run_data.options.block_size);
fprintf(fID, '\\end{itemize}\n');
fprintf(fID, '\n\n');

% Print algorithm configurations
n_alg = length(run_data.skel);
fprintf(fID, '\\section*{Algorithm Configurations}\n');
fprintf(fID, '\\begin{enumerate}[(1)]\n');
for i = 1:n_alg
    % Print algorithm
    lgd_str = split(run_data.lgd{i},') ');
    fprintf(fID, '\t\\item %s\n', lgd_str{2});

    % Extra temporary parameters to simplify logic
    musc = run_data.musc{i};
    skel = run_data.skel{i};

    % Set flags
    if isempty(musc)
        chol_flag = false;
    else
        chol_flag = contains(lower(musc), 'chol');
    end
    if isempty(skel)
        mp_flag = false;
    else
        chol_flag = chol_flag || ...
        strcmpi(skel, 'bcgs_iro_ls') || ...
        strcmpi(skel, 'bcgs_iro_2s') || ...
        strcmpi(skel, 'bcgs_iro_1s') || ...
        contains(lower(skel), 'bcgs_pi') || ...
        contains(lower(skel), 'cwy');
        mp_flag = contains(lower(skel), '_mp');
    end

    % Fill in missing parameters to simplify logic
    param = run_data.param{i};
    param = param_init(param);
    if mp_flag
        param = mp_param_init(param);
    end

    if isempty(skel)
    %% Muscle only
        % Print INTRAORTHO defaults only for the muscles that take more
        % than 'verbose' as a parameters.  Note that chol and rpltol are
        % mutually exclusive.

        % CGS_SRO(R)
        if strcmpi(musc, 'cgs_sro')
            fprintf(fID, '\t\\begin{itemize}\n');
            fprintf(fID, '\t\t\\item Replacement tolerance: %2.2e\n', 0);
            fprintf(fID, '\t\\end{itemize}\n');
        end
        if strcmpi(musc, 'cgs_sror')
            fprintf(fID, '\t\\begin{itemize}\n');
            fprintf(fID, '\t\t\\item Replacement tolerance: %2.2e\n', ...
                param.rpltol);
            fprintf(fID, '\t\\end{itemize}\n');
        end

        % Cholesky
        if chol_flag
            fprintf(fID, '\t\\begin{itemize}\n');
            switch param.chol
                case 'chol'
                    fprintf(fID, '\t\t\\item Cholesky: built-in \\texttt{chol}\n');
                case 'chol_nan'
                    fprintf(fID, '\t\t\\item Cholesky: \\texttt{chol\\_nan}\n');
                case 'chol_free'
                    fprintf(fID, '\t\t\\item Cholesky: \\texttt{chol\\_free}\n');
            end
            fprintf(fID, '\t\\end{itemize}\n');
        end

    else
    %% Skeleton-muscle
        % BCGSSR+R first (does not take Chol or MP options)
        if strcmpi(skel, 'bcgs_sror')
            fprintf(fID, '\t\\begin{itemize}\n');
            fprintf(fID, '\t\t\\item Replacement tolerance: %2.2e\n',...
                param.rpltol);
            fprintf(fID, '\t\\end{itemize}\n');
            break
        end

        if mp_flag || chol_flag
            fprintf(fID, '\t\\begin{itemize}\n');
        end

        % Cholesky
        if chol_flag
            switch param.chol
                case 'chol'
                    fprintf(fID, '\t\t\\item Cholesky: built-in \\texttt{chol}\n');
                case 'chol_nan'
                    fprintf(fID, '\t\t\\item Cholesky: \\texttt{chol\\_nan}\n');
                case 'chol_free'
                    fprintf(fID, '\t\t\\item Cholesky: \\texttt{chol\\_free}\n');
            end
        end

        % MP
        if mp_flag
            switch param.mp_package
                case 'advanpix'
                    fprintf(fID, '\t\t\\item MP Package: Advanpix\n');
                    fprintf(fID, '\t\t\\item MP Digits: %d\n', param.mp_digits);

                case {'symbolic math', 'symbolic toolbox', 'vpa'}
                    fprintf(fID, '\t\t\\item MP Package: Symbolic Math\n');
                    fprintf(fID, '\t\t\\item MP Digits: %d\n', param.mp_digits);

                otherwise
                    % Do nothing
            end
        end

        if mp_flag || chol_flag
            fprintf(fID, '\t\\end{itemize}\n');
        end
    end
end
fprintf(fID, '\\end{enumerate}\n');
fprintf(fID, '\n\n');

% Kappa Plots
fprintf(fID, '\\section*{$\\kappa$ Plots}\n');
plot_str = {'loss_ortho', 'rel_res', 'rel_chol_res'};
for i = 1:3
    fprintf(fID, '\\begin{figure}[H]\n');
    fprintf(fID, '\t\\begin{center}\n');
    if run_data.options.save_pdf
        fprintf(fID, '\t\t\\resizebox{.95\\textwidth}{!}{\\includegraphics{%s.pdf}}\n', plot_str{i});
    elseif run_data.options.save_eps
        fprintf(fID, '\t\t\\resizebox{.95\\textwidth}{!}{\\includegraphics{%s.eps}}\n', plot_str{i});
    end
    fprintf(fID, '\t\\end{center}\n');
    fprintf(fID, '\\end{figure}\n');
    fprintf(fID, '\n');
end

% End document
fprintf(fID, '\\end{document}\n');

fclose(fID);
end
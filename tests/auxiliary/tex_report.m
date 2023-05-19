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
fprintf(fID, '\\usepackage{epsfig, epstopdf}\n');
fprintf(fID, '\\usepackage{float}\n');
fprintf(fID, '\\usepackage{geometry, graphicx}\n');
fprintf(fID, '\\geometry{a4paper}\n');
fprintf(fID, '\n');
fprintf(fID, '\\title{$\\kappa$ Plots}\n');
fprintf(fID, '\\author{}\n');
fprintf(fID, '\\date{%s}\n',char(run_data.datetime,'dd-MMM-yyyy HH:mm:ss'));
fprintf(fID, '\n');

% Begin document
fprintf(fID, '\\begin{document}\n');
fprintf(fID, '\\maketitle\n');
fprintf(fID, '\n');

% Run details
fprintf(fID, '\\section{Test Details}\n');
fprintf(fID, '\\begin{itemize}\n');
fprintf(fID, '\t\\item Matrix type: \\texttt{%s}\n', ...
    run_data.options.mat_type);
fprintf(fID, '\t\\item Dimensions of $X$: $n = %d$, $p = %d$, $s = %d$\n', ...
    run_data.options.num_rows, ...
    run_data.options.num_partitions, ...
    run_data.options.block_size);
fprintf(fID, '\\end{itemize}\n');

% Print algorithm configurations
n_alg = length(run_data.skel);
for i = 1:n_alg
    lgd_str = split(run_data.lgd{i},') ');
    fprintf(fID, '\\begin{enumerate}[(1)]\n');
    fprintf(fID, '\t\\item %s', lgd_str{2});
    fprintf(fID, '\\end{enumerate}\n');
end

% Kappa Plots
plot_str = {'loss_ortho', 'rel_res', 'rel_chol_res'};
for i = 1:3
    save_str = sprintf('%s/%s', run_data.dir_str, plot_str{i});
    fprintf(fID, '\\begin{figure}\n');
    fprintf(fID, '\t\\begin{center}\n');
    fprintf(fID, '\t\t\\resizebox{.9\\textwidth}{!}{\\includegraphics{%s.eps}} \\\\\n', save_str);
    fprintf(fID, '\t\\end{tabular}\n');
    fprintf(fID, '\t\\end{center}\n');
    fprintf(fID, '\\end{figure}\n');
    fprintf(fID, '\n');
end

% End document
fprintf(fID, '\\end{document}\n');

fclose(fID);
end
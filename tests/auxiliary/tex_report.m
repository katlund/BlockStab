function tex_report(run_options, dir_str, skel, musc, param)
%


% Modified from https://gitlab.mpi-magdeburg.mpg.de/lund/low-sync-block-arnoldi/-/blob/main/main/tex_bfom.m
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
% Set-up
plot_str = {'loss_ortho', 'rel_res', 'rel_chol_res'};


% Open new file
save_str = sprintf('%s/figures.tex', run_data.save_str);
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
fprintf(fID, '\\section{Run Details}\n');
fprintf(fID, '\\begin{itemize}\n');
fprintf(fID, '\t\\item Matrix type: \\texttt{%s}\n', replace(options.matrix_type, '_', '\_') );
fprintf(fID, '\t\\item Dimensions of $X$: $%d \\times %d$\n', run_data.prob_size);
fprintf(fID, '\t\\item Block size: $%d$\n', run_data.block_size);
fprintf(fID, '\\end{itemize}\n');

% Kappa Plots
for i = 1
    for a = 1:num_alg
        if mod(a,2)
            if a ~= num_alg
                fprintf(fID, '\\begin{figure}[H]\n');
                fprintf(fID, '\t\\begin{tabular}{cc}\n');
                fprintf(fID, '\t\t\\resizebox{.45\\textwidth}{!}{\\includegraphics{%s_%s.eps}} &\n', alg_config{a}, conv_str{i});
            else % then a is the last one
                fprintf(fID, '\\begin{figure}[H]\n');
                fprintf(fID, '\t\\resizebox{.45\\textwidth}{!}{\\includegraphics{%s_%s.eps}}\n', alg_config{a}, conv_str{i});
                fprintf(fID, '\\end{figure}\n');
                fprintf(fID, '\n');
            end
        else
            fprintf(fID, '\t\t\\resizebox{.45\\textwidth}{!}{\\includegraphics{%s_%s.eps}} \\\\\n', alg_config{a}, conv_str{i});
            fprintf(fID, '\t\\end{tabular}\n');
            fprintf(fID, '\\end{figure}\n');
            fprintf(fID, '\n');
        end
    end
end

% End document
fprintf(fID, '\\end{document}\n');

fclose(fID);
end
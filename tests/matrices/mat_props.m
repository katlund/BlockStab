function matProps(XXdim)
% MATPROPS(XXdim) generates a LaTeX table containing properies of the
% matrices used in RUNTEST for dimensions XXdim = [m, p, s].

%%

m = XXdim(1); p = XXdim(2); s = XXdim(3);
sizestr = sprintf('m%d_p%d_s%d', m, p, s);
mat = {'rand_uniform', 'rand_normal', 'rank_def','laeuchli', 'monomial', 'stewart', 'stewart_extreme', 'hilbert', 's-step', 'newton'};
matstr = {'\randuniform', '\randnormal', '\rankdef','\laeuchli', '\monomial', '\stewart', '\stewartextreme', '\hilbert', '\sstep', '\newton'};
loadstr = mat;
nummat = length(mat);
for i = 1:nummat
    loadstr{i} = sprintf('%s_%s', mat{i}, sizestr);
end

fID = fopen('matrix_props','w');

fprintf(fID, '\\begin{tabular}{c|c|c|c|c}\n');
fprintf(fID, '    Matrix ID	& $\\sigma_{1}$	& $\\sigma_{n}$	& $\\kappa(\\XX)$ & $\\normF{\\XX}$ \\\\ \\hline\n');
for i = 1:nummat
    load(loadstr{i}, 'XXprops');
    fprintf(fID, '    \\texttt{%s}	& %0.2e & %0.2e & %0.2e & %0.2e',...
        matstr{i}, XXprops.sv(1), XXprops.sv(end), XXprops.cond, XXprops.normF);
    if i < nummat
        fprintf(fID, '\\\\ \\hline \n');
    else % end table
        fprintf(fID, '\n');
    end
end
fprintf(fID, '\\end{tabular}\n');

fclose(fID);
open matrix_props
end
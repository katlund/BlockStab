function gen_plots(run_data, new_dir_str, var_eps)
% GEN_PLOTS(run_data, new_dir_str, var_eps) is a subroutine that generates
% plots for RUNKAPPAPLOT given a run_data struct; if new_dir_str is not
% provided, then the one saved in the .mat file specified by save_str is
% used.
%
% Part of [BlockStab](https://github.com/katlund/BlockStab) package.  Check README
% for how to properly cite and reuse this file.

%%
% Load run_data
n_alg = length(run_data.lgd);
if nargin == 1
    var_eps = eps;
elseif nargin == 2
    run_data.dir_str = new_dir_str;
    var_eps = eps;
end

% Pull out options
options = run_data.options;

% Initialize figures and axes
fg = cell(1,3);
ax = cell(1,3);
for k = 1:3
    fg{k} = figure;
    fg{k}.Position(3) = 700; % slightly stretches horizontally
    ax{k} = gca;
    hold on;
end

% Plot data
for j = 1:n_alg
    plot(ax{1}, run_data.condXX, run_data.loss_ortho(:,j),...
        run_data.symb{j}, 'Color', run_data.alg_cmap(j,:), ...
        'MarkerSize', 10, 'LineWidth', 1);
    plot(ax{2}, run_data.condXX, run_data.rel_res(:,j),... 
        run_data.symb{j}, 'Color', run_data.alg_cmap(j,:), ...
        'MarkerSize', 10, 'LineWidth', 1);
    plot(ax{3}, run_data.condXX, run_data.rel_chol_res(:,j),...
        run_data.symb{j}, 'Color', run_data.alg_cmap(j,:), ...
        'MarkerSize', 10, 'LineWidth', 1);
end
% Plot comparison lines
plot(ax{1}, run_data.condXX, var_eps*run_data.condXX, 'k--', ...
    run_data.condXX, var_eps*(run_data.condXX.^2), 'k-')
plot(ax{3}, run_data.condXX, var_eps*run_data.condXX, 'k--', ...
    run_data.condXX, var_eps*(run_data.condXX.^2), 'k-')

% Make plots pretty and save them
plot_str = {'loss_ortho', 'rel_res', 'rel_chol_res'};
for k = 1:3
    % Aesthetics
    set(ax{k}, 'Yscale', 'log', 'Xscale', 'log',...
        'XGrid', 'on', 'YGrid', 'on',...
        'XMinorGrid', 'off', 'YMinorGrid', 'off',...
        'FontSize', 12);
    
    % X-axis label
    xlabel(ax{k}, '$\kappa(\mathcal{X})$', ...
        'Interpreter', 'Latex', ...
        'FontSize', 12)
    
    % Legends and titles
    if k == 1
        lgd_loo = run_data.lgd;
        lgd_loo{end+1} = '$O(\varepsilon) \kappa(\mathcal{X})$'; %#ok<*AGROW> 
        lgd_loo{end+1} = '$O(\varepsilon) \kappa^2(\mathcal{X})$';
        legend(ax{k}, lgd_loo, 'Location', 'BestOutside', ...
            'Interpreter', 'Latex', ...
            'FontSize', 14, ...
            'EdgeColor','none', ...
            'Color','none');
        title(ax{k}, ...
            'Loss of Orthogonality', ...
            '$\Vert I - \bar\mathcal{Q}^T \bar\mathcal{Q}\Vert$', ...
            'Interpreter', 'Latex', ...
            'FontSize', 14);
        movegui(fg{1},'northwest')
    elseif k == 2
        legend(ax{k}, run_data.lgd, 'Location', 'BestOutside', ...
            'Interpreter', 'Latex', ...
            'FontSize', 14, ...
            'EdgeColor', 'none', ...
            'Color', 'none');
        title(ax{k}, ...
            'Relative Residual', ...
            ['$\Vert \mathcal{X} - ' ...
            '\bar\mathcal{Q}\bar\mathcal{R}\Vert/\Vert X\Vert$'], ...
            'Interpreter', 'Latex', ...
            'FontSize', 14);
        movegui(fg{2},'northeast')
    elseif k == 3
        lgd_chol = run_data.lgd;
        lgd_chol{end+1} = '$O(\varepsilon) \kappa(\mathcal{X})$'; %#ok<*AGROW> 
        lgd_chol{end+1} = '$O(\varepsilon) \kappa^2(\mathcal{X})$';
        legend(ax{k}, lgd_chol, 'Location', 'BestOutside', ...
            'Interpreter', 'Latex', ...
            'FontSize', 14, ...
            'EdgeColor', 'none', ...
            'Color', 'none');
        title(ax{k}, ...
            'Relative Cholesky Residual', ...
            ['$\Vert \mathcal{X}^T \mathcal{X} - ' ...
            '\bar\mathcal{R}^T \bar\mathcal{R}\Vert/ ' ...
            '\Vert \mathcal{X}\Vert^2$'], ...
            'Interpreter', 'Latex', ...
            'FontSize', 14);
        movegui(fg{3},'south')
    end

    % Save figures
    save_str = sprintf('%s/%s', run_data.dir_str, plot_str{k});
    if options.save_eps
        saveas(fg{k}, save_str, 'epsc');
    end

    if options.save_pdf
        exportgraphics(fg{k}, [save_str '.pdf'],'BackgroundColor','none')
    end

    if options.save_fig
        savefig(fg{k}, save_str, 'compact');
    end
end

% Print where figures are saved
if options.save_eps
    fprintf('EPS files saved in %s\n', run_data.dir_str);
end
if options.save_pdf
    fprintf('PDF files saved in %s\n', run_data.dir_str);
end
if options.save_fig
    fprintf('FIG files saved in %s\n', run_data.dir_str);
end

end
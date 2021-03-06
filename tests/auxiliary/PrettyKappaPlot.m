function PrettyKappaPlot(fg, ax, lgd_str, folder_str)
% PRETTYKAPPAPLOT(fg, ax, savestr) takes cells of figure and axis handles,
% as well as a save string, for a loss of orthogonality or residual plot
% generated by KAPPAPLOT and sets default plot properties, to make the plot
% objectively prettier.

%%
for i = 1:3
    axi = ax{i};
    fgi = fg{i};
    
    % Aesthetics
    set(axi, 'Yscale', 'log', 'Xscale', 'log',...
        'XGrid', 'on', 'YGrid', 'on',...
        'XMinorGrid', 'off', 'YMinorGrid', 'off');
    
    % X-axis label
    xlabel(axi, '\kappa(X)')
    
    % Legends and titles
    if i == 1
        lgd_str_loo = lgd_str;
        lgd_str_loo{end+1} = 'O(\epsilon) \kappa(X)';
        lgd_str_loo{end+1} = 'O(\epsilon) \kappa^2(X)';
        legend(axi, lgd_str_loo, 'Location', 'BestOutside');
        title(axi, 'Loss of Orthogonality: ||I - Q''*Q||');
        save_str = sprintf('%s/loss_ortho', folder_str);
    else
        legend(axi, lgd_str, 'Location', 'BestOutside');
        if i == 2
            title(axi, 'Relative Residual: ||X - Q*R||/||X||');
            save_str = sprintf('%s/res', folder_str);
        elseif i == 3
            title(axi, 'Relative Cholesky Residual: ||X''*X - R''*R||/||X||^2');
            save_str = sprintf('%s/res_chol', folder_str);
        end
    end
    savefig(fgi, save_str, 'compact');
    saveas(fgi, save_str, 'epsc')
end
end
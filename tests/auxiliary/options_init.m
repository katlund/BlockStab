function options = options_init(options)
% options = OPTIONS_INIT(options) initializes options for RUNKAPPAPLOT
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
if nargin == 1
    options = struct( ...
        'matrix_type', 'default', ...
        'scale', -(1:16), ...
        'num_rows', 100, ...
        'num_partitions', 10, ...
        'block_size', 2, ...
        'save_eps', false, ...
        'save_fig', false, ...
        'tex_report', false);
else
    if ~isfield(options, 'matrix_type')
        options.matrix_type = 'default';
    end

    if ~isfield(options, 'scale')
        switch matrix_type
            case 'default'
                options.scale = -(1:16);
                
            case 'glued'
                options.scale = 1:8;
                
            case 'laeuchli'
                options.scale = logspace(-1, -16, 16);
                
            case 'monomial'
                options.scale = 2:2:12;

        end
    end

    if ~isfield(options, 'num_rows')
        options.num_rows = 100;
    end

    if ~isfield(options, 'num_partitions')
        options.num_partitions = 10;
    end

    if ~isfield(options, 'block_size')
        options.block_size = 2;
    end

    if ~isfield(options, 'save_eps')
        options.save_eps = false;
    end

    if ~isfield(options, 'save_fig')
        options.save_fig = false;
    end

    if ~isfield(options, 'tex_report')
        options.tex_report = false;
    end
end

end
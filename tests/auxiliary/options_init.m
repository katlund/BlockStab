function options = options_init(mat_type, options)
% options = OPTIONS_INIT(mat_type, options) initializes options for
% RUNKAPPAPLOT.
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
if nargin == 0
    options = struct( ...
        'scale', -(1:16), ...
        'num_rows', 100, ...
        'num_partitions', 10, ...
        'block_size', 2, ...
        'save_eps', true, ...
        'save_fig', true, ...
        'tex_report', true);
elseif nargin >= 1
    if isempty(mat_type)
        mat_type = 'default';
    end

    if nargin == 1
        options = struct( ...
            'save_eps', true, ...
            'save_fig', true, ...
            'save_pdf', true, ...
            'tex_report', true);
    end

    if ~isfield(options, 'scale')
        switch mat_type
            case 'default'
                options.scale = -(1:16);
            case 'glued'
                options.scale = 1:8;
            case 'laeuchli'
                options.scale = logspace(-1, -16, 16);
            case 'monomial'
                options.scale = 2:2:12;
            case 'piled'
                options.scale = 1:12;
        end
    end

    if ~isfield(options, 'num_rows')
        switch mat_type
            case 'monomial'
                options.num_rows = 200;
            otherwise
                options.num_rows = 100;
        end
    end

    if ~isfield(options, 'num_partitions')
        switch mat_type
            case 'monomial'
                options.num_partitions = 20;
            otherwise
                options.num_partitions = 10;
        end
    end

    if ~isfield(options, 'block_size')
        switch mat_type
            case 'monomial'
                options.block_size = 6;
            otherwise
                options.block_size = 2;
        end
    end

    if ~isfield(options, 'save_eps')
        options.save_eps = true;
    end

    if ~isfield(options, 'save_pdf')
        options.save_pdf = true;
    end

    if ~isfield(options, 'save_fig')
        options.save_fig = true;
    end

    if ~isfield(options, 'tex_report')
        options.tex_report = true;
    end
end

end
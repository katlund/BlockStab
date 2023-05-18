function str = alg_string(skel, musc, param)
% str = ALG_STRING(skel, musc, param) converts the algorithm identifiers
% given by the cells skel, musc, and param into unique, legible strings for
% plots and reports.

%%
% Assume all inputs have same length are formatted as cell arrays
if ischar(skel)
    skel = {skel};
end
if ischar(musc)
    musc = {musc};
end
if ischar(param)
    param = {param};
end
n_alg = length(skel);
if n_alg ~= length(musc) || n_alg ~= length(param)
    error('Cell arrays skel, musc, and param must all have the same length.')
end

% Allocate space for str
str = cell(1, n_alg);

% Loop through all configurations
for i = 1:n_alg
    if isempty(skel{i})
        % Only format the muscle
        str{i} = sprintf('(%d) %s', i, basic_strrep(musc{i}));

    else
        % Skeleton-muscle
        if isempty(musc{i})
            switch skel{i}
                case {'bcgs_iro_ls', 'bcgs_iro_ls_mp', 'bcgs_sror'}
                    str{i} = sprintf('(%d) %s$', ...
                        i, basic_strrep(skel{i}));

                otherwise
                    str{i} = sprintf('(%d) %s$\\circ$HouseQR', ...
                        i, basic_strrep(skel{i}));
            end
        else
            str{i} = sprintf('(%d) %s$\\circ$%s', ...
                i, basic_strrep(skel{i}), basic_strrep(musc{i}));
        end
    end
end
end

%% Auxiliary
function str = basic_strrep(str)
    str = lower(str);
    if ~contains(str, 'qr')
        str = strrep(str, '_iro_', 'i+');
        str = strrep(str, '_iro', 'i+');
        str = strrep(str, '_sro', 's+');
        str = strrep(str, '_ro_', '+');
        str = strrep(str, 'ro_', '+');
        str = strrep(str, '_ro', '+');
        str = strrep(str, 'ro', '+');
        str = strrep(str, '_', '-');
        str = upper(str);

    else
        switch str
            case 'cholqr'
                str = 'CholQR';
            case 'iter_cholqr'
                str = 'IterCholQR';
            case 'sh_cholqr_roro'
                str = 'ShCholQR++';
            case 'cholqr_ro'
                str = 'CholQR+';
            case 'houseqr'
                str = 'HouseQR';
            case 'cholqr_pinv'
                str = 'CholQR-pinv';
            case 'globalqr'
                str = 'GlobalQR';
        end
    end
end
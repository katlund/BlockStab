function str = alg_string(skel, musc)
% str = ALG_STRING(skel, musc) converts the algorithm identifiers given by
% the cells skel and musc into unique, legible strings for plots and
% reports.
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
% Assume all inputs have same length are formatted as cell arrays
if ischar(skel)
    skel = {skel};
end
if ischar(musc)
    musc = {musc};
end
n_alg = length(skel);
if n_alg ~= length(musc)
    error('Cell arrays skel and musc must all have the same length.')
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
                    str{i} = sprintf('(%d) %s', ...
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
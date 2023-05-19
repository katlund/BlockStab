function alg_list = basic_strrep(alg_list)
% str = BASIC_STRREP(str) performs a basic string replacement for
% algorithm names in str.  See ALG_STRING.

%%
if ischar(alg_list)
    alg_list = char(alg_list);
end
for i = 1:length(alg_list)
    alg_list{i} = lower(alg_list{i});
    if ~contains(alg_list{i}, 'qr')
        alg_list{i} = strrep(alg_list{i}, '_iro_', 'i+');
        alg_list{i} = strrep(alg_list{i}, '_iro', 'i+');
        alg_list{i} = strrep(alg_list{i}, '_sro', 's+');
        alg_list{i} = strrep(alg_list{i}, '_ro_', '+');
        alg_list{i} = strrep(alg_list{i}, 'ro_', '+');
        alg_list{i} = strrep(alg_list{i}, '_ro', '+');
        alg_list{i} = strrep(alg_list{i}, 'ro', '+');
        alg_list{i} = strrep(alg_list{i}, '_', '-');
        alg_list{i} = upper(alg_list{i});

    else
        switch alg_list{i}
            case 'cholqr'
                alg_list{i} = 'CholQR';
            case 'iter_cholqr'
                alg_list{i} = 'IterCholQR';
            case 'sh_cholqr_roro'
                alg_list{i} = 'ShCholQR++';
            case 'cholqr_ro'
                alg_list{i} = 'CholQR+';
            case 'houseqr'
                alg_list{i} = 'HouseQR';
            case 'cholqr_pinv'
                alg_list{i} = 'CholQR-pinv';
            case 'globalqr'
                alg_list{i} = 'GlobalQR';
        end
    end
end
end
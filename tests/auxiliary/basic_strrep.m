function alg = basic_strrep(alg)
% alg = BASIC_STRREP(alg) performs a basic string replacement for
% algorithm names in str.  See ALG_STRING.
%
% Part of [BlockStab](https://github.com/katlund) package.  Check README
% for how to properly cite and reuse this file.

%%
if ischar(alg)
    alg = lower(alg);
    if ~contains(alg, 'qr')
        alg = strrep(alg, '_iro_', 'i+');
        alg = strrep(alg, '_iro', 'i+');
        alg = strrep(alg, '_sro', 's+');
        alg = strrep(alg, '_ro_', '+');
        alg = strrep(alg, 'ro_', '+');
        alg = strrep(alg, '_ro', '+');
        alg = strrep(alg, 'ro', '+');
        alg = upper(alg);
        alg = strrep(alg, '_MP', '$^{\rm{MP}}$');
        alg = strrep(alg, '+MP', '+$^{\rm{MP}}$');
        alg = strrep(alg, '_', '-');

    else
        switch alg
            case 'cholqr'
                alg = 'CholQR';
            case 'iter_cholqr'
                alg = 'IterCholQR';
            case 'sh_cholqr_roro'
                alg = 'ShCholQR++';
            case 'cholqr_ro'
                alg = 'CholQR+';
            case 'houseqr'
                alg = 'HouseQR';
            case 'cholqr_pinv'
                alg = 'CholQR-pinv';
            case 'globalqr'
                alg = 'GlobalQR';
        end
    end
else
    for i = 1:length(alg)
        alg{i} = lower(alg{i});
        if ~contains(alg{i}, 'qr')
            alg{i} = strrep(alg{i}, '_iro_', 'i+');
            alg{i} = strrep(alg{i}, '_iro', 'i+');
            alg{i} = strrep(alg{i}, '_sro', 's+');
            alg{i} = strrep(alg{i}, '_ro_', '+');
            alg{i} = strrep(alg{i}, 'ro_', '+');
            alg{i} = strrep(alg{i}, '_ro', '+');
            alg{i} = strrep(alg{i}, 'ro', '+');
            alg{i} = upper(alg{i});
            alg{i} = strrep(alg{i}, '_MP', '$^{\rm{MP}}$');
            alg{i} = strrep(alg{i}, '+MP', '+$^{\rm{MP}}$');
            alg{i} = strrep(alg{i}, '_', '-');
    
        else
            switch alg{i}
                case 'cholqr'
                    alg{i} = 'CholQR';
                case 'iter_cholqr'
                    alg{i} = 'IterCholQR';
                case 'sh_cholqr_roro'
                    alg{i} = 'ShCholQR++';
                case 'cholqr_ro'
                    alg{i} = 'CholQR+';
                case 'houseqr'
                    alg{i} = 'HouseQR';
                case 'cholqr_pinv'
                    alg{i} = 'CholQR-pinv';
                case 'globalqr'
                    alg{i} = 'GlobalQR';
            end
        end
    end
end
end
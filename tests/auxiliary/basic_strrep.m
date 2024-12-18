function alg = basic_strrep(alg)
% alg = BASIC_STRREP(alg) performs a basic string replacement for
% algorithm names in str.  See ALG_STRING.
%
% Part of [BlockStab](https://github.com/katlund/BlockStab) package.  Check README
% for how to properly cite and reuse this file.

%%
if ischar(alg) || isstruct(alg)
    alg = {alg};
end
for i = 1:length(alg)
    if isstruct(alg{i})
        % Avoid formatting user-provided strings-- too many options
        alg{i} = lower(alg{i}.id);
    else
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
if length(alg) == 1
    % Extract string for single entries
    alg = alg{1};
end
end
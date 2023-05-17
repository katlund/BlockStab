function alg_str = alg_string(alg)
% ALG_STRING(musc) converts the algorithm identifier given by alg into a
% more legible string for plots

%%
if ischar(alg)
    alg = {alg};
end
n_alg = length(alg);

alg_str = lower(alg);
for i = 1:n_alg
    str = alg_str{i};
    if ~contains(str, 'qr')
        str = strrep(str, '_free', '');
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
            case 'cholqr_mp'
                str = 'CholQR-MP';
            case 'cholqr_vpa'
                str = 'CholQR-VPA';
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
        end
    end
    alg_str{i} = str;
end
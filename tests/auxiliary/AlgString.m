function alg_str = AlgString(alg)
% ALGSTRING(musc) converts the algorithm identifier given by alg into a
% more legible string for plots

%%
if ischar(alg)
    alg = {alg};
end
nalg = length(alg);

alg_str = lower(alg);
for i = 1:nalg
    str = alg_str{i};
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
            case 'cholqr_free'
                str = 'CholQR-FREE';
        end
    end
    alg_str{i} = str;
end
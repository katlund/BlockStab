function A = mp_switch(A, mp_package, mp_spec)
% A = MP_SWITCH(A, mp_package, mp_spec) is a subroutine for switching between
% multiprecision packages for casting A to a desired precision.
% - mp_package: 'advanpix', 'symbolic math', or 'none'
% - mp_spec: 'single', 'double', or 'quad'; note that mp_package will be
%   ignored except whem mp_spec = 'quad'
%
% Part of [BlockStab](https://github.com/katlund/BlockStab) package.  Check README
% for how to properly cite and reuse this file.
  
%%
% Switch -- first level based on mp_spec
switch mp_spec
    case 'single'
        A = single(A);

    case 'double'
        A = double(A);

    case 'quad'
        % Defaults for both are quad precision
        switch mp_package
            case 'advanpix'
                A = mp(A);

            case {'symbolic math', 'symbolic toolbox', 'vpa'}
                A = vpa(A);

            otherwise
                error('Not a valid toolbox name')
        end

    otherwise
        error('Only single, double, and quad precisions are configured.')
end

end
function A = mp_switch(A, mp_package, mp_spec)
% A = MP_SWITCH(A, mp_package, mp_spec) is a subroutine for switching between
% multiprecision packages for casting A to a desired precision.
% - mp_package: 'advanpix', 'symbolic math', or 'none'
% - mp_spec: 'single', 'double', or 'quad'; note that mp_package will be
%   ignored except whem mp_spec = 'quad'
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).
    
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
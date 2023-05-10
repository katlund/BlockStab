function A = mp_switch(A, param)
% A = MP_SWITCH(A, param) is a subroutine for switching between
% mixed-precision packages for casting A to a desired precision.  Quadruple
% precision is the default, but other precisions can be specified in the
% struct param.
%
% param should specify the following fields:
% - .mp_package: 'advanpix', 'symbolic math', or 'none'
% - .mp_digits: the desired number of digits, according to the chosen
%   package specifications.  The default for 'advanpix' is 34 (quad) and
%   for 'symbolic math' is 32 (quad). 
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%%
% Default
if isempty(param)
    param.mp_package = 'none';
end

% Switch
switch param.mp_package

    case 'advanpix'
        if ~isfield(param, 'mp_digits')
            mp.Digits(34);
        else
            mp.Digits(param.mp_digits)
        end
        A =  mp(A);

    case {'symbolic math', 'symbolic toolbox', 'vpa'}
        if ~isfield(param, 'mp_digits')
            digits(34);
        else
            digits(param.mp_digits);
        end
        A = vpa(A);

    case 'none'
        % default double -- do nothing
end
end
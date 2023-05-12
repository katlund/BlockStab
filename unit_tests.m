% Unit tests
install_blockstab;

% Clear workspace
clear all; %#ok<*CLALL> 

%% See which toolboxes are installed
toolboxes = matlab.addons.installedAddons().Name;
if any(contains(toolboxes, 'Advanpix'))
    advanpix = true;
else
    advanpix = false;
end
if any(contains(toolboxes, 'Symbolic Math'))
    symmath = true;
else
    symmath = false;
end

%% MP_SWITCH
% Checks to make sure that qp is set up correctly.  MP versions should put
% qp(x) in higher precision and thus make 1 + qp(x) > 1.  Regular version
% should include higher rounding error, thus making 1 + qp(x) = 1.
x = 1e-16;

% Advanpix
if advanpix
    param.mp_package = 'advanpix';
    param.mp_digits = 34;
    qpa = @(x) mp_switch(x, param);
    
    fprintf('MP_SWITCH, ADVANPIX: %d\n', double(1+qpa(x) > 1));
end

% Symbolic Math Toolbox
if symmath
    param = [];
    param.mp_package = 'symbolic math';
    param.mp_digits = 32;
    qps = @(x) mp_switch(x, param);
    fprintf('MP_SWITCH, SYMBOLIC MATH: %d\n', double(1+qps(x) > 1));
end

% None
param = [];
param.mp_package = 'none';
qpn = @(x) mp_switch(x, param);
fprintf('MP_SWITCH, NONE: %d\n', 1+qpn(x) == 1);
fprintf('\n');

%% CHOL_SWITCH
% Check that CHOL_NAN returns NAN where expected, CHOL_FREE not, and
% CHOL_FREE with MP works.
X = 2*eye(2);

% CHOL_NAN
fprintf('CHOL_NAN(+) = SQRT: %d\n', norm(chol_nan(X) - sqrt(X)) == 0);
fprintf('CHOL_NAN(-) = NAN: %d\n', prod(isnan(chol_nan(-X)),'all'));

% CHOL_NAN w/ Advanpix
if advanpix
    Y = qpa(X);
    fprintf('CHOL_NAN(+) = SQRT, ADVANPIX: %d\n',...
        norm(chol_nan(Y) - sqrt(Y) == 0));
    fprintf('CHOL_NAN(-) = NAN, ADVANPIX: %d\n',...
        prod(isnan(chol_nan(-Y)),'all'));
end

% CHOL_NAN w/ Symbolic Math
if symmath
    Y = qps(X);
    fprintf('CHOL_NAN(+) = SQRT, SYMBOLIC MATH: %d\n',...
        norm(chol_nan(Y) - sqrt(Y)) == 0);
    fprintf('CHOL_NAN(-) = NAN, SYMBOLIC MATH: %d\n',...
        prod(isnan(chol_nan(-Y)),'all'));
end

% CHOL_FREE
fprintf('CHOL_FREE(+) = SQRT: %d\n', norm(chol_free(X) - sqrt(X)) == 0);
fprintf('CHOL_FREE(-) = SQRT: %d\n', norm(chol_free(-X) - sqrt(-X)) == 0);

% CHOL_FREE w/ Advanpix
if advanpix
    Y = qpa(X);
    param = [];
    param.mp_package = 'advanpix';
    param.mp_digits = 34;
    fprintf('CHOL_FREE(+) = SQRT, ADVANPIX: %d\n',...
        norm(chol_free(Y,param) - sqrt(Y)) == 0);
    fprintf('CHOL_FREE(-) = SQRT, ADVANPIX: %d\n',...
        norm(chol_free(-Y,param) - sqrt(-Y)) == 0);
end

% CHOL_FREE w/ Symbolic Math
if symmath
    Y = qps(X);
    param = [];
    param.mp_package = 'symbolic math';
    param.mp_digits = 32;
    fprintf('CHOL_FREE(+) = SQRT, SYMBOLIC MATH: %d\n',...
        norm(chol_free(Y,param) - sqrt(Y)) == 0);
    fprintf('CHOL_FREE(-) = SQRT, SYMBOLIC MATH: %d\n',...
        norm(chol_free(-Y,param) - sqrt(-Y)) == 0);
end

fprintf('\n');

%% Extract list of all muscles
musc_list = {dir('main\muscles\').name};
musc_list(1:2) = [];
musc_list(end+1:end+2) = {'global.m', 'global-no-scale.m'};

%% INTRAORTHO & INNER PROD
% Run through all options; ignore .m
n = 10; s = 5;
rng(4); X = rand(n,s);
rng(5); Y = rand(n,s);

% Standard double precision
param = [];
param.verbose = 1;
for i = 1:length(musc_list)
    musc = musc_list{i}(1:end-2);
    XY = InnerProd(X, Y, musc);
    switch musc
        case 'global'
            fprintf('INNERPROD, global: %d\n', ...
                norm(trace(X'*Y)/s - XY(1,1)) == 0);
        case 'global-no-scale'
            fprintf('INNERPROD, global-no-scale: %d\n', ...
                norm(trace(X'*Y) - XY(1,1)) == 0);
        otherwise
            fprintf('INNERPROD, %s: %d\n', ...
                musc, norm(X'*Y - XY) == 0);
    end
    if i <= length(musc_list) - 2
        if contains(musc,'qr')
            fprintf('INTRAORTHO, %s\n', musc);
        else
            fprintf('INTRAORTHO, %s: Evaluate LOO and RelRes for irregularites.\n', musc);
        end
    end
    IntraOrtho(X, musc, param);
    fprintf('\n');
end

% Quadruple precision
param = [];
param.verbose = 1;

if advanpix
    param = [];
    param.verbose = 1;
    param.mp_package = 'advanpix';
    param.mp_digits = 34;
    for i = 1:length(musc_list)
        musc = musc_list{i}(1:end-2);
        XY = InnerProd(X, Y, musc);
        switch musc
            case 'global'
                fprintf('INNERPROD, global, advanpix: %d\n', ...
                    norm(trace(X'*Y)/s - XY(1,1)) == 0);
            case 'global-no-scale'
                fprintf('INNERPROD, global-no-scale, advanpix: %d\n', ...
                    norm(trace(X'*Y) - XY(1,1)) == 0);
            otherwise
                fprintf('INNERPROD, %s, advanpix: %d\n', ...
                    musc, norm(X'*Y - XY) == 0);
        end
        if i <= length(musc_list) - 2
            if contains(musc,'qr')
                fprintf('INTRAORTHO, %s, advanpix\n', musc);
            else
                fprintf('INTRAORTHO, %s, advanpix: Evaluate LOO and RelRes for irregularites.\n', musc);
            end
        end
        IntraOrtho(X, musc, param);
        fprintf('\n');
    end
end

if symmath
    param = [];
    param.verbose = 1;
    param.mp_package = 'symbolic math';
    param.mp_digits = 32;
    for i = 1:length(musc_list)
        musc = musc_list{i}(1:end-2);
        XY = InnerProd(X, Y, musc);
        switch musc
            case 'global'
                fprintf('INNERPROD, global, symbolic math: %d\n', ...
                    norm(trace(X'*Y)/s - XY(1,1)) == 0);
            case 'global-no-scale'
                fprintf('INNERPROD, global-no-scale, symbolic math: %d\n', ...
                    norm(trace(X'*Y) - XY(1,1)) == 0);
            otherwise
                fprintf('INNERPROD, %s, symbolic math: %d\n', ...
                    musc, norm(X'*Y - XY) == 0);
        end
        if i <= length(musc_list) - 2
            if contains(musc,'qr')
                fprintf('INTRAORTHO, %s, symbolic math\n', musc);
            else
                fprintf('INTRAORTHO, %s, symbolic math: Evaluate LOO and RelRes for irregularites.\n', musc);
            end
        end
        IntraOrtho(X, musc, param);
        fprintf('\n');
    end
end

%% BGS
skel_list = {dir('main\skeletons\').name};
skel_list(1:2) = [];

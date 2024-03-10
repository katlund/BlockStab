%% Unit tests
%
% Part of the BlockStab package documented in [Carson, et al.
% 2022](https://doi.org/10.1016/j.laa.2021.12.017).

%% Set-up
install_blockstab;

% Turn off warnings
warning('off')

% Clear workspace
clear all; %#ok<*CLALL>

% Open .txt file for writing report
report_name = 'unit_tests_report.txt';
fID = fopen(report_name,'w');

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
% Checks to make sure that precisions are set up correctly.
x = 1;

% Check single
sp = @(x) mp_switch(x, [], 'single');
fprintf(fID, 'MP_SWITCH, SINGLE: %d\n', isa(sp(x), 'single'));

% Check double
dp = @(x) mp_switch(x, [], 'double');
fprintf(fID, 'MP_SWITCH, DOUBLE: %d\n', isa(dp(x), 'double'));

% Check quad
% Advanpix
if advanpix
    param.mp_package = 'advanpix';
    qpa = @(x) mp_switch(x, param.mp_package, 'quad');
    fprintf(fID, 'MP_SWITCH, ADVANPIX, QUAD: %d\n', isa(qpa(x), 'mp'));
end

% Symbolic Math Toolbox
if symmath
    param = [];
    param.mp_package = 'symbolic math';
    qps = @(x) mp_switch(x, param.mp_package, 'quad');
    fprintf(fID, 'MP_SWITCH, SYMBOLIC MATH, QUAD: %d\n', isa(qps(x), 'sym'));
end

%% CHOL_SWITCH
% Check that CHOL_NAN returns NAN where expected, CHOL_FREE not, and
% CHOL_FREE with MP works.
X = 2*eye(2);

% CHOL_NAN
fprintf(fID, 'CHOL_NAN(+) = SQRT: %d\n', norm(chol_nan(X) - sqrt(X)) == 0);
fprintf(fID, 'CHOL_NAN(-) = NAN: %d\n', prod(isnan(chol_nan(-X)),'all'));

% CHOL_NAN w/ Advanpix
if advanpix
    Y = qpa(X);
    fprintf(fID, 'CHOL_NAN(+) = SQRT, ADVANPIX: %d\n',...
        norm(chol_nan(Y) - sqrt(Y) == 0));
    fprintf(fID, 'CHOL_NAN(-) = NAN, ADVANPIX: %d\n',...
        prod(isnan(chol_nan(-Y)),'all'));
end

% CHOL_NAN w/ Symbolic Math
if symmath
    Y = qps(X);
    fprintf(fID, 'CHOL_NAN(+) = SQRT, SYMBOLIC MATH: %d\n',...
        norm(chol_nan(Y) - sqrt(Y)) == 0);
    fprintf(fID, 'CHOL_NAN(-) = NAN, SYMBOLIC MATH: %d\n',...
        prod(isnan(chol_nan(-Y)),'all'));
end

% CHOL_FREE
fprintf(fID, 'CHOL_FREE(+) = SQRT: %d\n', norm(chol_free(X) - sqrt(X)) == 0);
fprintf(fID, 'CHOL_FREE(-) = SQRT: %d\n', norm(chol_free(-X) - sqrt(-X)) == 0);

% CHOL_FREE w/ Advanpix
if advanpix
    Y = qpa(X);
    param = [];
    param.mp_package = 'advanpix';
    param.mp_pair = {'double', 'quad'};
    fprintf(fID, 'CHOL_FREE(+) = SQRT, ADVANPIX: %d\n',...
        norm(chol_free(Y,param) - sqrt(Y)) == 0);
    fprintf(fID, 'CHOL_FREE(-) = SQRT, ADVANPIX: %d\n',...
        norm(chol_free(-Y,param) - sqrt(-Y)) == 0);
end

% CHOL_FREE w/ Symbolic Math
if symmath
    Y = qps(X);
    param = [];
    param.mp_package = 'symbolic math';
    param.mp_pair = {'double', 'quad'};
    fprintf(fID, 'CHOL_FREE(+) = SQRT, SYMBOLIC MATH: %d\n',...
        norm(chol_free(Y,param) - sqrt(Y)) == 0);
    fprintf(fID, 'CHOL_FREE(-) = SQRT, SYMBOLIC MATH: %d\n',...
        norm(chol_free(-Y,param) - sqrt(-Y)) == 0);
end

fprintf(fID, '\n');

%% Extract list of all muscles
musc_list = {dir('main/muscles/').name};
musc_list(1:2) = [];
musc_list(end+1:end+2) = {'global.m', 'global-no-scale.m'};

%% INTRAORTHO & INNER PROD
% Run through all options; ignore .m
n = 10; s = 3;
rng(4); X = rand(n,s);
rng(5); Y = rand(n,s);

% Standard double precision
param = [];
param.verbose = 1;
for i = 1:length(musc_list)
    musc = musc_list{i}(1:end-2);

    % InnerProd
    XY = InnerProd(X, Y, musc);
    switch musc
        case 'global'
            fprintf(fID, 'INNERPROD, GLOBAL: %d\n', ...
                norm(trace(X'*Y)/s - XY(1,1)) == 0);
        case 'global-no-scale'
            fprintf(fID, 'INNERPROD, GLOBAL-NO-SCALE: %d\n', ...
                norm(trace(X'*Y) - XY(1,1)) == 0);
        otherwise
            fprintf(fID, 'INNERPROD, %s: %d\n', ...
                upper(musc), norm(X'*Y - XY) == 0);
    end

    % IntraOrtho
    if i <= length(musc_list) - 2
        fprintf(fID, 'INTRAORTHO, %s\n', upper(musc));
        try
            [Q, R, T] = IntraOrtho(X, musc, param);
            fprintf(fID, 'PASS\n');
        catch ME
            fprintf(fID, 'FAIL: %s\n', ME.message);
        end
        fprintf(fID, '\n');
    end
end

% Multiprecision
musc_list(end-1:end) = []; % remove global and global-no-scale
param = [];
param.verbose = 1;

if advanpix
    param = [];
    param.verbose = 1;
    param.mp_package = 'advanpix';
    param.mp_pair = {'double', 'quad'};
    for i = 1:length(musc_list)
        musc = musc_list{i}(1:end-2);
        fprintf(fID, 'INTRAORTHO, %s, ADVANPIX\n', upper(musc));
        try
            [Q, R, T] = IntraOrtho(X, musc, param);
            fprintf(fID, 'PASS\n');
        catch ME
            fprintf(fID, 'FAIL: %s\n', ME.message);
        end
        fprintf(fID, '\n');
    end
end

if symmath
    param = [];
    param.verbose = 1;
    param.mp_package = 'symbolic math';
    param.mp_pair = {'double', 'quad'};
    for i = 1:length(musc_list)
        musc = musc_list{i}(1:end-2);
        fprintf(fID, 'INTRAORTHO, %s, SYMBOLIC MATH\n', upper(musc));
        try
            [Q, R, T] = IntraOrtho(X, musc, param);
            fprintf(fID, 'PASS\n');
        catch ME
            fprintf(fID, 'FAIL: %s\n', ME.message);
        end
        fprintf(fID, '\n');
    end
end

%% BGS
skel_list = {dir('main/skeletons/').name};
skel_list(1:2) = [];

n = 10; s = 2; p = 3;
rng(4); XX = rand(n,s*p);

% Look at Cartesian product between all skeletons and muscles (standard
% double precision)
param = [];
param.verbose = 1;
for j = 1:length(skel_list)
    skel = skel_list{j}(1:end-2);
    if strcmp(skel, 'bcgs_sror')
        % Skip unnecessary re-runs
        musc = 'cgs_sror';
        fprintf(fID, 'BGS, %s-%s:.\n', upper(skel), upper(musc));
        try
            [QQ, RR, TT] = BGS(XX, s, skel, musc, param);
            fprintf(fID, 'PASS\n');
        catch ME
            fprintf(fID, 'FAIL: %s\n', ME.message);
        end
        fprintf(fID, '\n');
    else
        if ~contains(skel, '_mp')
            for i = 1:length(musc_list)
                musc = musc_list{i}(1:end-2);
                fprintf(fID, 'BGS, %s-%s:\n', upper(skel), upper(musc));
                try
                    [QQ, RR, TT] = BGS(XX, s, skel, musc, param);
                    fprintf(fID, 'PASS\n');
                catch ME
                    fprintf(fID, 'FAIL: %s\n', ME.message);
                end
                fprintf(fID, '\n');
            end
        end
    end
end

% Look at Cartesian product between all skeletons and muscles
% multiprecision)
param = [];
param.verbose = 1;
if advanpix
    param.mp_package = 'advanpix';
    param.mp_pair = {'double', 'quad'};
    for j = 1:length(skel_list)
        skel = skel_list{j}(1:end-2);
        if strcmp(skel, 'bcgs_sror')
            % Skip unnecessary re-runs
            musc = 'cgs_sror';
            fprintf(fID, 'BGS, %s-%s, ADVANPIX:\n', upper(skel), upper(musc));
            try
                [QQ, RR, TT] = BGS(XX, s, skel, musc, param);
                fprintf(fID, 'PASS\n');
            catch ME
                fprintf(fID, 'FAIL: %s\n', ME.message);
            end
        else
            for i = 1:length(musc_list)
                musc = musc_list{i}(1:end-2);
                fprintf(fID, 'BGS, %s-%s, ADVANPIX:\n', upper(skel), upper(musc));
                try
                    [QQ, RR, TT] = BGS(XX, s, skel, musc, param);
                    fprintf(fID, 'PASS\n');
                catch ME
                    fprintf(fID, 'FAIL: %s\n', ME.message);
                end
                fprintf(fID, '\n');
            end
        end
    end
end

if symmath
    param.mp_package = 'symbolic math';
    param.mp_pair = {'double', 'quad'};
    for j = 1:length(skel_list)
        skel = skel_list{j}(1:end-2);
        if strcmp(skel, 'bcgs_sror')
            % Skip unnecessary re-runs
            musc = 'cgs_sror';
            fprintf(fID, 'BGS, %s-%s, SYMBOLIC MATH:\n', upper(skel), upper(musc));
            try
                [QQ, RR, TT] = BGS(XX, s, skel, musc, param);
                fprintf(fID, 'PASS\n');
            catch ME
                fprintf(fID, 'FAIL: %s\n', ME.message);
            end
            fprintf(fID, '\n');
        else
            for i = 1:length(musc_list)
                musc = musc_list{i}(1:end-2);
                fprintf(fID, 'BGS, %s-%s, SYMBOLIC MATH:\n', upper(skel), upper(musc));
                try
                    [QQ, RR, TT] = BGS(XX, s, skel, musc, param);
                    fprintf(fID, 'PASS\n');
                catch ME
                    fprintf(fID, 'FAIL: %s\n', ME.message);
                end
                fprintf(fID, '\n');
            end
        end
    end
end

%% Close file
fclose(fID);

%% Generate RunKappaPlot reports for all matrix types
cd tests
mat_type = {'default', 'glued', 'laeuchli', 'monomial', 'piled'};
for i = 1:4
    RunKappaPlot(mat_type{i});
    close all;
end

%% Basic MakeHeatmap test
MakeHeatmap([100 10 2], 'stewart', {'BCGS', 'BCGS_IRO', 'BCGS_SROR'}, {'CGS', 'HouseQR', 'CGS_SROR'}, 1, 1)

cd ..
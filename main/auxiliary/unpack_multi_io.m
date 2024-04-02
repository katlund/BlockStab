function [musc, musc_param] = unpack_multi_io(multi_io, param)
% [musc, musc_param] = UNPACK_MULTI_IO(multi_io, param) converts a multiIO
% struct into two cell arrays ordered according to IO_A, IO_1, and IO_2,
% respectively.  The second cell array corresponds to the musc parameters;
% they are augmented by param (which should correspond to the skel level)
% if provided.
%
% multi_io may have top-level fields io_a, io_1, and io_2.  io_a is always
% required, and either io_1 or io_2 (or both) must additionally be
% provided.  When one is missing, io_1 = io_2 is set as default.

%%
% Errors
if ~isfield(multi_io, 'io_a')
    error('io_a is missing')
end
if ~(isfield(multi_io, 'io_1') || isfield(multi_io, 'io_2'))
    error('either io_1 or io_2 must be provided')
end

% Defaults
if nargin == 1
    param = [];
end

% Initialize musc and param
musc = cell(1,3);
musc_param = cell(1,3);

% IO_A
f = fieldnames(multi_io.io_a);
musc{1} = f{1};

if ~isempty(multi_io.io_a.(musc{1})) && nargin == 2
    musc_param{1} = catstruct(multi_io.io_a.(musc{1}).param, param);
else
    musc_param{1} = param;
end

% IO_1
if isfield(multi_io, 'io_1')
    f = fieldnames(multi_io.io_1);
    musc{2} = f{1};
    
    if ~isempty(multi_io.io_1.(musc{2})) && nargin == 2
        musc_param{2} = catstruct(multi_io.io_1.(musc{2}).param, param);
    else
        musc_param{2} = param;
    end
    flag_1 = true;
else
    flag_1 = false;
end

% IO_2
if isfield(multi_io, 'io_2')
    f = fieldnames(multi_io.io_2);
    musc{3} = f{1};
    
    if ~isempty(multi_io.io_2.(musc{3})) && nargin == 2
        musc_param{3} = catstruct(multi_io.io_2.(musc{3}).param, param);
    else
        musc_param{3} = param;
    end
else
    % Copy IO_1
    musc{3} = musc{2};
    musc_param{3} = musc_param{2};
end

% Copy IO_2
if ~flag_1
    musc{2} = musc{3};
    musc_param{2} = musc_param{3};
end
end
function S = setFieldPath(S, path, value)
% SETFIELDPATH Set a (possibly nested) struct field using a dot-separated path.
%
%   S = SETFIELDPATH(S, "a.b.c", value) sets S.a.b.c = value, creating
%   intermediate structs as needed.
%
%   This helper is used for configuration overrides (sweeps/experiments).
%
%   2026 NIST/CTL Steve Blandino
%   This file is available under the terms of the NIST License.

if isstring(path)
    path = char(path);
end

if ~ischar(path) || isempty(path)
    error('setFieldPath:InvalidPath', 'PATH must be a non-empty char/string.');
end

parts = strsplit(path, '.');
S = iSetRec(S, parts, 1, value);

end

function S = iSetRec(S, parts, idx, value)
f = parts{idx};

if idx == numel(parts)
    S.(f) = value;
    return;
end

if ~isfield(S, f) || ~isstruct(S.(f))
    S.(f) = struct();
end

S.(f) = iSetRec(S.(f), parts, idx+1, value);
end

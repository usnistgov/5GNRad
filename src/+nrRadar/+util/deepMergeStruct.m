function out = deepMergeStruct(base, override)
% DEEPMERGESTRUCT Deep-merge two structs.
%
%   OUT = DEEPMERGESTRUCT(BASE, OVERRIDE) returns BASE with fields from
%   OVERRIDE applied. If a field exists in both and both values are structs,
%   merge recursively. Otherwise, OVERRIDE wins.
%
%   This is used to apply configuration overrides without changing original
%   field names.
%
%   2026 NIST/CTL Steve Blandino
%   This file is available under the terms of the NIST License.

if isempty(base)
    out = override;
    return;
end

if isempty(override)
    out = base;
    return;
end

if ~isstruct(base) || ~isstruct(override)
    error('deepMergeStruct:InvalidInput', 'BASE and OVERRIDE must be structs (or empty).');
end

out = base;
f = fieldnames(override);
for k = 1:numel(f)
    name = f{k};
    if isfield(base, name) && isstruct(base.(name)) && isstruct(override.(name))
        out.(name) = nrRadar.util.deepMergeStruct(base.(name), override.(name));
    else
        out.(name) = override.(name);
    end
end
end

function [simulation, target, prs, geometry, sens] = applyOverrides(...
    simulation, target, prs, geometry, sens, overrides)
% APPLYOVERRIDES Apply overrides to scenario configuration structs.
%
%   [SIMULATION, TARGET, PRS, GEOMETRY, SENS] = APPLYOVERRIDES(...)
%   applies user-provided overrides without requiring new scenario folders.
%
%   OVERRIDES forms supported:
%   1) N-by-2 cell array of {PATH, VALUE} rows, where PATH is a dot-separated
%      string like:
%         'sens.cfarThreshold'
%         'simulation.nMaxDrop'
%      Prefix aliases are accepted:
%         simConfig -> simulation
%         stConfig  -> target
%         prsConfig -> prs
%         sensConfig-> sens
%   2) name/value cell array:
%         { 'sens.cfarThreshold', 10, 'simulation.nMaxDrop', 20 }
%   3) struct with top-level fields simulation/target/prs/geometry/sens (or
%      simConfig/stConfig/prsConfig/sensConfig). These are deep-merged.
%
%   Notes
%   - This function does NOT rename any existing fields; it only sets/merges.
%
%   2026 NIST/CTL Steve Blandino
%   This file is available under the terms of the NIST License.

if isempty(overrides)
    return;
end

% ---- Struct form (deep merge) ----
if isstruct(overrides)
    if isfield(overrides, 'simConfig')
        overrides.simulation = overrides.simConfig;
        overrides = rmfield(overrides, 'simConfig');
    end
    if isfield(overrides, 'stConfig')
        overrides.target = overrides.stConfig;
        overrides = rmfield(overrides, 'stConfig');
    end
    if isfield(overrides, 'prsConfig')
        overrides.prs = overrides.prsConfig;
        overrides = rmfield(overrides, 'prsConfig');
    end
    if isfield(overrides, 'sensConfig')
        overrides.sens = overrides.sensConfig;
        overrides = rmfield(overrides, 'sensConfig');
    end

    if isfield(overrides, 'simulation')
        simulation = nrRadar.util.deepMergeStruct(simulation, overrides.simulation);
    end
    if isfield(overrides, 'target')
        target = nrRadar.util.deepMergeStruct(target, overrides.target);
    end
    if isfield(overrides, 'prs')
        prs = nrRadar.util.deepMergeStruct(prs, overrides.prs);
    end
    if isfield(overrides, 'geometry')
        geometry = nrRadar.util.deepMergeStruct(geometry, overrides.geometry);
    end
    if isfield(overrides, 'sens')
        sens = nrRadar.util.deepMergeStruct(sens, overrides.sens);
    end
    return;
end

% ---- Cell form (path/value) ----
if iscell(overrides)
    if isvector(overrides) && mod(numel(overrides),2) == 0
        % name/value vector -> make N-by-2
        overrides = reshape(overrides, 2, []).';
    end

    if size(overrides,2) ~= 2
        error('applyOverrides:InvalidCell', ...
            'Cell overrides must be N-by-2 {PATH,VALUE} or name/value vector.');
    end

    for k = 1:size(overrides,1)
        path = overrides{k,1};
        value = overrides{k,2};
        [simulation, target, prs, geometry, sens] = iApplyOne(...
            simulation, target, prs, geometry, sens, path, value);
    end
    return;
end

error('applyOverrides:InvalidType', ...
    'Overrides must be a struct, a cell array, or empty.');

end

function [simulation, target, prs, geometry, sens] = iApplyOne(...
    simulation, target, prs, geometry, sens, path, value)

if isstring(path)
    path = char(path);
end

if ~ischar(path) || isempty(path)
    error('applyOverrides:InvalidPath', 'Override path must be a non-empty char/string.');
end

% Split only on the first dot
dotIdx = find(path=='.', 1, 'first');
if isempty(dotIdx)
    error('applyOverrides:MissingPrefix', ...
        'Override path must include a prefix, e.g., "sens.cfarThreshold".');
end

prefix = path(1:dotIdx-1);
rest   = path(dotIdx+1:end);

% Accept common aliases (keep original variable naming in the codebase)
switch prefix
    case 'simConfig'
        prefix = 'simulation';
    case 'stConfig'
        prefix = 'target';
    case 'prsConfig'
        prefix = 'prs';
    case 'sensConfig'
        prefix = 'sens';
end

switch prefix
    case 'simulation'
        simulation = nrRadar.util.setFieldPath(simulation, rest, value);
    case 'target'
        target = nrRadar.util.setFieldPath(target, rest, value);
    case 'prs'
        prs = nrRadar.util.setFieldPath(prs, rest, value);
    case 'geometry'
        geometry = nrRadar.util.setFieldPath(geometry, rest, value);
    case 'sens'
        sens = nrRadar.util.setFieldPath(sens, rest, value);
    otherwise
        error('applyOverrides:UnknownPrefix', ...
            'Unknown override prefix "%s". Use simulation/target/prs/geometry/sens.', prefix);
end

end

function runs = makeRunsGrid(paths, valueCells, varargin)
% MAKERUNSGRID Convenience helper to build a grid sweep across parameters.
%
%   RUNS = MAKERUNSGRID(PATHS, VALUECELLS) creates the Cartesian product of
%   the values provided in VALUECELLS.
%
%   Inputs
%     paths      - cellstr/string array, each entry is a parameter path
%                 (e.g., {'sens.cfarThreshold','sens.nmsRadius'})
%     valueCells - cell array where valueCells{i} is a vector/cell array of values for paths{i}
%
%   NAME-VALUE
%     'TagPrefix' - prefix for tags (default: 'run')
%
%   Example
%     runs = nrRadar.exp.makeRunsGrid( ...
%         {'sens.cfarThreshold','sens.cfarTrnCellRange'}, ...
%         {10.^((0:10)/10), [4 8 12]}, ...
%         'TagPrefix','grid');
%
%   2026 NIST/CTL Steve Blandino
%   This file is available under the terms of the NIST License.

p = inputParser;
addParameter(p, 'TagPrefix', 'run');
parse(p, varargin{:});
opt = p.Results;

paths = cellstr(paths);
if ~iscell(valueCells)
    error('makeRunsGrid:InvalidInput', 'valueCells must be a cell array.');
end
if numel(paths) ~= numel(valueCells)
    error('makeRunsGrid:SizeMismatch', 'paths and valueCells must have same length.');
end

% Convert each value list to a cell array (so values can be non-scalar, strings, etc.)
vals = cell(size(valueCells));
for i = 1:numel(valueCells)
    vi = valueCells{i};
    if iscell(vi)
        vals{i} = vi(:);
    else
        vals{i} = num2cell(vi(:));
    end
end

% Build Cartesian product indices
nPer = cellfun(@numel, vals);
idxGrids = cell(1, numel(vals));

% ndgrid expects per-dimension index vectors
idxCells = cell(1, numel(vals));
for i = 1:numel(vals)
    idxCells{i} = 1:nPer(i);
end
[idxGrids{:}] = ndgrid(idxCells{:});

nRuns = numel(idxGrids{1});
runs = repmat(struct('tag', "", 'set', []), nRuns, 1);

for k = 1:nRuns
    setCell = cell(numel(paths), 2);
    tagParts = strings(1, numel(paths));
    for i = 1:numel(paths)
        idx = idxGrids{i}(k);
        v = vals{i}{idx};
        setCell{i,1} = paths{i};
        setCell{i,2} = v;

        tagParts(i) = iValToTag(paths{i}, v);
    end
    runs(k).set = setCell;
    runs(k).tag = string(opt.TagPrefix) + "_" + join(tagParts, "__");
end

end

function s = iValToTag(path, v)
% Build a short tag component that remains filesystem-friendly.
parts = strsplit(string(path), '.');
name = parts{end};

if ischar(v) || isstring(v)
    vv = string(v);
else
    try
        vv = string(mat2str(v));
    catch
        vv = "val";
    end
end

vv = regexprep(vv, '\s+', '');
vv = regexprep(vv, '[\[\]\(\)\,\;]', '_');
vv = regexprep(vv, '[^A-Za-z0-9_\-\.]', '');
s = name + "=" + vv;
end

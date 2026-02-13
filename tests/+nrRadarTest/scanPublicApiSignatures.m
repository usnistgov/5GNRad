function api = scanPublicApiSignatures(srcRoot)
%SCANPUBLICAPISIGNATURES Return signatures for public functions under +nrRadar.
% Public = everything under src/+nrRadar except +internal.
%
% Output: struct array with fields:
%   - name: fully qualified function name (e.g., "nrRadar.rx.getRxWaveform")
%   - file: relative file path
%   - narginMin/narginMax: inferred from signature parsing (best-effort)
%   - nargoutMin/nargoutMax: best-effort (based on bracketed outputs)
%   - raw: raw function line

repoRoot = fileparts(fileparts(fileparts(srcRoot)));
mfiles = dir(fullfile(srcRoot, "**", "*.m"));

keep = true(size(mfiles));
for i = 1:numel(mfiles)
    rel = string(fullfile(mfiles(i).folder, mfiles(i).name));
    rel = erase(rel, repoRoot + filesep);
    % Exclude internal package
    if contains(rel, "+nrRadar/+internal/")
        keep(i) = false;
    end
end
mfiles = mfiles(keep);

api = struct("name", {}, "file", {}, "narginMin", {}, "narginMax", {}, "nargoutMin", {}, "nargoutMax", {}, "raw", {});

for i = 1:numel(mfiles)
    fpath = fullfile(mfiles(i).folder, mfiles(i).name);

    % Skip classdef files
    if localIsClassdef(fpath)
        continue;
    end

    fline = localFirstFunctionLine(fpath);
    if strlength(fline) == 0
        continue;
    end

    [name, nmin, nmax, om, ox] = localParseSignature(fline, fpath, repoRoot);
    if strlength(name) == 0
        continue;
    end

    relFile = string(erase(fpath, repoRoot + filesep));
    api(end+1) = struct( ...
        "name", name, ...
        "file", relFile, ...
        "narginMin", nmin, ...
        "narginMax", nmax, ...
        "nargoutMin", om, ...
        "nargoutMax", ox, ...
        "raw", fline); %#ok<AGROW>
end

% Stable ordering
[~,idx] = sort(string({api.name}));
api = api(idx);
end

function tf = localIsClassdef(fpath)
txt = fileread(fpath);
tf = ~isempty(regexp(txt, "^\s*classdef\b", "once", "lineanchors"));
end

function fline = localFirstFunctionLine(fpath)
lines = splitlines(string(fileread(fpath)));
fline = "";
for k = 1:numel(lines)
    L = strtrim(lines(k));
    if startsWith(L, "function ")
        fline = L;
        return;
    end
    % stop early if we hit non-comment, non-empty before function
    if strlength(L) > 0 && ~startsWith(L, "%")
        return;
    end
end
end

function [fqName, narginMin, narginMax, nargoutMin, nargoutMax] = localParseSignature(fline, fpath, repoRoot)
% Best-effort parsing of: function [a,b]=pkg.fun(x,y,varargin)

fqName = "";
narginMin = NaN; narginMax = NaN;
nargoutMin = NaN; nargoutMax = NaN;

% Extract outputs (optional)
outTok = regexp(fline, "^function\s*(?<outs>\[[^\]]*\]|[^\s=]+)?\s*=\s*(?<rest>.+)$", "names");
if isempty(outTok)
    % Maybe no '=' form: function foo(x)
    outTok2 = regexp(fline, "^function\s*(?<rest>.+)$", "names");
    if isempty(outTok2), return; end
    outs = "";
    rest = strtrim(outTok2.rest);
else
    outs = strtrim(string(outTok.outs));
    rest = strtrim(string(outTok.rest));
end

% rest begins with function name and args
nameTok = regexp(rest, "^(?<fname>[A-Za-z]\w*)\s*(\((?<args>.*)\))?\s*$", "names");
if isempty(nameTok), return; end

funName = string(nameTok.fname);
argsStr = "";
if isfield(nameTok,"args") && ~isempty(nameTok.args)
    argsStr = string(nameTok.args);
end

% Determine fully qualified name from folder structure
rel = string(erase(fpath, repoRoot + filesep));
% rel like src/+nrRadar/+rx/getRxWaveform.m
parts = split(rel, filesep);
% Find +nrRadar and subsequent +pkg parts
idx = find(parts == "src", 1, "first");
if isempty(idx), return; end
% Locate +nrRadar
idxNR = find(contains(parts, "+nrRadar"), 1, "first");
if isempty(idxNR), return; end
pkgParts = parts(idxNR:end-1); % exclude file
pkg = "nrRadar";
for p = 2:numel(pkgParts)
    s = pkgParts(p);
    if startsWith(s, "+")
        pkg = pkg + "." + extractAfter(s, 1);
    end
end
fqName = pkg + "." + funName;

% nargin parsing (best-effort)
if strlength(argsStr) == 0
    narginMin = 0; narginMax = 0;
else
    args = strtrim(split(argsStr, ","));
    args(args=="") = [];
    hasVarargin = any(args == "varargin");
    hasVarargout = any(args == "varargout"); %#ok<NASGU>

    narginMin = sum(args ~= "varargin");
    narginMax = narginMin;
    if hasVarargin
        narginMax = Inf;
    end
end

% nargout parsing (best-effort)
if strlength(outs) == 0
    nargoutMin = 0; nargoutMax = 0;
elseif startsWith(outs,"[")
    inner = strtrim(extractBetween(outs, "[", "]"));
    if isempty(inner)
        nargoutMin = 0; nargoutMax = 0;
    else
        o = strtrim(split(inner, ","));
        o(o=="")=[];
        nargoutMin = numel(o);
        nargoutMax = nargoutMin;
    end
else
    nargoutMin = 1; nargoutMax = 1;
end
end

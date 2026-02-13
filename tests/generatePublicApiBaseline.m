function generatePublicApiBaseline(varargin)
% GENERATEPUBLICAPIBASELINE Generate the public API signature baseline.
%
% This writes a JSON file under tests/compat/publicApiSignatures.json that
% records function signatures for the public MATLAB API in src/+nrRadar,
% excluding src/+nrRadar/+internal.
%
% Usage:
%   cd(repoRoot); setup; addpath("tests");
%   generatePublicApiBaseline
%
% Optional name-value:
%   'OutFile' : output JSON path
%
% 2026 NIST/CTL
% This file is available under the terms of the NIST License.

p = inputParser;
addParameter(p, 'OutFile', "");
parse(p, varargin{:});

repoRoot = fileparts(fileparts(mfilename('fullpath')));
if strlength(p.Results.OutFile) == 0
    outFile = fullfile(repoRoot, "tests", "compat", "publicApiSignatures.json");
else
    outFile = p.Results.OutFile;
end

srcRoot = fullfile(repoRoot, "src", "+nrRadar");
api = nrRadarTest.scanPublicApiSignatures(srcRoot);

meta = struct();
meta.generatedAt = string(datetime("now"));
meta.matlabVersion = string(version);

payload = struct("meta", meta, "api", api);

txt = jsonencode(payload);
fid = fopen(outFile, "w");
fwrite(fid, txt, "char");
fclose(fid);

fprintf("[generatePublicApiBaseline] Wrote %s\n", outFile);
end

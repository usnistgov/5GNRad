function makeBaselines(varargin)
% MAKEBASELINES Generate Layer C (regression) and Layer D (performance) baselines.
%
% This script runs a small set of representative scenarios with fixed
% overrides (typically 1 drop, serial execution) and saves compact baseline
% summaries under tests/baselines/.
%
% Usage:
%   cd(repoRoot); setup;
%   addpath("tests");
%   makeBaselines
%
% Optional name-value:
%   'Scenarios' : string array of scenario folder names under /examples
%   'OutDir'    : output folder (default: tests/baselines)
%
% Notes:
% - Baselines should be regenerated intentionally (e.g., when changing
%   algorithms). Commit updated baseline .mat files along with a short note.
% - Baselines store *summaries*, not large cubes, to reduce brittleness.
%
% 2026 NIST/CTL
%
% This file is available under the terms of the NIST License.

p = inputParser;
addParameter(p, 'Scenarios', ["uma_trp1_fulldigital_rx_singletarget","uma_trp13_hybrid"]);
addParameter(p, 'OutDir', "");
parse(p, varargin{:});
scenarios = string(p.Results.Scenarios);

repoRoot = fileparts(fileparts(mfilename('fullpath')));
if strlength(p.Results.OutDir) == 0
    outDir = fullfile(repoRoot, "tests", "baselines");
else
    outDir = p.Results.OutDir;
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

setupPath = fullfile(repoRoot, "setup.m");
if exist(setupPath, 'file') == 2
    run(setupPath);
end

rng(0, "twister");

baselineIndex = struct();
baselineIndex.generatedAt = string(datetime("now"));
baselineIndex.matlabVersion = string(version);
baselineIndex.scenarios = scenarios;

for i = 1:numel(scenarios)
    scenarioName = scenarios(i);
    scenarioPath = fullfile("examples", scenarioName);

    fprintf("[makeBaselines] Scenario %d/%d: %s\n", i, numel(scenarios), scenarioName);

    overrides = {"simulation.nMaxDrop", 1};

    [simConfig, stConfig, prsConfig, geometry, sensConfig, bg, tgtCh] = ...
        nrRadar.cfg.configureScenario(scenarioPath, "Overrides", overrides);

    % Determinism: serial
    runFcn = @() nrRadar.run(simConfig, stConfig, prsConfig, geometry, sensConfig, bg, tgtCh, [], "Parallel", "off");

    % Warmup (JIT, caching)
    try
        runFcn();
    catch ME
        warning("Warmup failed for %s: %s", scenarioName, ME.message);
    end

    tStart = tic;
    [results, detStats, detectionOutput, info] = runFcn();
    elapsed_s = toc(tStart);

    baseline = struct();
    baseline.meta = struct();
    baseline.meta.scenarioName = scenarioName;
    baseline.meta.generatedAt  = string(datetime("now"));
    baseline.meta.matlabVersion = string(version);

    baseline.regression = summarizeRegression(results, detStats, detectionOutput, info);
    baseline.performance = summarizePerformance(results, detStats, detectionOutput, info, elapsed_s);

    outFile = fullfile(outDir, scenarioName + "_drop1.mat");
    save(outFile, "baseline");
    fprintf("[makeBaselines] Wrote %s\n", outFile);

    baselineIndex.files(i) = string(outFile);
end

save(fullfile(outDir, "baselineIndex.mat"), "baselineIndex");
fprintf("[makeBaselines] Wrote baseline index.\n");

end

function reg = summarizeRegression(results, detStats, detectionOutput, info)
% Compact, stable summary for backward-compat checks.

reg = struct();

% DetStats (use first drop if vectorized)
try
    reg.detStats = struct( ...
        "TP", detStats.truePositive(1), ...
        "FN", detStats.falseNegative(1), ...
        "FP", detStats.falsePositve(1), ...
        "FPR", detStats.falseAlarmProb(1));
catch
    reg.detStats = detStats; %#ok<STRNU>
end

% First payload detection (schema-safe)
try
    payload = detectionOutput(1).detectionPayload;
catch
    payload = detectionOutput(1);
end

reg.nDetections = 0;
reg.firstEstimate = [];
reg.firstVel = [];

if isfield(payload, "stPositionEstimate") && ~isempty(payload.stPositionEstimate)
    reg.nDetections = size(payload.stPositionEstimate, 1);
    reg.firstEstimate = payload.stPositionEstimate(1,:);
end
if isfield(payload, "velParallel") && ~isempty(payload.velParallel)
    reg.firstVel = payload.velParallel(1,:);
end

% Common scalar error stats if present
fields = ["rangeError","velocityError","positionErrorH","positionErrorV","snr"];
for k = 1:numel(fields)
    f = fields(k);
    if isfield(results, f) && ~isempty(results.(f))
        x = results.(f);
        reg.(f) = struct("mean", mean(x(:), "omitnan"), "median", median(x(:), "omitnan"));
    end
end

% Info light-touch (avoid large fields)
if nargin >= 4 && isstruct(info)
    reg.info = struct();
    if isfield(info,"nDrop"), reg.info.nDrop = info.nDrop; end
    if isfield(info,"cpi"), reg.info.cpi = info.cpi; end
end

end

function perf = summarizePerformance(results, detStats, detectionOutput, info, elapsed_s)
% Stable performance summary that avoids flakiness:
% - elapsed time for a single warmed run
% - output sizes (bytes) as a proxy for memory blow-ups

perf = struct();
perf.elapsed_s = elapsed_s;

% Bytes of selected outputs (proxy for memory regressions)
perf.bytes = struct();
perf.bytes.results = localBytes(results);
perf.bytes.detStats = localBytes(detStats);
perf.bytes.detectionOutput = localBytes(detectionOutput);

% Also store a couple of key array sizes if available
perf.sizes = struct();
if isfield(results,"rdmPeak")
    perf.sizes.rdmPeak = size(results.rdmPeak);
end

end

function b = localBytes(x)
try
    s = whos("x");
    b = s.bytes;
catch
    b = NaN;
end
end

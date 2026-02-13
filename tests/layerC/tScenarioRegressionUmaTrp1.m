classdef tScenarioRegressionUmaTrp1 < nrRadarTest.BaseTestCase
% TSCENARIOREGRESSIONUMATRP1
% Layer C regression test: runs a small reference scenario and compares a
% stable summary of outputs against a stored baseline (.mat).
%
% Baseline generation:
%   setup; addpath("tests"); makeBaselines
%
% Tags:
%   layerC, regression

    properties (Constant)
        ScenarioName = "uma_trp1_fulldigital_rx_singletarget";
    end

    methods (TestClassSetup)
        function setupOnce(tc)
            % Ensure project paths are configured
            % setup();
            rng(0,"twister");
        end
    end

    methods (Test, TestTags={'regression'})
        function compareAgainstBaseline(tc)
            repoRoot     = fileparts(tc.RepoRoot);
            scenarioPath = fullfile(repoRoot, "examples", tc.ScenarioName);
            baselineFile = fullfile(repoRoot, "tests", "baselines", tc.ScenarioName + "_drop1.mat");

            % Skip cleanly if baseline not present on this runner
            tc.assumeTrue(isfile(baselineFile), ...
                "Baseline missing. Run: setup; addpath('tests'); makeBaselines");

            tc.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture( ...
             "nrRadar:IO:TargetChannelMissing"));

            % Skip if required toolboxes/licenses are not available
            nrRadarTest.assumeToolboxes(tc);

            % Optional: suppress known benign warning (recommended once you add a warning ID)
            % tc.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture("nrRadar:IO:TargetChannelMissing"));

            S = load(baselineFile, "baseline");
            baseline = S.baseline;

            overrides = {"simulation.nMaxDrop", 1};
            [simConfig, stConfig, prsConfig, geometry, sensConfig, bg, tgtCh] = ...
                nrRadar.cfg.configureScenario(scenarioPath, "Overrides", overrides);

            % Run (serial for determinism)
            [results, detStats, detectionOutput, info] = ...
                nrRadar.run(simConfig, stConfig, prsConfig, geometry, sensConfig, bg, tgtCh, [], "Parallel","off");

            got = localSummarizeRegression(results, detStats, detectionOutput, info);

            % Baseline format: baseline.regression.<fields>
            tc.verifyTrue(isfield(baseline,"regression"), "Baseline missing 'regression' summary field.");
            base = baseline.regression;

            % ---- Stable comparisons (tolerant) ----
            tc.verifyEqual(got.detStats, base.detStats);
            tc.verifyEqual(got.nDetections, base.nDetections);

            tc.verifyEqual(got.firstEstimate, base.firstEstimate, ...
                "AbsTol", 1e-6, "RelTol", 1e-6);
            tc.verifyEqual(got.firstVel, base.firstVel, ...
                "AbsTol", 1e-6, "RelTol", 1e-6);

        end
    end
end

function reg = localSummarizeRegression(results, detStats, detectionOutput, info)
%LOCALSUMMARIZEREGRESSION Create a stable, compact summary for regression checks.

reg = struct();

% Detection stats (keep as small struct)
try
    reg.detStats = struct( ...
        "TP", detStats.truePositive(1), ...
        "FN", detStats.falseNegative(1), ...
        "FP", detStats.falsePositve(1), ...
        "FPR", detStats.falseAlarmProb(1));
catch
    % Fall back if fields are different
    reg.detStats = detStats;
end

% Normalize payload shape (supports both "payload" and direct struct forms)
try
    payload = detectionOutput(1).detectionPayload;
catch
    payload = detectionOutput(1);
end

reg.nDetections  = 0;
reg.firstEstimate = [];
reg.firstVel      = [];

if isfield(payload,"stPositionEstimate") && ~isempty(payload.stPositionEstimate)
    reg.nDetections  = size(payload.stPositionEstimate,1);
    reg.firstEstimate = payload.stPositionEstimate(1,:);
end
if isfield(payload,"velParallel") && ~isempty(payload.velParallel)
    reg.firstVel = payload.velParallel(1,:);
end

% Stable error/SNR summaries (mean/median)
fields = ["rangeError","velocityError","positionErrorH","positionErrorV","snr"];
for k = 1:numel(fields)
    f = fields(k);
    if isfield(results,f) && ~isempty(results.(f))
        x = results.(f);
        reg.(f) = struct( ...
            "mean",   mean(x(:),"omitnan"), ...
            "median", median(x(:),"omitnan"));
    end
end

% Minimal info snapshot (optional)
reg.info = struct();
if isstruct(info)
    if isfield(info,"nDrop"), reg.info.nDrop = info.nDrop; end
    if isfield(info,"cpi"),   reg.info.cpi   = info.cpi;   end
end
end

function localCompareOptionalStat(tc, got, base, fieldName, tol)
%LOCALCOMPAREOPTIONALSTAT Compare optional stats if present in both structs.
if isfield(got, fieldName) && isfield(base, fieldName)
    tc.verifyEqual(got.(fieldName).mean, base.(fieldName).mean, "AbsTol", tol, "RelTol", tol);
    tc.verifyEqual(got.(fieldName).median, base.(fieldName).median, "AbsTol", tol, "RelTol", tol);
end
end

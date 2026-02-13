classdef tPerformanceBudget < nrRadarTest.BaseTestCase
    % Layer D: performance regression test (non-flaky).
    %
    % Strategy:
    % - Use warmed single-run timing (not timeit) with generous tolerance.
    % - Compare output byte sizes to catch memory blow-ups deterministically.
    % - Runs only when a performance baseline exists.
    %
    % Baseline generation:
    %   makeBaselines
    %
    methods (Test, TestTags=["performance","layerD"])
        function runtimeAndMemoryWithinBudget(tc)
            repoRoot = fileparts(tc.RepoRoot);

            baselineIndexFile = fullfile(repoRoot, "tests", "baselines", "baselineIndex.mat");
            tc.assumeTrue(isfile(baselineIndexFile), ...
                "Missing baselines. Run tests/makeBaselines.m and commit tests/baselines outputs.");

            S = load(baselineIndexFile, "baselineIndex"); %#ok<NASGU>
            scenarios = ["uma_trp1_fulldigital_rx_singletarget","uma_trp13_hybrid"];

            for i = 1:numel(scenarios)
                scenarioName = scenarios(i);
                baselineFile = fullfile(repoRoot, "tests", "baselines", scenarioName + "_drop1.mat");
                tc.assumeTrue(isfile(baselineFile), "Missing baseline for " + scenarioName);

                B = load(baselineFile, "baseline");
                basePerf = B.baseline.performance;

                scenarioPath = fullfile("..\..\examples", scenarioName);
                overrides = {"simulation.nMaxDrop", 1};
                tc.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture( ...
    "nrRadar:IO:TargetChannelMissing"));

                [sim, st, prs, geom, sens, bg, tgtCh] = nrRadar.cfg.configureScenario(scenarioPath, "Overrides", overrides);

                runFcn = @() nrRadar.run(sim, st, prs, geom, sens, bg, tgtCh, [], "Parallel", "off");

                % Warmup
                try, runFcn(); catch, end

                t0 = tic;
                [results, detStats, detectionOutput] = runFcn();
                elapsed = toc(t0);

                bytes.results = localBytes(results);
                bytes.detStats = localBytes(detStats);
                bytes.detectionOutput = localBytes(detectionOutput);

                % --- Budget rules (generous to avoid flakiness) ---
                % Runtime must be within 2.5x baseline (single run timing can vary)
                tc.verifyLessThanOrEqual(elapsed, 2.5 * basePerf.elapsed_s, ...
                    sprintf("Runtime regression in %s: %.3fs vs baseline %.3fs", scenarioName, elapsed, basePerf.elapsed_s));

                % Output bytes should not exceed 1.5x baseline (deterministic proxy)
                tc.verifyLessThanOrEqual(bytes.results, 1.5 * basePerf.bytes.results + 1e6);
                tc.verifyLessThanOrEqual(bytes.detStats,  basePerf.bytes.detStats );
                tc.verifyLessThanOrEqual(bytes.detectionOutput, 1.5 * basePerf.bytes.detectionOutput + 1e6);
            end
        end
    end
end

function b = localBytes(x)
s = whos("x");
b = s.bytes;
end

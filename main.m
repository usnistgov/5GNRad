%%  Run PRS-based radar simulations across multiple 3GPP scenarios
%   MAIN executes the end-to-end 5G NR PRS radar simulation
%   pipeline for a list of scenario folders and exports results for each
%   scenario. The script loads scenario-specific configuration files, runs
%   the simulation (optionally using parallel workers), and saves outputs to
%   disk.
%
%   Workflow:
%     1) SETUP initializes required paths/toolboxes and project settings.
%     2) CONFIGSCENARIO loads:
%        - simConfig: simulation/system parameters (carrier, bandwidth, etc.)
%        - stConfig : ground-truth target states
%        - prsConfig: PRS configuration
%        - geometry : TX/RX/TRP geometry for the scenario
%        - sensConfig: sensing / processing configuration
%        - backgroundChannel: background/environment channel model
%        - targetChannel    : target channel model
%     3) RUN5GNRAD executes the PRS radar processing chain and returns:
%        - results, detStats, detectionOutput
%     4) EXPORTRESULTS writes scenario results to disk.
%
%   Configuration parameters (edit in the script):
%     scenarios       - String array of scenario folder paths.
%     desiredWorkers  - Number of parallel workers requested by RUN5GNRAD.
%     parallelMode    - "on"/"off" toggle for parallel execution.
%
%
%   2025-2026 NIST/CTL Steve Blandino
%
%   This file is available under the terms of the NIST License.

scenarios = [
    "examples/uma_trp1_3gpp"
    "examples/uma_trp9_3gpp_08lambda"
    "examples/uma_trp13_hybrid"
    "examples/uma_trp20_3gpp"    
    "examples/uma_trp1_fulldigital_rx_singletarget"    
    ];


parallelMode   = "on";
desiredWorkers = 10;        % choose based on RAM

%% Set path
setup();

for i = 1:numel(scenarios)
    scenarioPath = scenarios(i);
    fprintf("[%d/%d] Scenario: %s\n", i, numel(scenarios), scenarioPath);

    try
        %% Load configs
        [simConfig, stConfig, prsConfig, geometry, ...
            sensConfig,backgroundChannel,targetChannel] = ...
            nrRadar.cfg.configureScenario(scenarioPath);

        %% Run
        [results, detStats, detectionOutput] = ...
            nrRadar.run(simConfig, stConfig, prsConfig, geometry, ...
            sensConfig,backgroundChannel,targetChannel,desiredWorkers, 'Parallel', parallelMode);

        %% Store Results
        nrRadar.io.exportResults(scenarioPath, results,detStats,detectionOutput)
    catch ME
        fprintf(2, "Simulation Failed");
        fprintf(2, "    %s\n\n", ME.getReport('extended','hyperlinks','off'));
    end
end

fprintf('Simulation Complete\n');

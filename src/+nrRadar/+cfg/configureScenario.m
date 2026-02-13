function [simulation, target, prs,geometry,sens,backgroundChannel,targetChannel] =  configureScenario(scenarioPath, varargin)
% CONFIGSCENARIO Load all configuration components for a 5GNR radar scenario
%   [SIMULATION, TARGET, PRS, GEOMETRY, SENS, BACKGROUNDCHANNEL] = CONFIGSCENARIO(SCENARIOPATH)
%   loads all configuration and input data required to simulate a PRS-based
%   radar scenario. The function aggregates configuration from a set of
%   supporting files located within the specified scenario path.
%
%   Inputs:
%     scenarioPath - Path to the scenario folder (must contain /Input/)
%
%   Outputs:
%     SIMULATION        - Structure with system-level simulation parameters
%     TARGET            - Structure with target positions and velocities
%     PRS               - Structure with PRS configuration parameters
%     GEOMETRY          - Structure with TX and RX geometry
%     SENS              - Structure with sensing parameters (window, FFT sizes, etc.)
%     BACKGROUNDCHANNEL - Structure or cell array with background channel information
%
%   Each configuration is loaded using a dedicated helper function:
%     - configSimulation.m
%     - loadTarget.m
%     - configPrs.m
%     - loadGeometry.m
%     - loadSensConfig.m
%     - loadBackgroundChannel.m
%
%   Example:
%     [sim, tgt, prs, geom, sens, bg] = configScenario('examples/UMi-Av25');

%   2025 NIST/CTL Steve Blandino

%   This file is available under the terms of the NIST License.

import nrRadar.io.*
import nrRadar.cfg.*

p = inputParser;
addParameter(p, 'Overrides', []);
parse(p, varargin{:});
overrides = p.Results.Overrides;

simulation = configureSimulation(scenarioPath);
prs = configurePrs(scenarioPath);
sens = configureSensing(scenarioPath);
target = loadTarget(scenarioPath);
geometry = loadGeometry(scenarioPath);
backgroundChannel = loadBackgroundChannel(scenarioPath);
targetChannel = loadTargetChannel(scenarioPath);

% Optional: apply runtime overrides (useful for sweeps/experiments).
% Overrides can be provided as either:
%   - N-by-2 cell array: { 'sens.cfarThreshold', 10; 'simulation.nMaxDrop', 20; ... }
%   - name/value cell array: { 'sens.cfarThreshold', 10, 'simulation.nMaxDrop', 20, ... }
%   - struct with top-level fields: simulation/target/prs/geometry/sens (or simConfig/stConfig/prsConfig/sensConfig)
if ~isempty(overrides)
    [simulation, target, prs, geometry, sens] = nrRadar.util.applyOverrides(...
        simulation, target, prs, geometry, sens, overrides);
end
fprintf('Configuration Complete\n');

end


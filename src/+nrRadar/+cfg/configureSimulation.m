function [params] =  configureSimulation(scenarioPath, varargin)
% CONFIGSIMULATION Load system-level simulation parameters for 5GNR radar
%   PARAMS = CONFIGSIMULATION(SCENARIOPATH) loads the simulation configuration
%   for the scenario located at SCENARIOPATH. The function reads the
%   tab-delimited configuration file:
%       SCENARIOPATH/Input/simulationConfig.txt
%   and imports the antenna models:
%       SCENARIOPATH/Input/txAntennaModel.json
%       SCENARIOPATH/Input/rxAntennaModel.json
%
%   The configuration file must contain two tab-delimited columns:
%     (1) parameter name (string)
%     (2) parameter value (string)
%   The parameters are converted/validated using FIELDTONUM and default
%   values are applied when missing or invalid.
%
%
%   Output
%     PARAMS - Structure containing validated simulation parameters and
%              imported antenna models. Common fields include:
%       * systemNF                    - Noise figure [dB] (default: 7)
%       * systemBw                    - System bandwidth [Hz] (default: 100e6)
%       * channelScenario             - Channel scenario string:
%                                      {'UMiAV','UMaAV','RMaAV'} (default: 'UMiAV')
%       * antennaCouplingEfficiency   - Coupling efficiency in [0,1] (default: 1)
%       * carrierSubcarrierSpacing    - Subcarrier spacing [kHz] (default: 120)
%       * carrierNSizeGrid            - NR grid size (default: 66)
%       * nStDrop                     - Number of states per drop (default: 1)
%       * maxRangeInterest            - Max range of interest [m] (default: 400)
%       * trpYawDeg                   - TRP yaw rotation [deg] (default: 0)
%       * txPower                     - Transmit power [dBm] (default: 52)
%       * txAntenna                   - Imported TX antenna/beamforming model
%       * rxAntenna                   - Imported RX antenna/beamforming model
%       * systemFc                    - Carrier frequency [Hz]. If not provided,
%                                      defaults to PARAMS.txAntenna.meta.fc_Hz
%
%
%   See also FIELDTONUM, IMPORTBFFROMJSON.
%
%   2026 NIST/CTL Steve Blandino
%
%   This file is available under the terms of the NIST License.

import nrRadar.util.fieldToNum

cfgPath = fullfile(scenarioPath, 'Input/simulationConfig.txt');
txAntennaPath = fullfile(scenarioPath, 'Input/txAntennaModel.json');
rxAntennaPath = fullfile(scenarioPath, 'Input/rxAntennaModel.json');

requiredFiles = { ...
    'simulationConfig.txt', cfgPath; ...
    'txAntennaModel.json',  txAntennaPath; ...
    'rxAntennaModel.json',  rxAntennaPath ...
    };

missing = requiredFiles(~cellfun(@isfile, requiredFiles(:,2)), :);
if ~isempty(missing)
    msg = sprintf( ...
        "Simulation configuration is incomplete for scenarioPath:\n  %s\n\n" + ...
        "Missing required file(s):\n%s\n\n" + ...
        "Expected under:\n  %s\n" + ...
        "Fix: ensure the 'Input' folder contains these files.", ...
        string(scenarioPath), ...
        join("  - " + string(missing(:,1)) + " (looked for: " + string(missing(:,2)) + ")", newline), ...
        string(fullfile(scenarioPath,'Input')));

    ME = MException('NIST5GNRad:MissingFile', '%s', msg);
    throwAsCaller(ME);
end

fprintf('Loading Simulation Configuration\n');
paramsList = readtable(cfgPath,'Delimiter','\t', 'Format','%s %s' );
paramsCell = (table2cell(paramsList))';
params = cell2struct(paramsCell(2,:), paramsCell(1,:), 2);
%% Check validity
params = fieldToNum(params, 'systemNF', [0 20], 'step', eps, 'defaultValue', 7);
params = fieldToNum(params, 'systemBw', [10 20 40 50 60 80 100 200 400]*1e6, 'defaultValue', 100e6);
params = fieldToNum(params, 'channelScenario', {'UMiAV', 'UMaAV', 'RMaAV'}, 'defaultValue', 'UMiAV');
params = fieldToNum(params, 'antennaCouplingEfficiency', [0 1], 'step', eps, 'defaultValue', 1);
params = fieldToNum(params, 'carrierSubcarrierSpacing', 15*2.^(0:4), 'defaultValue', 120);
params = fieldToNum(params, 'carrierNSizeGrid', [1,275], 'step', 1, 'defaultValue', 66);
params = fieldToNum(params, 'nStDrop', [0,10], 'step', 1, 'defaultValue', 1);
params = fieldToNum(params, 'maxRangeInterest', [10,2000], 'step', eps, 'defaultValue', 400);
params = fieldToNum(params, 'trpYawDeg', [0,360], 'step', eps, 'defaultValue', 0);
params = fieldToNum(params, 'txPower', [-300 300], 'step', eps, 'defaultValue', 52);
params = fieldToNum(params, 'nMaxDrop', [1 inf], 'step', eps, 'defaultValue', 50);
params.txAntenna = nrRadar.io.importBFfromJSON(txAntennaPath);
params.rxAntenna = nrRadar.io.importBFfromJSON(rxAntennaPath);
params = fieldToNum(params, 'systemFc', [1e9 60e9], 'step', eps, 'defaultValue', params.txAntenna.meta.fc_Hz);


end

function [params] =  configSimulation(scenarioPath, varargin)
% CONFIGSIMULATION Load system-level simulation parameters for 5GNR radar
%   PARAMS = CONFIGSIMULATION(SCENARIOPATH) loads the simulation configuration 
%   from the file 'Input/simulationConfig.txt' under the specified SCENARIOPATH. 
%   It parses name-value pairs, validates them against defined ranges or allowed 
%   sets, and assigns default values if missing or invalid.
%
%   The configuration file must be tab-delimited with two columns: parameter name and value.
%
%   Output:
%     PARAMS - Structure containing validated simulation parameters:
%         * systemFc                  - Carrier frequency [Hz] (default: 30e9)
%         * systemNF                  - Noise figure [dB] (default: 7)
%         * systemBw                  - System bandwidth [Hz], allowed: 20–400 MHz (default: 100e6)
%         * channelScenario           - Channel scenario string: {'UMiAV','UMaAV','RMaAV'} (default: 'UMiAV')
%         * antennaNumH               - Number of antenna elements (horizontal) (default: 32)
%         * antennaNumV               - Number of antenna elements (vertical) (default: 32)
%         * antennaCouplingEfficiency - TX coupling efficiency (0–1) (default: 0.7)
%         * carrierSubcarrierSpacing  - Subcarrier spacing in kHz, values: {15, 30, 60, 120, 240} (default: 120)
%         * carrierNSizeGrid          - Number of resource blocks (1–275) (default: 66)
%
%   If the file is not found, PARAMS is returned as an empty struct.
%
%   Example:
%     sim = configSimulation('examples/UMi-Av25');
%
%   See also: FIELDTONUM

%   2025 NIST/CTL Steve Blandino
%   This file is available under the terms of the NIST License.

cfgPath = fullfile(scenarioPath, 'Input/simulationConfig.txt');

if isfile(cfgPath)

    paramsList = readtable(cfgPath,'Delimiter','\t', 'Format','%s %s' );
    paramsCell = (table2cell(paramsList))';
    params = cell2struct(paramsCell(2,:), paramsCell(1,:), 2);
    %% Check validity
    params = fieldToNum(params, 'systemFc', [1e9 60e9], 'step', eps, 'defaultValue', 30e9);
    params = fieldToNum(params, 'systemNF', [0 20], 'step', eps, 'defaultValue', 7);
    params = fieldToNum(params, 'systemBw', [10 20 40 50 60 80 100 200 400]*1e6, 'defaultValue', 100e6);
    params = fieldToNum(params, 'channelScenario', {'UMiAV', 'UMaAV', 'RMaAV'}, 'defaultValue', 'UMiAV');
    params = fieldToNum(params, 'antennaNumH', [1 1024], 'step', 1, 'defaultValue', 32);
    params = fieldToNum(params, 'antennaNumV', [1 1024], 'step', eps, 'defaultValue', 32);
    params = fieldToNum(params, 'antennaCouplingEfficiency', [0 1], 'step', eps, 'defaultValue', 0.7);
    params = fieldToNum(params, 'carrierSubcarrierSpacing', 15*2.^(0:4), 'defaultValue', 120);
    params = fieldToNum(params, 'carrierNSizeGrid', [1,275], 'step', 1, 'defaultValue', 66);

else
    params = [];
end
end

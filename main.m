%% 5GNR Radar (5NRad)
% This script simulates a radar system based on 5G New Radio (NR) using
% Positioning Reference Signals (PRS). It performs the following:
%
%   - Loads a predefined scenario configuration
%   - Generates and transmits PRS OFDM waveforms
%   - Simulates monostatic sensing through a CDL channel model
%   - Processes the received signal to compute range-Doppler maps
%   - Estimates target position and velocity over multiple time instances
%   - Saves results to /Output/error.csv within the scenario folder
%
%
% Example:
%   Run the script from the project root directory with the following:
%       scenarioNameStr = 'examples/UMi-Av25';
%
%   Output:
%       Results table written to:
%           ./examples/UMi-Av25/Output/error.csv

%   2025 NIST/CTL Steve Blandino
%
%   This file is available under the terms of the NIST License.



% scenarioNameStrVector = {'examples3GPP/UMa-Av200-8x8-120',...
%     'examples3GPP/UMa-Av200-8x8-30', ...
%     'examples3GPP/UMa-Av200-8x8-30-24RB', 'examples3GPP/UMa-Av200-8x8-30-51RB',...
%     'examples3GPP/UMa-Av200-8x8-30-133RB','examples3GPP/UMa-Av200-8x8-60', ...
%     'examples3GPP/UMa-Av200-8x8-120', 'examples3GPP/UMa-Av200-8x8-128ss', ...
%     'examples3GPP/UMa-Av200-8x8-192ss'};%, 'examples/UMa-Av200-8x4', 'examples/UMa-Av200-8x8','examples/UMa-Av200-16x8','examples/UMa-Av200-16x16','examples/UMa-Av200-24x24'};
scenarioNameStrVector = {'examples/UMi-Av25'}
%% Set path
rootFolderPath = pwd;
fprintf('--------5G NR Radar --------\n');
fprintf('Current root folder:\n\t%s\n',rootFolderPath);
[path,folderName] = fileparts(rootFolderPath);
addpath(genpath(fullfile(rootFolderPath,'src')));

for i = 1:length(scenarioNameStrVector)
    scenarioNameStr = scenarioNameStrVector{i};
    fprintf('Use customized scenario: %s.\n',scenarioNameStr);

    % %% Run
    [simConfig, stConfig, prsConfig, geometry, sensConfig,backgroundChannel,targetChannel] = configScenario(scenarioNameStr);
    results = run5GNRad(simConfig, stConfig, prsConfig, geometry, sensConfig,backgroundChannel,targetChannel);

    %% Store Results
    resultsTab = struct2table(results);
    outputPath = fullfile(scenarioNameStr, 'Output');

    if ~isfolder(outputPath)
        mkdir(outputPath);
    else
        rmdir(outputPath, 's');
        mkdir(outputPath);
    end
    writetable(resultsTab, fullfile(outputPath, 'error.csv'))

end
function params = configSens(scenarioPath)
% CONFIGSENS Load and validate sensing configuration for 5GNR radar
%   PARAMS = CONFIGSENS(SCENARIOPATH) reads sensing parameters from the
%   file 'Input/sensConfig.txt' within the specified SCENARIOPATH directory.
%   It converts values to the correct types, validates them against expected
%   ranges, and applies defaults where necessary.
%
%   The configuration file must be tab-delimited, with two columns: parameter
%   name and value.
%
%   Output:
%     PARAMS - Structure with the following validated fields:
%       * dopplerFftLen        - Doppler FFT size (default: 64)
%       * window               - Doppler window type: 'rect', 'hamming',
%                                'blackmanharris', or 'gaussian' (default: 'blackmanharris')
%       * windowLen            - Window length (default: dopplerFftLen)
%       * windowOverlap        - Overlap ratio [0, 1) (default: 0.5)
%       * numberSensingSymbols - OFDM symbols per CPI (default: 256)
%       * cfarGrdCellRange     - CFAR guard cells in range (default: 0)
%       * cfarGrdCellVelocity  - CFAR guard cells in Doppler (default: 0)
%       * cfarTrnCellRange     - CFAR training cells in range (default: 0)
%       * cfarTrnCellVelocity  - CFAR training cells in Doppler (default: 0)
%       * cfarThreshold        - CFAR detection threshold (default: 3)
%       * isCfar               - Flag indicating CFAR configuration present
%
%   If the file is missing, PARAMS is returned as an empty struct.
%
%   Example:
%       sens = CONFIGSENS('examples/UMi-Av25');
%
%   See also: FIELDTONUM

%   2025 NIST/CTL Steve Blandino
%   This file is available under the terms of the NIST License.


%% Load params
cfgPath = fullfile(scenarioPath, 'Input/sensConfig.txt');
if isfile(cfgPath)
    paramsList = readtable(cfgPath,'Delimiter','\t', 'Format','%s %s' );
    paramsCell = (table2cell(paramsList))';
    params = cell2struct(paramsCell(2,:), paramsCell(1,:), 2);

    %% Check validity
    params = fieldToNum(params, 'dopplerFftLen', [1 4096], 'step', eps, 'defaultValue',64);
    params = fieldToNum(params, 'window', {'rect' ,'hamming','blackmanharris','gaussian'}, 'defaultValue', 'blackmanharris');
    params = fieldToNum(params, 'windowLen', [1 params.dopplerFftLen], 'step', 1, 'defaultValue', params.dopplerFftLen);
    params = fieldToNum(params, 'windowOverlap', [0 1-eps], 'step', eps, 'defaultValue', 0.5);
    params = fieldToNum(params, 'numberSensingSymbols', [1 4096], 'step', 1, 'defaultValue',256);
    params.isCfar = any(cellfun(@(x) startsWith(x,'cfar'), fieldnames(params)));
    params = fieldToNum(params, 'cfarGrdCellRange', [0 1e4], 'step', 1, 'defaultValue', 0);
    params = fieldToNum(params, 'cfarGrdCellVelocity', [0 1e4], 'step', 1, 'defaultValue', 0);
    params = fieldToNum(params, 'cfarTrnCellRange', [0 1e4], 'step', 1, 'defaultValue', 0);
    params = fieldToNum(params, 'cfarTrnCellVelocity', [0 1e4], 'step', 1, 'defaultValue', 0);
    params = fieldToNum(params, 'cfarThreshold', [0 1e4], 'step', eps, 'defaultValue', 3);  
    params = fieldToNum(params, 'angleEstimationMethod', {'ideal' ,'nearest','scan'}, 'defaultValue', 'ideal');


else
    params = [];
end




function params = configPrs(scenarioPath)
%CONFIGPRS Load PRS configuration parameters from file.
%   PARAMS = CONFIGPRS(SCENARIOPATH) loads a PRS configuration from a 
%   'prsConfig.txt' file located under the Input folder of the specified 
%   SCENARIOPATH. The function reads the configuration parameters, converts 
%   them to the correct types, checks if they are within the expected ranges, 
%   and assigns default values where necessary.
%
%   The loaded PRS parameters include:
%
%     - PRSResourceSetPeriod  : PRS resource set slot periodicity and slot offset 
%                               specified as [TPRSPeriod TPRSOffset] in slots.
%     - PRSResourceOffset     : Slot offset of each PRS resource relative to 
%                               PRS resource set slot offset (0...511).
%     - PRSResourceRepetition : PRS resource repetition factor (1, 2, 4, 6, 8, 16, 32).
%     - PRSResourceTimeGap    : Slot offset between consecutive repeated instances 
%                               of a PRS resource (1, 2, 4, 8, 16, 32).
%     - NumRB                 : Number of PRBs allocated for PRS (0...275).
%     - RBOffset              : Starting PRB index relative to carrier resource grid (0...274).
%     - CombSize              : Comb size of PRS resources (2, 4, 6, 12).
%     - REOffset              : Starting resource element offset in the first OFDM symbol (0...(CombSize-1)).
%     - NPRSID                : Sequence identity of PRS resources (0...4095).
%     - NumPRSSymbols         : Number of OFDM symbols allocated per PRS resource (0...12).
%     - SymbolStart           : Starting OFDM symbol index in a slot for PRS resource (0...13).
%
%   If the configuration file is missing, an empty struct is returned.
%
%   See also nrPRSConfig, nrPRS, nrPRSIndices.
%
%   2025 NIST/CTL Steve Blandino
%   This file is available under the terms of the NIST License.



%% Load params
cfgPath = fullfile(scenarioPath, 'Input/prsConfig.txt');

if isfile(cfgPath)
    paramsList = readtable(cfgPath,'Delimiter','\t', 'Format','%s %s' );
    paramsCell = (table2cell(paramsList))';
    params = cell2struct(paramsCell(2,:), paramsCell(1,:), 2);

    %% Check validity
    params = fieldToNum(params, 'PRSResourceSetPeriod', [0 1e3], 'step', 1, 'defaultValue', [1 0]);
    params = fieldToNum(params, 'PRSResourceOffset', [0 511], 'step', 1, 'defaultValue', 0);
    params = fieldToNum(params, 'PRSResourceRepetition', 2.^(0:5), 'step', 1, 'defaultValue', 1);
    params = fieldToNum(params, 'PRSResourceTimeGap', 2.^(0:5), 'step', 1, 'defaultValue', 1);
    params = fieldToNum(params, 'NumRB', [0 275], 'step', 1, 'defaultValue', 66);
    params = fieldToNum(params, 'RBOffset', [0 274], 'step', 1, 'defaultValue', 0);
    params = fieldToNum(params, 'CombSize', [2 4 6 12], 'step', 1, 'defaultValue', 4);
    params = fieldToNum(params, 'REOffset', [0 params.CombSize-1], 'step', 1, 'defaultValue', 0);
    params = fieldToNum(params, 'NPRSID', [0 4095], 'step', 1, 'defaultValue', 4);
    params = fieldToNum(params, 'NumPRSSymbols', [0 12], 'step', 1, 'defaultValue', 4);
    params = fieldToNum(params, 'SymbolStart', [0 13], 'step', 1, 'defaultValue', 1);

else
    params = [];
end
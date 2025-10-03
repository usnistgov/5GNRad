function backgroundData = loadBackgroundChannel(scenarioPath, fc)
% LOADBACKGROUNDCHANNEL Load background propagation channel data
%   BKGDATA = LOADBACKGROUNDCHANNEL(SCENARIOPATH) loads background
%   multipath components from a JSON file based on the specified scenario
%   path SCENARIOPATH. The function supports scenario names containing:
%   * 'umi' → loads '30GHz_3GPP_38.901_UMa_NLOS'
%   * 'uma' → loads '30GHz_3GPP_38.901_UMi_NLOS'
%
%   The JSON file must exist in a subfolder named 'backgroundChannel' in
%   the current working directory and follow the naming convention:
%   'backgroundChannel_<scenarioType>.json'.
%
%   The function returns a structure BACKGROUNDDATA containing the decoded
%   contents of the background channel file.
%
%   Example:
%       data = loadBackgroundChannel('umi');
%
%   Errors:
%       - Raises an error if the scenario type is unrecognized or
%         unsupported.
%       - Raises an error if the corresponding JSON file is not found.

%   2025 NIST/CTL Steve Blandino

%   This file is available under the terms of the NIST License.

    % Get the current folder
    rootPath = pwd;

    % Determine scenario type (UMi, UMa, RMa)
    if contains(lower(scenarioPath), 'uma')
        scenarioType = '3GPP_38.901_UMa_NLOS';
    elseif contains(lower(scenarioPath), 'umi')
        scenarioType = '3GPP_38.901_UMi_NLOS';
    elseif contains(lower(scenarioPath), 'rma')
        scenarioType = '3GPP_38.901_RMa_NLOS';
    else
        error('Unknown scenario type in scenarioPath: %s', scenarioPath);
    end

    fcGHz = fc/1e9;
    % Build path to backgroundChannel.json
    backgroundFile = fullfile(rootPath, 'channel', ['backgroundChannel_',num2str(fcGHz),'GHz_', scenarioType, '.json']);

    % Load the JSON file
    if exist(backgroundFile, 'file')
        fid = fopen(backgroundFile, 'r');
        raw = fread(fid, inf);
        str = char(raw');
        fclose(fid);
        backgroundData = jsondecode(str);
    else
        error('Background channel file not found: %s', backgroundFile);
    end
end
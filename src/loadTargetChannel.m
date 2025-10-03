function targetData = loadTargetChannel(scenarioPath, fc)
% LOADTARGETCHANNEL Load target channel configuration file
%   TARGETDATA = LOADTARGETCHANNEL(SCENARIOPATH, FC) loads the target
%   channel data stored in a JSON file for the given scenario and
%   carrier frequency FC. The function automatically determines the
%   scenario type (UMa, UMi, or RMa) based on the SCENARIOPATH string
%   and constructs the expected JSON filename accordingly.
%
%   The JSON filename is constructed as:
%       targetChannel_<fcGHz>GHz_<scenarioType>.json
%   where:
%       <fcGHz>        = carrier frequency in GHz (derived from FC in Hz)
%       <scenarioType> = 'UMaAV', 'UMiAV', or 'RMaAV'
%
%   Example:
%       data = loadTargetChannel('path/to/uma/scenario', 28e9);
%       % Loads file: targetChannel_28GHz_UMaAV.json
%
%   Input Arguments:
%       SCENARIOPATH - Path string used to infer the scenario type.
%       FC           - Carrier frequency in Hz.
%
%   Output Arguments:
%       TARGETDATA    - Struct containing decoded JSON data if the file
%                       exists. Returns [] if the file cannot be found.
%
%   See also JSONDECODE, FULLFILE

%   2025 NIST/CTL Steve Blandino
%   This file is available under the terms of the NIST License.

% Get the current folder
rootPath = pwd;

% Determine scenario type (UMi, UMa, RMa)
if contains(lower(scenarioPath), 'uma')
    scenarioType = 'UMaAV';
elseif contains(lower(scenarioPath), 'umi')
    scenarioType = 'UMiAV';
elseif contains(lower(scenarioPath), 'rma')
    scenarioType = 'RMaAV';
else
    error('Unknown scenario type in scenarioPath: %s', scenarioPath);
end

fcGHz = fc/1e9;
% Build path to backgroundChannel.json
backgroundFile = fullfile(rootPath, 'channel', ['targetChannel_',num2str(fcGHz),'GHz_', scenarioType, '.json']);

% Load the JSON file
if exist(backgroundFile, 'file')
    fid = fopen(backgroundFile, 'r');
    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);
    targetData = jsondecode(str);
else
    targetData = [];
end
end
function targetData = loadTargetChannel(scenarioPath)
% LOADTARGETCHANNEL Load 3GPP target channel data from scenario folder
%   TARGETDATA = LOADTARGETCHANNEL(SCENARIOPATH) loads target channel
%   parameters from the JSON file 'targetChannel.json' located under:
%     fullfile(pwd, SCENARIOPATH, 'Input', 'targetChannel.json')
%
%   If the file exists, this function decodes the JSON content and returns
%   it as a MATLAB struct/array in TARGETDATA.
%
%   If the file does not exist, TARGETDATA is returned as empty and 
%   simulation will proceed with a single MPC target.
%
%   Input
%     SCENARIOPATH - String or char path to the scenario root folder.
%
%   Output
%     TARGETDATA   - Decoded JSON content,
%                    or [] if the file is not found.
%
%   2026 NIST/CTL Steve Blandino
%
%   This file is available under the terms of the NIST License.


fprintf('Loading Target Channel\n');

% Build path to targetChannel.json
channelFile = fullfile(scenarioPath, ...
    'Input','targetChannel.json');

% Load the JSON file
if exist(channelFile, 'file')
    fid = fopen(channelFile, 'r');
    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);
    targetData = jsondecode(str);
else
    targetData = [];
    warning('nrRadar:IO:TargetChannelMissing', ['3GPP Target channel file not found. ' ...
        'Proceeding simulation with single MPC target.'])
end

end
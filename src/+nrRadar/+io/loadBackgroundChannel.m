function backgroundChannel = loadBackgroundChannel(scenarioPath)
%LOADBACKGROUNDCHANNEL Load background channel data from scenario JSON file
%   BACKGROUNDCHANNEL = LOADBACKGROUNDCHANNEL(SCENARIOPATH) loads the
%   background channel structure from a JSON file located at:
%       fullfile(pwd, SCENARIOPATH, 'Input', 'backgroundChannel.json')
%
%   If the file exists, the function reads and decodes the JSON content
%   using JSONDECODE and returns it as BACKGROUNDCHANNEL. If the file does
%   not exist, BACKGROUNDCHANNEL is returned as an empty array [].
%
%   Input:
%     SCENARIOPATH - Scenario folder path (relative to the current working
%                    directory), used to locate the Input folder.
%
%   Output:
%     BACKGROUNDCHANNEL - Background channel data decoded from JSON. If the
%                        file is not found, returns [].
%
%   Example:
%     bg = loadBackgroundChannel("scenarios/UMa");
%
%   2026 NIST/CTL Steve Blandino
%
%   This file is available under the terms of the NIST License.


fprintf('Loading Background Channel\n');

channelFile = fullfile(scenarioPath, 'Input','backgroundChannel.json');

% Load the JSON file
if exist(channelFile, 'file')
    fid = fopen(channelFile, 'r');
    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);
    backgroundChannel = jsondecode(str);
else
    backgroundChannel = [];
    warning('Background channel file not found. Proceeding without background channel.')
end

end
function target = loadTarget(scenarioPath)
% LOADTARGET Load target position and velocity
%   TARGET = LOADTARGET(SCENARIOPATH) loads the target configuration from
%   the file 'Input/targetConfig.txt' located in the specified SCENARIOPATH.
%
%   The configuration file must be a numeric matrix where:
%     * Columns 1-3 specify the 3D position [x, y, z]
%     * Columns 4-6 specify the 3D velocity [vx, vy, vz]
%
%   The returned TARGET is a structure with the following fields:
%     * TARGET.position - N-by-3 matrix of positions
%     * TARGET.velocity - N-by-3 matrix of velocities
%
%   Example:
%       target = loadTarget('/path/to/scenario');

%   2025 NIST/CTL Steve Blandino

%   This file is available under the terms of the NIST License.

cfgPath = fullfile(scenarioPath, 'Input/targetConfig.txt');

targetMatrix = readmatrix(cfgPath);
target.position = targetMatrix(:,1:3);
target.velocity = targetMatrix(:,4:6);
end
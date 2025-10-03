function geometry = loadGeometry(scenarioPath)
% LOADGEOMETRY Load transmitter and receiver geometry
%   GEOMETRY = LOADGEOMETRY(SCENARIOPATH) reads the base station (BS)
%   configuration file located at 'Input/bsConfig.txt' under the specified
%   SCENARIOPATH. The function assumes that the same geometry applies to
%   both transmitter and receiver.
%
%   The returned GEOMETRY is a structure with the following fields:
%     * GEOMETRY.tx - matrix of transmitter positions
%     * GEOMETRY.rx - matrix of receiver positions (equal to tx)
%
%   Example:
%       geometry = loadGeometry('/path/to/scenario');
%
%   The file 'bsConfig.txt' must be a numeric matrix compatible with
%   MATLAB's readmatrix function.

%   2025 NIST/CTL Steve Blandino

bsFile = fullfile(scenarioPath, 'Input/bsConfig.txt');

geometry.tx = readmatrix(bsFile);
geometry.rx = geometry.tx;

end
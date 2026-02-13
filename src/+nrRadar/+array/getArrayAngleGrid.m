function [azGrid, elGrid, vis]= getArrayAngleGrid(numPtsVert, numPtsHorz, dVert, dHorz, lambda)
% GETARRAYANGLEGRID Generate azimuth/elevation grids from array geometry
%   [AZGRID, ELGRID, VIS] = GETARRAYANGLEGRID(NUMPTSVERT, NUMPTSHORZ, DVERT, DHORZ, LAMBDA)
%   builds uniformly sampled direction-cosine axes based on the array
%   element spacings DVERT (vertical) and DHORZ (horizontal), maps them to
%   elevation (ELGRID) and azimuth (AZGRID), and returns a visibility mask
%   (unit disk in the direction-cosine plane).
%
%   Inputs:
%       NUMPTSVERT  - Number of samples along the vertical (z) axis.
%       NUMPTSHORZ  - Number of samples along the horizontal (y) axis.
%       DVERT       - Vertical element spacing in meters.
%       DHORZ       - Horizontal element spacing in meters.
%       LAMBDA      - Wavelength in meters.
%
%   Outputs:
%       AZGRID  - Azimuth grid in degrees, size [NUMPTSVERT x NUMPTSHORZ].
%       ELGRID  - Elevation grid in degrees, same size as AZGRID.
%       VIS     - Logical mask equal to 1 inside the visible region where
%                 u_y^2 + u_z^2 <= 1.
%
%   Details:
%       * Direction cosines are defined as:
%             u_y = (λ / DHORZ) * f_y
%             u_z = (λ / DVERT) * f_z
%         where f_y and f_z are normalized spatial frequencies (cycles per
%         element). Here we construct centered, uniformly spaced samples
%         over [-0.5, 0.5) for (f_y, f_z) and then scale by λ/D.
%
%       * If DHORZ = DVERT = λ/2, then u_y and u_z span [-1, 1), and VIS
%         corresponds exactly to the physical visible unit disk.
%
%       * Angle mapping:
%             el = asind(u_z)
%             az = atan2d(u_y, sqrt(max(0, 1 - u_y^2 - u_z^2)))
%
%   Example:
%       lambda = 0.01;                 % ~30 GHz
%       dHorz = lambda/2; dVert = lambda/2;
%       [az, el, vis] = getArrayAngleGrid(129, 257, dVert, dHorz, lambda);
%       imagesc(az(1,:), el(:,1), vis); axis xy
%       xlabel('Azimuth (°)'); ylabel('Elevation (°)'); title('Visible region')
%
%   See also: NDGRID, ATAN2D, ASIND.
%
%   2025 NIST/CTL Steve Blandino
%   This file is available under the terms of the NIST License.

%  Centered, uniformly spaced sample indices
kz = (-numPtsVert/2:numPtsVert/2-1);                       % centered bins (fftshift)
ky = (-numPtsHorz/2:numPtsHorz/2-1);

% Normalized spatial frequency (cycles per element)
fz = kz / numPtsVert;                                     % cycles per element
fy = ky / numPtsHorz;

% Direction cosines
uz = (lambda/dVert) * fz;                                % = sin(el)          if dz=λ/2 -> 2*fz
uy = (lambda/dHorz) * fy;                                % = cos(el)sin(az)   if dy=λ/2 -> 2*fy

% Grids: rows = z (vertical), cols = y (horizontal)
[uzGrid, uyGrid] = ndgrid(uz, uy);

% Visibility mask (unit disk in (uy,uz))
vis = (uzGrid.^2 + uyGrid.^2) <= 1 + eps;

% Map to angles 
elGrid = asind( max(-1,min(1, uzGrid)) );             % [-90, +90] deg
azGrid = atan2d( uyGrid, sqrt(max(0, 1 - uzGrid.^2 - uyGrid.^2)) );  % [-90, +90] deg
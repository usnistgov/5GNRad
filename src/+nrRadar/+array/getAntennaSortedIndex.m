function idxVec = getAntennaSortedIndex(pos)
% GETANTENNASORTEDINDEX Return element indices sorted by spatial position
%   IDX = GETANTENNASORTEDINDEX(ARRAY) returns the linear indices of the 
%   elements in the array ARRAY, sorted according to their physical 
%   position in the array plane. The elements are ordered by increasing 
%   z-coordinate (rows) and then by increasing y-coordinate (columns), as 
%   determined from the element positions returned by GETELEMENTPOSITION.
%
%   The input ARRAY must be a valid phased array or antenna array 
%   System object supporting the GETELEMENTPOSITION method. The positions 
%   are assumed to be expressed in meters in a 3×Nrx matrix, where rows 
%   correspond to [x; y; z] coordinates and columns to elements.
%
%   Example:
%       % Get element ordering for a 4x8 URA
%       h = phased.URA([4 8], [0.5 0.5]);
%       idxVec = getAntennaSortedIndex(h);
%
%   Notes:
%       * The function automatically detects the unique y- and z-coordinates 
%         of the array elements and builds a 2D grid indexed as (z, y).
%       * If duplicate coordinates are found (e.g., due to dual-polarized
%         elements), only the first index is selected.
%       * The output IDX is a column vector of indices corresponding to the 
%         element positions ordered row-wise (increasing z) and then 
%         column-wise (increasing y).
%
%   See also GETELEMENTPOSITION, PHASED.URA, PHASED.ULA.

%   2025 NIST/CTL Steve Blandino
%
%   This file is available under the terms of the NIST License.

% pos = getElementPosition(array);   % 3 x Nrx  (meters)
ypos = pos(2,:).';                        % Nrx x 1
zpos = pos(3,:).';

% Unique coordinate lists (sorted)
ylist = unique(ypos);
zlist = unique(zpos);

Ny = numel(ylist);                        % columns  (along y)
Nz = numel(zlist);                        % rows     (along z)
assert(Ny*Nz == numel(ypos) || Nz*Ny == numel(ypos), ...
       'Element grid size mismatch: Ny×Nz must equal Nrx.');

% Map each (z,y) to its element index (choose first if duplicates, e.g., dual-pol)
idxGrid = nan(Nz, Ny);                    % rows increase with z, cols with y
for iz = 1:Nz
    for iy = 1:Ny
        idx = find( abs(zpos - zlist(iz)) < 1e-9 & abs(ypos - ylist(iy)) < 1e-9 );
        if isempty(idx)
            error('Missing element at z=%.6g, y=%.6g', zlist(iz), ylist(iy));
        end
        idxGrid(iz,iy) = idx(1);          % if multiple (e.g., dual-pol), take the first; or sum later
    end
end
idxVec = idxGrid(:);    
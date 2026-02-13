function idxElemZIncreasing = getRfChainSortedIndex(array, F_RF)
% GETRFCHAINSORTEDINDEX Order RF chains (columns of F_RF) by spatial position.
% Uses getAntennaSortedIndex(array) for element ordering and maps each RF
% chain to a representative position (centroid of its connected elements).
%
% Output:
%   rfOrderZIncreasing: [NRF x 1] permutation of RF chain indices.

idxElemZIncreasing = nrRadar.array.getAntennaSortedIndex(getElementPosition(array));  

if ~isempty(F_RF)
pos = getElementPosition(array).';                  % [Nelem x 3] [x y z]

[Nelem, NRF] = size(F_RF);
assert(numel(idxElemZIncreasing) == Nelem, "Element count mismatch.");

% Compute a representative position per RF chain
rfPos = nan(3,NRF);
tol = 1e-12;

for k = 1:NRF
    w = F_RF(:,k);
    idx = find(abs(w) > tol);
    assert(~isempty(idx), "RF chain %d has no connected elements.", k);

    % Weighted centroid (magnitude)
    ww = abs(w(idx));
    ww = ww / sum(ww);
    rfPos(:,k) = ww.' * pos(idx,:);
end
idxElemZIncreasing = getAntennaSortedIndex(rfPos);  
end
end
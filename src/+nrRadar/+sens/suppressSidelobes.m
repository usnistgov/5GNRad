function pkOut = suppressSidelobes(pk, tolNear, sizes, wrapD, wrapAz, wrapEl, tolSame)
% SUPPRESSSIDELOBES  Merge redundant or side-lobe peaks in 4-D detection space.
%   PKOUT = SUPPRESSSIDELOBES(PK) merges redundant detections in a
%   4-D peak list PK = [r d el az val], where each row represents a
%   detection in range, Doppler, elevation, and azimuth bins with its
%   corresponding value (e.g., power or score).
%
%   PKOUT = SUPPRESSSIDELOBES(PK, TOLNEAR, SIZES, WRAPD, WRAPAZ, WRAPEL, TOLSAME)
%   merges peaks that are close in at least three of the four dimensions.
%
%   Inputs
%   ------
%   pk        : [K x 5] array, each row [r d el az val] with 1-based indices.
%   tolNear   : Scalar, maximum allowed difference (in bins) on the single
%               differing dimension to consider two peaks duplicates.
%               Default = 4.
%   sizes     : [R D El Az] array, used for wrap-aware differences on
%               Doppler, Azimuth, and Elevation (use Inf for non-circular
%               axes). Default = [Inf Inf Inf Inf].
%   wrapD     : Logical, true if Doppler axis is circular (default = true).
%   wrapAz    : Logical, true if Azimuth axis is circular (default = true).
%   wrapEl    : Logical, true if Elevation axis is circular (default = true).
%   tolSame   : Scalar, tolerated difference (bins) on the 3 “equal” dims
%               (default = 2).
%
%   Behavior
%   --------
%   Two peaks are considered duplicates if at least 3 of the 4 coordinates
%   differ by ≤ tolSame bins (after wrap adjustment), and the remaining
%   dimension differs by ≤ tolNear. The weaker of the two peaks is removed.
%
%   Output
%   ------
%   pkOut     : [K' x 5] array of merged peaks, with stronger representatives
%               retained and near-duplicates removed.
%
%   Example
%   --------
%   pkOut = suppressSidelobes(pk, 4, [256 128 32 64], true, true, true, 2);
%
%   See also: nms4dGreedyOnV, clusterPeaks4D
%
%   2025 NIST/CTL Steve Blandino
%
%   This file is available under the terms of the NIST License.

    if nargin < 2 || isempty(tolNear), tolNear = 4; end
    if nargin < 3 || isempty(sizes),   sizes   = [inf inf inf inf]; end
    if nargin < 4 || isempty(wrapD),   wrapD   = true; end
    if nargin < 5 || isempty(wrapAz),  wrapAz  = true; end
    if nargin < 6 || isempty(wrapEl),  wrapEl  = true; end
    if nargin < 7 || isempty(tolSame), tolSame = 2; end

    if isempty(pk), pkOut = pk; return; end

    % Sort by power so we keep the strongest representative
    [~,ord] = sort(pk(:,5), 'descend');
    pk = pk(ord,:);

    K = size(pk,1);
    keep = true(K,1);

    Rsz = sizes(1); Dsz = sizes(2); Elsz = sizes(3); Azsz = sizes(4); %#ok<NASGU>

    for i = 1:K
        if ~keep(i), continue; end
        ri = pk(i,1); di = pk(i,2); eli = pk(i,3); azi = pk(i,4);

        for j = i+1:K
            if ~keep(j), continue; end

            % raw diffs
            dr  = abs(pk(j,1) - ri);
            dd  = abs(pk(j,2) - di);
            del = abs(pk(j,3) - eli);
            da  = abs(pk(j,4) - azi);

            % wrap-aware on D, Az, El (if sizes are finite)
            if wrapD  && isfinite(Dsz),  dd  = min(dd,  Dsz  - dd);  end
            if wrapAz && isfinite(Azsz), da  = min(da,  Azsz - da);  end
            if wrapEl && isfinite(Elsz), del = min(del, Elsz - del); end

            diffs = [dr dd del da];

            % Count dims that are "same enough" (<= tolSame)
            nearSame = (diffs <= tolSame);
            nNear    = sum(nearSame);

            if nNear >= 3
                    keep(j) = false; % suppress weaker duplicate
                % end
            end
        end
    end

    pkOut = pk(keep,:);
end

function peaks = pick_peaks_nms(rdm_pow, noiseEst, alpha, supRad)
% rdm_pow: power RDM (not magnitude)
% noiseEst: CFAR noise estimate (same size)
% alpha: CFAR multiplier (linear)
% supRad: [range_halfwidth, doppler_halfwidth] suppression radius (bins)

% CFAR-thresholded binary map
B = rdm_pow > alpha .* noiseEst;

localMax = rdm_pow >= movmax(movmax(rdm_pow, 3, 1), 3, 2);
cands = localMax & B;

% Sort candidates by strength (desc)
[rows, cols] = find(cands);
vals = rdm_pow(cands);
[vals, ord] = sort(vals, 'descend'); rows = rows(ord); cols = cols(ord);

% NMS: keep a peak, suppress neighbors within supRad
[H, W] = size(rdm_pow);
suppressed = false(size(rdm_pow));


% Pre-allocate peaks array for speed (guess size, trim later)
    peaks = zeros(numel(vals), 3); 
    count = 0;
    
    % Cache suppression radii
    r_rad = supRad(1);
    c_rad = supRad(2);

    for k = 1:numel(vals)
        r = rows(k); 
        c = cols(k);

        if ~suppressed(r,c)
            % Record the peak
            count = count + 1;
            peaks(count, :) = [r, c, vals(k)];

            r_min = max(1, r - r_rad);
            r_max = min(H, r + r_rad);
            c_min = max(1, c - c_rad);
            c_max = min(W, c + c_rad);

            % Mark neighborhood as suppressed
            suppressed(r_min:r_max, c_min:c_max) = true;
        end
    end
    
    % Trim unused rows
    peaks = peaks(1:count, :);
    
end
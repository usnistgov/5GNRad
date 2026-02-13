function out = cluster_peaks_4d(pk, dims, opts)
% CLUSTER_PEAKS_4D  Discard weak peaks and cluster similar ones (DBSCAN).
% pk   : Nx5 matrix, rows = [r d el az val]
% dims : struct with cube dimensions and (optional) wrap flags:
%        .R, .D, .El, .Az    (sizes)
%        .wrapD = true/false (default: true, Doppler circular)
%        .wrapAz = true/false (default: true, Azimuth circular)
% opts : struct of options:
%        .minVal      : discard peaks with val < minVal (default: 0)
%        .norm        : [σr σd σel σaz] for spatial normalization (default [1 1 1 1])
%        .valWeight   : weight of value in distance (0 = ignore val; default 0)
%        .valScale    : scale for value feature (default = robust MAD)
%        .eps         : DBSCAN neighborhood radius (default: 2.5)
%        .minPts      : DBSCAN min points for core (default: 2)
%
% Returns:
%   out.labels        : Nx1 cluster label (0=noise)
%   out.clusters      : table per-cluster with centroid, max, sum, size
%   out.members{k}    : indices of pk rows in cluster k
%   out.pk_kept       : filtered pk after minVal (same order as labels)

    if isempty(pk)
        out = empty_out();
        return;
    end

    % ---------- defaults ----------
    if ~isfield(dims,'wrapD'),  dims.wrapD  = true;  end
    if ~isfield(dims,'wrapAz'), dims.wrapAz = true;  end
    if ~isfield(opts,'minVal'),    opts.minVal    = 0;     end
    if ~isfield(opts,'norm'),      opts.norm      = [1 1 1 1]; end
    if ~isfield(opts,'valWeight'), opts.valWeight = 0;     end
    if ~isfield(opts,'eps'),       opts.eps       = 2.5;   end
    if ~isfield(opts,'minPts'),    opts.minPts    = 2;     end

    % ---------- 1) discard weak peaks ----------
    keep = pk(:,5) >= opts.minVal;
    pk2  = pk(keep,:);
    if isempty(pk2)
        out = empty_out();
        return;
    end

    % ---------- 2) build normalized feature space ----------
    r  = pk2(:,1); d  = pk2(:,2); el = pk2(:,3); az = pk2(:,4); v = pk2(:,5);

    % scales for spatial axes
    sig = opts.norm(:)';               % [σr σd σel σaz]
    sig(sig<=0) = 1;

    % value scale (robust)
    if ~isfield(opts,'valScale') || isempty(opts.valScale)
        medv = median(v);
        madv = median(abs(v - medv)) + eps;
        vScale = 1.4826*madv;          % approx std if Gaussian
    else
        vScale = max(opts.valScale, eps);
    end
    wV = max(opts.valWeight, 0);       % 0.. (0 ignores value in distance)

    % normalized features (value optional)
    % Distance uses custom wrap-aware metric; we store raw coords and pass scales.
    X = [r d el az v];  %#ok<NASGU>

    % ---------- 3) DBSCAN with wrap-aware distance ----------
    labels = dbscan_wrap(pk2(:,1:4), v, dims, sig, wV, vScale, opts.eps, opts.minPts);

    % ---------- 4) summarize clusters ----------
    K = max(labels);
    clusters = zeros(K, 11); % [r_c d_c el_c az_c val_max val_sum size r_med d_med el_med az_med]
    members = cell(K,1);
    for k = 1:K
        idx = find(labels==k);
        members{k} = idx;
        rr = r(idx); dd = d(idx); ee = el(idx); aa = az(idx); vv = v(idx);

        % weighted centroid (by val), rounded to nearest bin for indices
        wsum = sum(vv) + eps;
        w = vv / wsum;

        % Linear (non-circular) centroids
        rc  = sum(w .* rr);
        elc = sum(w .* ee);

        % circular centroids for Doppler/Azimuth 
        if dims.wrapD
            dc = circ_mean_w_bins(dd, w, dims.D);
        else
            dc = sum(w .* dd);
        end
        if dims.wrapAz
            azc = circ_mean_w_bins(aa, w, dims.Az);
        else
            azc = sum(w .* aa);
        end

        % Medians 
        rmed  = median(rr);
        elmed = median(ee);

        % Circular medians for Doppler/Azimuth 
        if dims.wrapD
            dmed = circ_median_w_bins(dd, vv, dims.D);
        else
            dmed = median(dd);
        end
        if dims.wrapAz
            azmed = circ_median_w_bins(aa, vv, dims.Az);
        else
            azmed = median(aa);
        end

        clusters(k,:) = [rc dc elc azc, max(vv), sum(vv), numel(idx), ...
                         rmed dmed elmed azmed];
    end

    clustersTbl = array2table(clusters, 'VariableNames', ...
        {'r_cent','d_cent','el_cent','az_cent','val_max','val_sum','size', ...
         'r_med','d_med','el_med','az_med'});

    out = struct();
    out.labels   = labels;
    out.clusters = clustersTbl;
    out.members  = members;
    out.pk_kept  = pk2;   % same order as labels
end

% DBSCAN 
function labels = dbscan_wrap(Xsp, V, dims, sig, wV, vScale, epsRad, minPts)
% Xsp: Nx4 (r,d,el,az), integer-like positions
% V  : Nx1 values (power)
% Distance: sqrt( (Δr/σr)^2 + (Δd/σd)^2 + (Δel/σel)^2 + (Δaz/σaz)^2 + (wV*Δv/vScale)^2 )

    N = size(Xsp,1);
    labels = zeros(N,1);    % 0=noise
    visited = false(N,1);
    clusterId = 0;

    % precompute scales
    sr = sig(1); sd = sig(2); se = sig(3); sa = sig(4);
    R = dims.R; D = dims.D; E = dims.El; A = dims.Az;

    for i = 1:N
        if visited(i), continue; end
        visited(i) = true;

        % region query
        nbrs = region_query(i);

        if numel(nbrs) < minPts
            labels(i) = 0;           % noise (may later become border)
        else
            clusterId = clusterId + 1;
            labels(i) = clusterId;
            % expand cluster
            S = nbrs(:)';
            j = 1;
            while j <= numel(S)
                q = S(j);
                if ~visited(q)
                    visited(q) = true;
                    nbrs2 = region_query(q);
                    if numel(nbrs2) >= minPts
                        S = [S, setdiff(nbrs2(:)', S)]; %#ok<AGROW>
                    end
                end
                if labels(q) == 0
                    labels(q) = clusterId;
                end
                j = j + 1;
            end
        end
    end

    % --- nested: region query with wrap-aware metric ---
    function nb = region_query(idx)
        dr = abs(Xsp(:,1) - Xsp(idx,1));
        dd = abs(Xsp(:,2) - Xsp(idx,2));
        de = abs(Xsp(:,3) - Xsp(idx,3));
        da = abs(Xsp(:,4) - Xsp(idx,4));

        % circular wrap for Doppler/Azimuth if requested
        if dims.wrapD
            dd = min(dd, D - dd);
        end
        if dims.wrapAz
            da = min(da, A - da);
        end

        dv = abs(V - V(idx));

        % normalized distance
        dist2 = (dr/sr).^2 + (dd/sd).^2 + (de/se).^2 + (da/sa).^2;
        nb = find(dist2 <= (epsRad^2));
    end
end

function out = empty_out()
    out = struct('labels',zeros(0,1), 'clusters', ...
        array2table(zeros(0,11), 'VariableNames', ...
        {'r_cent','d_cent','el_cent','az_cent','val_max','val_sum','size','r_med','d_med','el_med','az_med'}), ...
        'members',{cell(0,1)}, 'pk_kept', zeros(0,5));
end

function m = circ_mean_w_bins(x, w, P)
% Weighted circular mean for bin-like variable x in [0,P) with period P.
% Returns m in [0,P).
    x = x(:); w = w(:); w = w / sum(w);
    ang = 2*pi * (mod(x, P) / P);
    z = sum(w .* exp(1j*ang));
    if abs(z) < 1e-12
        % ill-defined mean → fall back to weighted median
        m = circ_median_w_bins(x, w, P);
        return;
    end
    a = angle(z);                       % [-pi, pi)
    m = mod(P * (a / (2*pi)), P);
end

function m = circ_median_w_bins(x, w, P)
% Weighted circular median in bin-domain with period P.
% Strategy: unwrap around circular mean anchor → weighted linear median → wrap back.
    x = mod(x(:)-1, P)+1;
    w = w(:); w = w / sum(w);

    % anchor from unweighted circular mean for stability
    ang = 2*pi * (x / P);
    z = mean(exp(1j*ang));
    if abs(z) < 1e-12
        [~,ix] = max(w); anchor = x(ix);
    else
        anchor = mod(P * (angle(z)/(2*pi)), P);
    end

    % unwrap around anchor into (-P/2, P/2]
    y = x - anchor;
    y = y - P .* round(y / P);

    % weighted linear median
    [ys, ix] = sort(y);
    ws = w(ix) / sum(w);
    cs = cumsum(ws);
    j = find(cs >= 0.5, 1, 'first');
    ym = ys(j);

    m = mod(ym + anchor-1, P)+1;
end
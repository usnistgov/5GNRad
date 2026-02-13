function exportBFtoJSON(filename, cfg)
% cfg.meta: struct
% cfg.beamformer: struct with fields:
%   .type = "full_digital" or "hybrid"
%   .wElem = complex vector (Nelem x 1)  [REQUIRED]
%   .F_RF  = complex matrix (Nelem x Nrf) [OPTIONAL]
%   .wBB   = complex vector (Nrf x 1)     [OPTIONAL]
% cfg.combiner: same fields as beamformer (optional)

    out = struct();

    % ---- meta ----
    out.meta = cfg.meta;

    % ---- beamformer ----
    out.beamformer = packSide(cfg.beamformer);

    % ---- combiner (optional) ----
    if isfield(cfg,'combiner') && ~isempty(cfg.combiner)
        out.combiner = packSide(cfg.combiner);
    end

    txt = jsonencode(out, "PrettyPrint", true);
    fid = fopen(filename, 'w');
    assert(fid>0, "Cannot open %s for writing", filename);
    fwrite(fid, txt, 'char');
    fclose(fid);
end

function sideOut = packSide(sideIn)
    sideOut = struct();
    sideOut.type  = string(sideIn.type);

    % canonical effective weights (required)
    assert(isfield(sideIn,'wElem') && ~isempty(sideIn.wElem), ...
        "side must include wElem");
    sideOut.wElem = packComplex(sideIn.wElem);

    % optional factorization
    if isfield(sideIn,'F_RF') && ~isempty(sideIn.F_RF) && ...
       isfield(sideIn,'wBB')  && ~isempty(sideIn.wBB)
        sideOut.factorization = struct();
        sideOut.factorization.F_RF = packComplex(sideIn.F_RF);
        sideOut.factorization.wBB  = packComplex(sideIn.wBB);
    end
end

function s = packComplex(X)
% Store complex array generically with shape + flattened re/im (row-major in JSON).
    sz = size(X);
    s = struct();
    s.shape = sz(:).';
    s.re = real(X(:)).';   % flatten column-major; shape tells how to reshape
    s.im = imag(X(:)).';
end

function cfg = importBFfromJSON(filename)
    txt = fileread(filename);
    in  = jsondecode(txt);

    cfg = struct();
    cfg.meta = in.meta;

    cfg.beamformer = unpackSide(in.beamformer);

    if isfield(in,'combiner')
        cfg.combiner = unpackSide(in.combiner);
    end
end

function side = unpackSide(sideIn)
    side = struct();
    side.type = string(sideIn.type);

    side.wElem = unpackComplex(sideIn.wElem);
    if isfield(sideIn,'manifold')
        side.manifold = unpackComplex(sideIn.manifold);
    end
    if isfield(sideIn,'factorization')
        side.F_RF = unpackComplex(sideIn.factorization.F_RF);
        side.wBB  = unpackComplex(sideIn.factorization.wBB);
    else
        side.F_RF = [];
        side.wBB  = [];
    end
end

function X = unpackComplex(s)
    re = s.re(:);
    im = s.im(:);
    X = complex(re, im);
    X = reshape(X, s.shape(:)');
end

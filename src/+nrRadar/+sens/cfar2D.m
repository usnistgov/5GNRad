function [y, noiseEst] = cfar2D(rdm, guardRadius, trainingRadius, threshold)
% CFAR2D_FAST  Vectorized 2-D CA-CFAR (Cell-Averaging) with variable window size.
%
%   [Y, NOISEEST] = CFAR2D_FAST(RDM, GUARDRADIUS, TRAININGRADIUS, THRESHOLD)
%   performs 2-D CA-CFAR detection on the power map RDM (linear scale).
%
%   Unlike the strict "valid-only" version, this implementation adapts the
%   number of training cells near edges, so detections are still possible
%   even when the full training window does not fit inside RDM.
%
%   INPUTS
%     RDM            : M×N real, nonnegative matrix (linear power).
%     GUARDRADIUS    : [gR gC] half-sizes of guard window (in cells).
%     TRAININGRADIUS : [tR tC] half-sizes of training window (>= guard).
%     THRESHOLD      : Scalar multiplier in linear domain. For +12 dB use
%                      THRESHOLD = 10^(12/10).
%
%   OUTPUTS
%     Y         : Same size as RDM, equals RDM at detections and 0 elsewhere.
%     NOISEEST  : Estimated local clutter power at each pixel (same size).
%
%   EXAMPLE
%     thr_dB = 12;
%     thr = 10^(thr_dB/10);
%     [y,nEst] = cfar2D_fast(RDM, [2 2], [6 6], thr);
%     imagesc(y>0); axis image; title('Detections incl. edges');
%
%   2025 Steve Blandino / NIST CTL
%   This file is available under the NIST License.

% ---- Radii and checks ----
gR = guardRadius(1); gC = guardRadius(2);
tR = trainingRadius(1); tC = trainingRadius(2);
assert(tR>=gR && tC>=gC, 'Training radii must be >= guard radii.');

[M,N] = size(rdm);

% ---- Box kernels ----
kerT = ones(2*tR+1, 2*tC+1, 'like', rdm);  % training square
kerG = ones(2*gR+1, 2*gC+1, 'like', rdm);  % guard square

% % ---- Pad data and unity mask (replicate to avoid zero bias) ----
idxR = [ones(1, tR), 1:M, M*ones(1, tR)];
idxC = [ones(1, tC), 1:N, N*ones(1, tC)];

% Apply padding via indexing
rdmPad = rdm(idxR, idxC);
onePad = ones(length(idxR), length(idxC), 'like', rdm);

% ---- Training sums/counts (VALID -> MxN) ----
sumT = conv2(rdmPad, kerT, 'valid');   % size MxN
cntT = conv2(onePad, kerT, 'valid');   % size MxN

% ---- Guard sums/counts (VALID -> larger, then center-crop to MxN) ----
sumG_full = conv2(rdmPad, kerG, 'valid');       % size M+2*(tR-gR) by N+2*(tC-gC)
cntG_full = conv2(onePad, kerG, 'valid');

dr = (tR - gR);
dc = (tC - gC);
rows = (1:M) + dr;   % [dr+1 : dr+M]
cols = (1:N) + dc;   % [dc+1 : dc+N]

sumG = sumG_full(rows, cols);    % now MxN
cntG = cntG_full(rows, cols);    % now MxN

% ---- Training-ring sums and adaptive K per CUT ----
ringSum  = sumT - sumG;
Kmap     = max(cntT - cntG, 1);          % avoid divide-by-zero
noiseEst = ringSum ./ Kmap;

% ---- Threshold test everywhere (linear domain) ----
detect = rdm > noiseEst * threshold;

% ---- Output map: keep original power at detections ----
y = zeros(size(rdm), 'like', rdm);
y(detect) = rdm(detect);
end

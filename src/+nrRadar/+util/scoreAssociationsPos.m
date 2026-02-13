function [metrics, assoc] = scoreAssociationsPos(gtPos, detPos,gtVelVec, detVel, sensorPos, varargin)
% SCOREASSOCIATIONSPOS  Associate detections to ground truth in 3-D and score results.
%   [METRICS, ASSOC] = SCOREASSOCIATIONSPOS(GTPOS, DETPOS, GTVEL, DETVEL, SENSORPOS)
%   computes a gated Mahalanobis-distance cost between each ground-truth
%   position (rows of GTPOS) and each detection position (rows of DETPOS),
%   performs Hungarian assignment (Munkres) on the gated cost matrix, and
%   returns summary metrics and bookkeeping of the associations. Detected
%   positions/velocities are aligned to the order of GTPOS and passed to
%   EVALTP_POSVEL to produce per-target error statistics.
%
%   [...] = SCOREASSOCIATIONSPOS(..., Name, Value) specifies options using
%   one or more name-value arguments:
%     'SigmaXYZ'  — Position standard deviation(s) used to build a diagonal
%                   covariance when 'Cov' is not given. Scalar or 3-vector.
%                   Default: 6  (same units as GTPOS/DETPOS, e.g., meters)
%     'Cov'       — 3×3 position-error covariance matrix (overrides
%                   'SigmaXYZ' when provided). Default: []
%     'PGate'     — Gating probability for a χ² gate with 3 DoF. The gate
%                   value is gate2 = chi2inv(PGate,3). Default: 0.9973
%                   (≈ 3-σ in 3-D).
%     'BigCost'   — Large finite cost used to represent "out of gate"
%                   pairs in the assignment matrix. Default: 1e9
%
%   Inputs
%   ------
%   GTPOS     : [Ngt×3] double. Ground-truth Cartesian positions.
%   DETPOS    : [Ndet×3] double. Detection Cartesian positions.
%   GTVEL     : [Ngt×1] or [Ngt×?] velocity info for EVALTP_POSVEL
%               (pass [] if not used by your implementation).
%   DETVEL    : [Ndet×1] velocity info for EVALTP_POSVEL
%               (pass [] if not available).
%   SENSORPOS : [1×3] or [3×1] sensor/reference position used by
%               EVALTP_POSVEL. 
%
%   Outputs
%   -------
%   METRICS : struct with summary counts and rates
%       .TP, .FN, .FP                 — true/false positives/negatives
%       .TPR, .FNR, .FPR             — rates computed w.r.t. Ngt/Ndet
%       .stats                        — struct returned by EVALTP_POSVEL
%                                       
%
%   ASSOC   : struct with association bookkeeping
%       .gtToDet   — [Ngt×1] detection index matched to each GT (0 if none)
%       .detToGt   — [Ndet×1] GT index matched to each detection (0 if none)
%       .costMat   — [Ngt×Ndet] gated cost matrix (Mahalanobis² or BIG)
%       .gate2     — χ² gate threshold (3 DoF)
%       .SigmaXYZ  — Sigma used (echo of input)
%       .Cov       — 3×3 covariance used after defaults/regularization
%       .InvCov    — inverse of Cov used for Mahalanobis²
%
%   Method
%   ------
%   * Build Mahalanobis-squared distance D² = (x−y)ᵀ Σ⁻¹ (x−y) for all
%     GT–detection pairs using either the provided covariance 'Cov' or a
%     diagonal covariance from 'SigmaXYZ'.
%   * Apply χ²(3) gating at probability 'PGate'; out-of-gate pairs receive
%     cost = 'BigCost'.
%   * Solve the assignment with Munkres/Hungarian on the gated costs.
%   * Align DETPOS/DETVEL to GT order (NaN where unmatched) and call
%     EVALTP_POSVEL to compute detailed error statistics.
%
%   Notes
%   -----
%   * Requires the function MUNKRES (Hungarian assignment)
%   * Positions must be [N×3]. Units are user-defined but must be
%     consistent across inputs and 'SigmaXYZ'/'Cov'.
%
%   Example
%   -------
%       gtPos  = [0 0 0; 10 0 0];
%       detPos = [0.5 0 0; 9.6 0.2 0];
%       [metrics, assoc] = scoreAssociationsPos(gtPos, detPos, [], [], []);
%       metrics.TP   % -> 2
%       assoc.gtToDet
%
%   See also MUNKRES, CHI2INV, EVALTP_POSVEL.
%
%   2025 NIST/CTL Steve Blandino
%   This file is available under the terms of the NIST License.

% -------- parameters
p = inputParser;
p.addParameter('SigmaXYZ', 6, @(x)isnumeric(x)&&(isscalar(x)||(isvector(x)&&numel(x)==3)));
p.addParameter('Cov', [], @(x)isnumeric(x)&&isequal(size(x),[3 3]));
p.addParameter('PGate', 0.9973, @(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
p.addParameter('BigCost', 1e9, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.parse(varargin{:});
sigXYZ = p.Results.SigmaXYZ;
Cov    = p.Results.Cov;
pGate  = p.Results.PGate;
BIG    = p.Results.BigCost;

% -------- validate inputs
gtPos  = double(gtPos);
detPos = double(detPos);
assert(size(gtPos,2)==3 && size(detPos,2)==3, 'gtPos/detPos must be [N x 3].');
Ngt  = size(gtPos,1);
Ndet = size(detPos,1);


gtVel = sum((-sensorPos + gtPos) ./ vecnorm(sensorPos - gtPos, 2, 2) .* gtVelVec, 2);


if Ngt==0
    metrics = struct('TP',0,'FN',0,'FP',Ndet,'TPR',0,'FNR',0,'FPR',1, ...
        'absErrorXYZ',zeros(0,3),'absError3D',zeros(0,1), ...
        'rmseXYZ',[NaN NaN NaN],'rmse3D',NaN,'maeXYZ',[NaN NaN NaN], ...
        'mae3D',NaN,'biasXYZ',[NaN NaN NaN]);
    assoc = struct('gtToDet', zeros(0,1), 'detToGt', zeros(Ndet,1), ...
                   'costMat', zeros(0,Ndet), 'gate2', NaN, ...
                   'SigmaXYZ', sigXYZ, 'Cov', Cov, 'InvCov', []);
    return;
end

% -------- covariance / inverse covariance
if ~isempty(Cov)
    CovUse = Cov;
else
    if isscalar(sigXYZ), sigXYZ = [sigXYZ sigXYZ sigXYZ]; end
    CovUse = diag(sigXYZ(:).^2);
end
% regularize
CovUse = (CovUse + CovUse.')/2;
% Eigen decomposition
[U,S] = eig(CovUse,'vector');
S = max(S, 1e-12*max(S));      % clip tiny eigenvalues
InvCov = U * diag(1 ./ S) * U';

% -------- build Mahalanobis-squared distances (3D)
costMat = zeros(Ngt, Ndet);
for j = 1:Ndet
    dX = gtPos - detPos(j,:);                % [Ngt x 3]
    % Mahalanobis squared via (dX*InvCov) dot dX row-wise
    t  = dX * InvCov;                        % [Ngt x 3]
    costMat(:,j) = sum(t .* dX, 2);          % [Ngt x 1]
end

% -------- gating: chi-square with dof=3
gate2 = chi2inv(pGate, 3);
costMat_gated = costMat;
costMat_gated(costMat > gate2) = BIG;

% -------- Hungarian assignment
assign  = nrRadar.util.munkres(costMat_gated);                 % [1 x Ndet] (col indices)
gtToDet = assign(:);
gtToDet(~isfinite(gtToDet)) = 0;
gtToDet = round(gtToDet);
gtToDet(gtToDet < 1 | gtToDet > Ndet) = 0;

detToGt = zeros(Ndet,1);
for i = 1:Ngt
    j = gtToDet(i);
    if j>=1 && j<=Ndet && detToGt(j)==0
        detToGt(j) = i;
    end
end

% -------- gating validity per assigned pair
valid = false(Ngt,1);
rows  = find(gtToDet > 0);                         % only rows with an assignment
if ~isempty(rows)
    cols   = gtToDet(rows);
    linIdx = sub2ind([Ngt, Ndet], rows, cols);
    valid(rows) = isfinite(costMat_gated(linIdx)) & (costMat_gated(linIdx) < BIG);
end

TP = sum(valid);
FN = Ngt - TP;
FP = Ndet - TP;

% -------- build detPos/detVel ALIGNED to gtPos (NaN where no valid match)
detPosAligned = nan(Ngt, 3);
detVelAligned = nan(Ngt, 1);

assigned = (gtToDet > 0) & valid;
if any(assigned)
    jj = gtToDet(assigned);                         % det indices matched to these GT rows
    detPosAligned(assigned, :) = detPos(jj, :);
    if exist('detVel','var') && ~isempty(detVel)
        detVelAligned(assigned) = detVel(jj);
    end
end

errStruct = nrRadar.util.evalTP_PosVel(gtPos, detPosAligned, gtVel, detVelAligned, sensorPos);

% -------- summary metrics
metrics = struct();
metrics.TP  = TP;
metrics.FN  = FN;
metrics.FP  = FP;
metrics.TPR = TP / max(Ngt,1);
metrics.FNR = FN / max(Ngt,1);
metrics.FPR = FP / max(Ndet,1);
metrics.stats = errStruct;
metrics.info.gtVel = gtVel;

assoc = struct('gtToDet', gtToDet, 'detToGt', detToGt, ...
               'costMat', costMat_gated, 'gate2', gate2, ...
               'SigmaXYZ', sigXYZ, 'Cov', CovUse, 'InvCov', InvCov);
end

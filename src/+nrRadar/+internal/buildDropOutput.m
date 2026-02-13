function drop = buildDropOutput(q, metrics, snrvarVec, ...
                                          stPosition, stPositionEstimate_k, ...
                                          velParallel_k, peakRda, ...
                                          nStDrop)
%RUN5GNRAD_BUILDDROPOUTPUT Package simulation results for a single drop.
%
%   DROP = RUN5GNRAD_BUILDDROPOUTPUT(Q, METRICS, SNRVARVEC, ...)
%   collects error statistics, detection metrics, and payload data into
%   a single struct for the current drop Q.
%
%   Inputs:
%       q                  - (Scalar) Current drop index
%       metrics            - Struct containing calculated stats (pos, vel, TP/FN/FP)
%       snrvarVec          - [N_targets x 1] Vector of SNR values (dB)
%       stPosition         - [N_targets x 3] Ground truth positions for this drop
%       stPositionEstimate - [N_det x 3] Estimated positions
%       gtVr               - [N_targets x 1] Ground truth radial velocities
%       velParallel_k      - [N_det x 1] Estimated radial velocities
%       peakRda            - [N_det x 1] Peak RDA values for detections
%       simConfig          - Simulation config (for nStDrop)
%
%   Output:
%       drop               - Struct containing all drop-specific results.

    %% Initialize & Basic Indices
    targetIdx = (q-1)*nStDrop + (1:nStDrop);
    
    % Initialize output struct
    drop = struct();

    % Metadata
    drop.q         = q;
    drop.targetIdx = targetIdx(:);
    drop.timeIndex = repmat(q, nStDrop, 1);

    %% Position & Velocity Errors
    errXYZ = metrics.stats.pos.errXYZ;

    % Cartesian Errors
    drop.positionErrorX = errXYZ(:,1);
    drop.positionErrorY = errXYZ(:,2);
    drop.positionErrorZ = errXYZ(:,3);

    % Horizontal/Vertical Split
    drop.positionErrorH = vecnorm(errXYZ(:,1:2), 2, 2);
    drop.positionErrorV = errXYZ(:,3);

    % Polar/Spherical Errors
    drop.rangeError     = metrics.stats.pos.range_err;
    drop.velocityError  = metrics.stats.vel.vr_err;
    drop.azimuthError   = metrics.stats.pos.az_err_deg;
    drop.elevationError = metrics.stats.pos.el_err_deg;

    %% Detection Performance Metrics
    drop.TP  = metrics.TP;
    drop.FN  = metrics.FN;
    drop.FP  = metrics.FP;
    drop.FPR = metrics.FPR;
    
    % Detection Count (Debug)
    drop.nDetectedTarget = size(stPositionEstimate_k, 1);

    %% Signal Quality
    % Store the SNR slice corresponding to these targets
    drop.snrSlice = snrvarVec(:); 

    %% Detailed Payload (for Debug/Visualization)
    % This mirrors the structure used for detailed post-processing
    drop.detectionPayload = struct( ...
        'q',                  q, ...
        'targetIdx',          targetIdx(:), ...
        'gtPosition',         stPosition, ...
        'stPositionEstimate', stPositionEstimate_k, ...
        'gtVel',              metrics.info.gtVel, ...
        'velParallel',        velParallel_k(:), ...
        'peakRda',            peakRda(:) ...
    );

end
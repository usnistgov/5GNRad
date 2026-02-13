function cache = dropWorker_rdmOnly(q, txWaveform, staticParams, simConfig,...
    sensConfig, stPosition, velocityVectors, backgroundChannel, targetChannel)
% DROPWORKER_RDMONLY Compute range-Doppler cube for one drop (no detection).
%
%   CACHE = DROPWORKER_RDMONLY(...) runs the channel + receiver processing
%   chain up to the range-Doppler cube computation, then returns the cube
%   and the metadata required to run rdmDetection later.
%
%   This is intended for detection-parameter sweeps where the expensive
%   channel/receiver processing is identical and only detection thresholds
%   (CFAR, NMS, clustering, DOA method, etc.) are varied.
%
%   Returned fields
%     cache.q              - drop index
%     cache.stPosition     - ground-truth positions for this drop
%     cache.velocityVectors- ground-truth velocities for this drop
%     cache.stSNR          - per-snapshot SNR from the channel builder (dB)
%     cache.syncOffset     - range sync offset used by rdmDetection
%     cache.receiveArray   - phased array object for DOA estimation
%     cache.analogAngle    - analog beam angle estimate (if applicable)
%     cache.rdmCube        - range-Doppler cube (complex)
%
%   2026 NIST/CTL Steve Blandino
%   This file is available under the terms of the NIST License.

%% Channel
[H, stSNR, syncOffset, channelLen, receiveArray, analogAngle] = ...
    nrRadar.channel.buildCdlChannel(targetChannel, backgroundChannel, ...
    simConfig.txAntenna, simConfig.rxAntenna, simConfig.trpYawDeg, ...
    velocityVectors, staticParams, 'bsPos', staticParams.txPos, ...
    'tgtPos', stPosition);

%% Receive signal
rxWaveform  = nrRadar.rx.getRxWaveform(txWaveform, H, staticParams);

%% Channel estimation
channelEstimate = nrRadar.rx.estimateChannel(rxWaveform, staticParams, 'channelLength', channelLen);
clear rxWaveform

%% Range Doppler cube
rdmCube = nrRadar.sens.getRangeDoppler(channelEstimate, staticParams);
clear channelEstimate

%% Package
cache = struct();
cache.q               = q;
cache.stPosition      = stPosition;
cache.velocityVectors = velocityVectors;
cache.stSNR           = stSNR;
cache.syncOffset      = syncOffset;
cache.receiveArray    = receiveArray;
cache.analogAngle     = analogAngle;
cache.rdmCube         = rdmCube;

end

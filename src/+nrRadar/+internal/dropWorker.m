function drop = dropWorker(q, txWaveform, staticParams, simConfig,...
    sensConfig, stPosition, velocityVectors, backgroundChannel, targetChannel)
%RUN5GNRAD_DROPWORKER Run one "drop" (q-th snapshot group) of the 5GNRad pipeline.
%
%   DROP = RUN5GNRAD_DROPWORKER(Q, PRE, SIMCONFIG, SENSCONFIG, STPOSITION,
%   VELOCITYVECTORS, BACKGROUNDCHANNEL, TARGETCHANNEL)
%
%   Designed to be parfor-safe: no shared writes, no fprintf. Everything
%   returned in DROP is self-contained for post-processing/concatenation.
%
%   Required fields in PRE (from run5GNRad_precompute):
%     c, fc, lambda, txPos, rxPos
%     carrier, prs, ofdmInfo, ofdmFftLen, cpLengths, ofdmTs, sampleRate
%     numSlots, numSymPerSlot, nSymTot, symbolIndices, nProcSym
%     ofdmGrid, txWaveform, txSymbolStride
%     startIdx, endIdx, numberSubcarriers
%     rangeFFTLen, rangeWindow, rangeBinDestgrd, dopplerFftLen, dopplerWindow3D, velocityBin
%     Nfft_v, Nfft_h, W, azGrid, elGrid
%     prsIndPerSlot, symIndSlotPerSlot, slotIdxList, isEvenCombSym
%     rows, snrvar (linear)
%
%   Notes:
%     - backgroundChannel(q) is passed through to getSensingCdl as in your code.
%     - targetChannel is optional; if empty, it is not passed.

%% Get target+baground channel
[H, stSNR,syncOffset, channelLen,receiveArray,analogAngle] = ...
    nrRadar.channel.buildCdlChannel(targetChannel, backgroundChannel, ...
    simConfig.txAntenna, simConfig.rxAntenna,simConfig.trpYawDeg,...
    velocityVectors, staticParams, 'bsPos', staticParams.txPos, ...
    'tgtPos', stPosition);

%% Receive Signal
rxWaveform  = nrRadar.rx.getRxWaveform(txWaveform, H, staticParams);

%% Channel estimation
channelEstimate = nrRadar.rx.estimateChannel(rxWaveform, staticParams, 'channelLength', channelLen);
clear rxWaveform

%% Range Doppler Map
rangeDopplerMap = nrRadar.sens.getRangeDoppler(channelEstimate, staticParams);
clear channelEstimate

%% Detection and parameter estimate
[stPositionEstimate, velParallel, stPostProcSNR] = nrRadar.sens.rdmDetection...
    (rangeDopplerMap, sensConfig, simConfig, receiveArray, syncOffset, staticParams, 'angleEstimate',analogAngle);

%% Metrics 
metrics = nrRadar.util.scoreAssociationsPos(stPosition, stPositionEstimate, velocityVectors, velParallel, staticParams.txPos);

%% Package 
drop = nrRadar.internal.buildDropOutput(q, metrics, stSNR, ...
                                          stPosition, stPositionEstimate, ...
                                          velParallel, stPostProcSNR, ...
                                          simConfig.nStDrop);

end

function [H, snr, syncOffset, channelLen,receiveArray, angleEstimate] = ...
    buildCdlChannel(targetChannel, backgroundChannel, txAntenna, rxAntenna,trpYawDeg,...
                                            velocityVectors, staticParams, varargin)
% BUILDLCHANNEL Build full sensing CDL channel response and derived metrics.
%   [H,SNR,SYNCOFFSET,CHANNELLEN,RECEIVEARRAY,ANGLEESTIMATE] = ...
%   BUILDLCHANNEL(TARGETCHANNEL,BACKGROUNDCHANNEL,TXANTENNA,RXANTENNA,TRPYAWDEG,...
%   VELOCITYVECTORS,STATICPARAMS) constructs the full channel impulse response
%   matrix H for a sensing scenario using a clustered delay line (CDL)
%   channel model. The function wraps NR radar channel generation and
%   returns additional quantities such as SNR, timing synchronization offset,
%   effective channel length, the instantiated receive array, and optional
%   angle estimates.
%
%   Inputs:
%     TARGETCHANNEL       - Target channel structure. If empty ([]), the
%                           target channel is omitted and only the provided
%                           background channel is used.
%     BACKGROUNDCHANNEL   - Background (environment) channel structure.
%     TXANTENNA           - Transmit antenna/array configuration passed to
%                           the channel generator.
%     RXANTENNA           - Receive antenna/array configuration passed to
%                           the channel generator.
%     TRPYAWDEG           - TRP yaw rotation angle in degrees used to rotate
%                           angles into the local TRP frame.
%     VELOCITYVECTORS     - Velocity vectors used by the sensing channel
%                           generator.
%     STATICPARAMS        - Structure containing static parameters:
%                           * fc            : carrier frequency (Hz)
%                           * ofdmTs        : OFDM sample time (s)
%                           * ofdmInfo      : struct with field SampleRate
%                           * symbolIndices : OFDM symbol indices
%                           * snrvar        : noise variance (linear)
%
%   Name-Value Pairs:
%     'angleEstimation'   - Angle estimation method. One of the supported
%                           methods in validAngleEstimation. Default 'ideal'.
%     'bsPos'             - Base-station / TRP position (implementation-defined).
%                           Default [].
%     'tgtPos'            - Target position (implementation-defined).
%                           Default [].
%
%   Outputs:
%     H               - Full channel impulse response (time-domain) as
%                           returned by nrRadar.channel.getSensingCdl.
%     SNR                 - Signal-to-noise ratio in dB computed as
%                           10*log10(STATICPARAMS.snrvar*TgtPwr), where TgtPwr
%                           is returned by the channel generator.
%     SYNCOFFSET          - Timing synchronization offset returned by the
%                           channel generator.
%     CHANNELLEN          - Effective channel length, equal to size(HFULL,1).
%     RECEIVEARRAY        - Receive array object/configuration returned by
%                           the channel generator (may differ from RXANTENNA
%                           depending on generator behavior).
%     ANGLEESTIMATE       - Estimated angles (format depends on the selected
%                           angle estimation method).
%
%   See also NRRADAR.CHANNEL.GETSENSINGCDL
%
%   2026 NIST/CTL Steve Blandino
%
%   This file is available under the terms of the NIST License.

p = inputParser;

addParameter(p, 'angleEstimation', 'ideal', @(x) any(validatestring(x, validAngleEstimation)));
addParameter(p, 'bsPos', []);
addParameter(p, 'tgtPos', []);
% s

parse(p, varargin{:});

angleEstimation      = p.Results.angleEstimation;
bsPos      = p.Results.bsPos;
tgtPos      = p.Results.tgtPos;

    %% Unpack Static Parameters
    % Extract constants needed for the base arguments
    fc            = staticParams.fc;
    ofdmTs        = staticParams.ofdmTs;
    sampleRate    = staticParams.ofdmInfo.SampleRate;
    symbolIndices = staticParams.symbolIndices;

    %% Construct Arguments for getSensingCdl
    
    % Base positional/velocity arguments
    baseArgs = { ...
        velocityVectors, ...
        fc, ...
        ofdmTs ...
    };

    % Configuration Options
    opts = struct( ...
        'bandwidth',            sampleRate, ...
        'transmitArray',        txAntenna, ...
        'receiveArray',         rxAntenna, ...
        'backgroundChannel',    backgroundChannel, ...
        'trpYawDeg',            trpYawDeg, ...
        'symbolIndices',        symbolIndices, ...
        'angleEstimation',      angleEstimation,...
        'bsPos',                bsPos, ...
        'tgtPos',               tgtPos ...
    );

    % Append target channel only if it exists/is non-empty
    if ~isempty(targetChannel)
        opts.targetChannel = targetChannel;
    end

    %% Call Generator
    
    fullArgs = [baseArgs, namedargs2cell(opts)];
    
    [H, tgtPwr, syncOffset, receiveArray, angleEstimate] = ...
        nrRadar.channel.getSensingCdl(fullArgs{:});
    snr = 10*log10(staticParams.snrvar*tgtPwr);
    channelLen = size(H,1);

end
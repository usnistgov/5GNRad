function [outCirAntennaPort, tgtPower, syncOffset,rxArrayObj, out] = getSensingCdl(tgtVel, fc, pri, varargin)
% GETSENSINGCDL Generate sensing CIR
%
% Supports:
%   - Single TRP (default)
%   - Multiple co-located TRPs with different array yaw (boresight) angles
%     computed efficiently: MPCs generated once, then per-TRP beamforming applied.
%
% Multi-TRP I/O:
%   outCir:            [nSamples x nRealization x NtrpSelected]
%   outCirAntennaPort: [nSamples x nRealization x NRF x NtrpSelected]
%   out:               1xNtrpSelected struct array
%
% Name-value additions:
%   'trpYawDeg'      : yaw angles for each TRP (deg). Default 0.
%                      Can be scalar or vector (length Ntrp).
%   'trpSelect'      : indices of TRPs to compute. Default 1:numel(trpYawDeg).
%
% Existing NV (kept):
%   'bandwidth','scenario','backgroundChannel','targetChannel','angleEstimation','scanvector',
%   'transmitArray','receiveArray','wElemTx','wElemRx','rxAnalogCombiner', ...

% ---------------------------- Input parser --------------------------------
p = inputParser;
validAngleEstimation = {'ideal','nearest','scan'};

addParameter(p, 'bandwidth', []);
addParameter(p, 'bistaticAngle', 0);
addParameter(p, 'aspectAngle', 0);
addParameter(p, 'nRcsSamples', 1);
addParameter(p, 'nRealization', 1);
addParameter(p, 'transmitArray', phased.ULA(64));
addParameter(p, 'receiveArray', phased.ULA(64));
addParameter(p, 'firstScatterRangeMax', 30);
addParameter(p, 'scanvector', 0);
addParameter(p, 'scenario', 'UMaAV');
addParameter(p, 'angleEstimation', 'ideal', @(x) any(validatestring(x, validAngleEstimation)));
addParameter(p, 'backgroundChannel', []);
addParameter(p, 'targetChannel', []);
addParameter(p, 'returnElementCIR', true);
addParameter(p, 'symbolIndices', []);
addParameter(p, 'customAntennaPattern', true);
addParameter(p, 'tgtRayPowerThr', 40);
addParameter(p, 'bsPos', []);
addParameter(p, 'tgtPos', []);

% Multi-TRP controls (co-located)
addParameter(p, 'trpYawDeg', 0);          % scalar or vector [Ntrp x 1]
addParameter(p, 'trpSelect', []);         % indices into trpYawDeg (default all)

parse(p, varargin{:});

nRealization      = p.Results.nRealization;
bandwidth         = p.Results.bandwidth;
transmitArray     = p.Results.transmitArray;
receiveArrayConfig      = p.Results.receiveArray;
scenario          = p.Results.scenario;
angleEstimation   = p.Results.angleEstimation;
backgroundChannel = p.Results.backgroundChannel;
targetChannel     = p.Results.targetChannel;
symbolIndices       = p.Results.symbolIndices;
tgtRayPowerThr = p.Results.tgtRayPowerThr;
bsPos = p.Results.bsPos;
tgtPos = p.Results.tgtPos;

if isempty(symbolIndices)
    symbolIndices = 1:nRealization;
else
    nRealization = symbolIndices(end);
end

trpYawDegAll      = p.Results.trpYawDeg(:).';
if isempty(trpYawDegAll)
    trpYawDegAll = 0;
end
if isempty(p.Results.trpSelect)
    trpSelect = 1:numel(trpYawDegAll);
else
    trpSelect = p.Results.trpSelect(:).';
end
NtrpSel = numel(trpSelect);

% ---------------------------- Parameters ----------------------------------
c = 299702547;
lambda = c/fc;

if ~isempty(bandwidth)
    sampleRate = bandwidth;
else
    sampleRate = 100e6;
end
timeVector = 0:pri:pri*(nRealization-1);

% ---------------------------- Target model --------------------------------
if isempty(targetChannel)
    tgtRcs = nrRadar.channel.getSigmaRCS('uav-small', 'returnLargeScale', 1);

    tgtRange   = norm(tgtPos - bsPos);
    tgtDelay   = 2*tgtRange/c;
    bs2tgtVec  = (tgtPos - bsPos)/tgtRange;
    vRel       = dot(tgtVel, bs2tgtVec,2);
    fD         = 2*(vRel/lambda);

    [PG, HasLOSCluster] = nrRadar.channel.computePathLoss(scenario, fc, bsPos, tgtPos, 0);
    tgtPg = -(2*PG - tgtRcs + 10*log10(c^2/(4*pi*fc^2)));

    dopplerPhase = exp(1j * 2 * pi * fD * timeVector);

    vx = bs2tgtVec(:,1); vy = bs2tgtVec(:,2); vz = bs2tgtVec(:,3);
    aoaAzTgt = atan2d(vy, vx);
    aoaElTgt = atan2d(vz, sqrt(vx.^2 + vy.^2));

    aodAzTgt = aoaAzTgt;
    aodElTgt = aoaElTgt;

    aoaAzLOS = aoaAzTgt(1);
    aoaElLOS = aoaElTgt(1);
    tgtPower = 10^(tgtPg/10);

else

    tgtDelay  = vertcat(targetChannel.PathDelays);

    aoaAzLOS  = arrayfun(@(s) s.AnglesAoA(1), targetChannel).';
    aoaAzTgt  = wrapTo180(vertcat(targetChannel.AnglesAoA));

    aoaElLOS  = 90 - arrayfun(@(s) s.AnglesZoA(1), targetChannel).';
    aoaElTgt  = 90 - vertcat(targetChannel.AnglesZoA);

    aodAzTgt  = wrapTo180(vertcat(targetChannel.AnglesAoD));
    aodElTgt  = 90 - vertcat(targetChannel.AnglesZoD);

    dodVector = nrRadar.util.angle2vector(aodAzTgt, 90-aodElTgt, 1);
    dodVector = dodVector ./ vecnorm(dodVector,2,2);

    doaVector = nrRadar.util.angle2vector(aoaAzTgt, 90-aoaElTgt, 1);
    doaVector = doaVector ./ vecnorm(doaVector,2,2);

    nMpc = arrayfun(@(s) size(s.AnglesZoA,1), targetChannel);
    nt = length(nMpc);
    idx  = repelem((1:nt).', nMpc);                 % [Ntot x 1], target id for each MPC row
    vRel = sum((dodVector + doaVector) .* tgtVel(idx,:), 2);  % [Ntot x 1]
    fD = vRel/lambda;
    polarizationPhase = vertcat(targetChannel.InitialPhases);
    dopplerPhase = exp(1j * 2 * pi * fD * timeVector) .* exp(1j*polarizationPhase(:,1));
    tgtPg2pol = (vertcat(targetChannel.AveragePathGains));
    [~, keepMaxPolId] = max(sum(tgtPg2pol));
    tgtPg = 10*log10(tgtPg2pol(:,keepMaxPolId));
    HasLOSCluster = vertcat(targetChannel.HasLOSCluster);
    tgtPower = nrRadar.channel.getTargetPower(targetChannel);

end

% ---------------------------- Background model ----------------------------
if isempty(backgroundChannel)
    aoaAz = [];
    aodAz = [];
    aoaEl = [];
    aodEl = [];
    delays = [];
    phases = [];
    envPgLin = [];
else
    aoaAz = backgroundChannel.AnglesAoA.';
    aodAz = backgroundChannel.AnglesAoD.';
    aoaEl = backgroundChannel.AnglesZoA.';
    aodEl = backgroundChannel.AnglesZoD.';
    delays = backgroundChannel.PathDelays.';
    phases = backgroundChannel.InitialPhases.';
    envPgLin = backgroundChannel.AveragePathGains.';
end

envPg = 10*log10(abs(envPgLin));
keepTime = false(nRealization,1);
keepTime(symbolIndices) = true;
keepPath = (tgtPg > (max(tgtPg) - tgtRayPowerThr));

% Combine (target + background)
path_gain_combined = [tgtPg(keepPath)', envPg];
delays_combined    = [tgtDelay(keepPath)', delays];

aoaAz_combined     = [aoaAzTgt(keepPath)', aoaAz];
aoaEl_combined     = [aoaElTgt(keepPath)', aoaEl];

aodAz_combined     = [aodAzTgt(keepPath)', aodAz];
aodEl_combined     = [aodElTgt(keepPath)', aodEl];

phases_combined = [angle(dopplerPhase(keepPath,:)); repmat(phases, nRealization, 1).'].'; % [T x P]
phases_combined = phases_combined(keepTime,:);
nRealization = size(phases_combined,1);
% Sort by delay
[delays_sorted, sortIdx] = sort(delays_combined);

path_gain_sorted = path_gain_combined(sortIdx);
aoaAz_sorted = aoaAz_combined(sortIdx);
aoaEl_sorted = aoaEl_combined(sortIdx);
aodAz_sorted = aodAz_combined(sortIdx);
aodEl_sorted = aodEl_combined(sortIdx);
phases_sorted = phases_combined(:,sortIdx);

% Build per-path CIR (uninterpolated)
nSamples = ceil((delays_sorted(end) - delays_sorted(1))*sampleRate) + 10;
timeSampling = (0:nSamples-1).'/sampleRate;

cirP = (sqrt(10.^(path_gain_sorted/10)) .* exp(1j*phases_sorted)).';  % [P x T]

% Common pack (reused across TRPs)
common = struct();
common.c = c;
common.fc = fc;
common.sampleRate = sampleRate;
common.bandwidth = bandwidth;
common.nRealization = size(phases_combined,1);
common.nSamples = nSamples;
common.timeSampling = timeSampling;
common.delays_sorted = delays_sorted;
common.cirP = cirP; % [P x T]

common.aoaAz = aoaAz_sorted; common.aoaEl = aoaEl_sorted;
common.aodAz = aodAz_sorted; common.aodEl = aodEl_sorted;

common.aoaAzLOS = aoaAzLOS; common.aoaElLOS = aoaElLOS;
common.aoaAzTgt = aoaAzTgt; common.aoaElTgt = aoaElTgt;
common.aodAzTgt = aodAzTgt; common.aodElTgt = aodElTgt;

common.HasLOSCluster = HasLOSCluster;


% ---------------------------- Apply per-TRP BF ----------------------------

nDigitalChains = receiveArrayConfig.meta.array.Mprime*receiveArrayConfig.meta.array.Nprime;

outCirAntennaPort = zeros(nSamples, nRealization, nDigitalChains, NtrpSel);
out = zeros(NtrpSel,2);

for ii = 1:NtrpSel
    k = trpSelect(ii);
    yawDeg = trpYawDegAll(k);

    [outCirAntennaPort(:,:,:,ii), out(ii,:), rxArrayObj] = applyTrpBeamforming( ...
        common, yawDeg, ...
        'transmitArray', transmitArray, ...
        'receiveArray', receiveArrayConfig, ...
        'angleEstimation', angleEstimation);
end

cirCropStart = max(0,floor((min(tgtDelay(keepPath))-delays_sorted(1))*sampleRate)-5);
cirCropEnd = min(nSamples-1,floor((max(tgtDelay(keepPath))-delays_sorted(1))*sampleRate)+5);
lenCir = cirCropEnd-cirCropStart+1;
% Crop in case of very long CIR
if lenCir >1024
    cirCropEnd = cirCropStart+1024-1;
end
outCirAntennaPort([1:cirCropStart-1, cirCropEnd+1:end],:,:,:) = [];
syncOffset = (delays_sorted(1)+cirCropStart/sampleRate)*c/2; % common


% -------------------------------------------------------------------------
% ---------------------------- Local helpers ------------------------------
% -------------------------------------------------------------------------

% function Xk = pickTRP(X, k)
%     % Allows passing either:
%     %   - [] (empty)
%     %   - numeric matrix (shared across TRPs)
%     %   - cell array {Ntrp} with per-TRP entries
%     if isempty(X)
%         Xk = [];
%         return;
%     end
%     if iscell(X)
%         if k > numel(X) || isempty(X{k})
%             Xk = [];
%         else
%             Xk = X{k};
%         end
%     else
%         Xk = X;
%     end
% end

end % end main function


function [outCirAntennaPort, angleEstimate, receiveArray] = applyTrpBeamforming(common, yawDeg, varargin)
%APPLYTRPBEAMFORMING TRP-specific yaw rotation and beamforming on a common MPC set.
%   OUTCIRANTENNAPORT = APPLYTRPBEAMFORMING(COMMON, YAWDEG) rotates the global azimuth
%   angles in COMMON into the local TRP array frame by applying a yaw
%   rotation of YAWDEG degrees about the z-axis, then applies transmit and
%   receive beamforming to the per-path CIR. It returns the beamformed 
%   channel impulse response (CIR) at the antenna port / RF-chain outputs 
%   as a sampled CIR cube.
%
%   OUTCIRANTENNAPORT = APPLYTRPBEAMFORMING(..., 'transmitArray', TXCFG)
%   specifies the transmit array configuration TXCFG. The function builds a
%   phased.URA using:
%     - TXCFG.meta.array.M, TXCFG.meta.array.N
%     - TXCFG.meta.array.dV_lambda, TXCFG.meta.array.dH_lambda
%   and uses TXCFG.beamformer.wElem as element-domain TX beamforming weights.
%
%   OUTCIRANTENNAPORT = APPLYTRPBEAMFORMING(..., 'receiveArray', RXCFG)
%   specifies the receive array configuration RXCFG. RXCFG.meta.name selects
%   the receive architecture:
%     - 'analog'      : single RF chain (NRF = 1) with analog combining
%     - 'hybrid'      : replicated subarray / subarray steering (NRF > 1)
%     - 'full-digital': per-element receive steering then optional RF-chain
%                       reordering using RXCFG.beamformer.F_RF
%
%   OUTCIRANTENNAPORT = APPLYTRPBEAMFORMING(..., 'angleEstimation', METHOD)
%   sets the method used to estimate the RX angle when RXCFG.meta.name is
%   'analog'. METHOD is passed to estimateRxAngle. Default is 'ideal'.
%
%   Inputs
%     COMMON  - Structure containing MPC/CIR information and sampling parameters:
%              * fc             : Carrier frequency (Hz)
%              * aodAz, aodEl    : AoD azimuth/elevation angles (deg)
%              * aoaAz, aoaEl    : AoA azimuth/elevation angles (deg)
%              * cirP            : Element-domain received CIR per path [P x T]
%              * delays_sorted   : Path delays (seconds)
%              * nSamples        : Number of time samples in output CIR
%              * nRealization    : Number of realizations/snapshots (T)
%              * timeSampling    : Sample period (seconds)
%              * bandwidth       : (Optional) bandwidth for sinc interpolation
%              * HasLOSCluster   : LOS flag(s)
%              * aoaAzLOS, aoaElLOS : LOS angles used for analog angle estimation
%              * aod/aoa*Tgt fields as passed through to OUT
%     YAWDEG  - TRP yaw rotation in degrees (positive rotation about z-axis).
%
%   Output
%     OUTCIRANTENNAPORT - Sampled, beamformed CIR at antenna port / RF chains:
%                         size [nSamples x nRealization x NRF], where NRF is
%                         the number of RF chains (or 1 for analog).
%
%
%   See also phased.URA, phased.ReplicatedSubarray, phased.SteeringVector,
%   estimateRxAngle.
%
%   2026 NIST/CTL Steve Blandino
%
%   This file is available under the terms of the NIST License.

%% Input processing
p = inputParser;
addParameter(p, 'transmitArray', []);
addParameter(p, 'receiveArray', []);
addParameter(p, 'angleEstimation', 'ideal');
parse(p, varargin{:});

transmitArrayConfig   = p.Results.transmitArray;
receiveArrayConfig    = p.Results.receiveArray;
angleEstimation = p.Results.angleEstimation;
angleEstimate = nan(1,2);
fcHz = common.fc;
lambda = 299792458/fcHz;

% Rotate GLOBAL azimuth angles into TRP LOCAL array frame (yaw about z)
aodAz_loc = mod(common.aodAz - yawDeg + 180, 360) - 180;
aoaAz_loc = mod(common.aoaAz - yawDeg + 180, 360) - 180;
aodEl_loc = common.aodEl;
aoaEl_loc = common.aoaEl;

%% Define antenna arrays
transmitArray = phased.URA([transmitArrayConfig.meta.array.M transmitArrayConfig.meta.array.N], ...
    'ElementSpacing',[lambda*transmitArrayConfig.meta.array.dV_lambda,lambda*transmitArrayConfig.meta.array.dH_lambda]);
transmitArray.Element = phased.NRAntennaElement;
txSV = phased.SteeringVector('SensorArray', transmitArray, 'IncludeElementResponse', true);
txPV = txSV(fcHz, [aodAz_loc; aodEl_loc]); % [Ntx x P]

rxMode = receiveArrayConfig.meta.name;
isHybrid = strcmp(rxMode,'hybrid');
isAnalog = strcmp(rxMode,'analog');
isFullDigital = strcmp(rxMode,'full_digital');

if isHybrid
    numElemPerSub = receiveArrayConfig.meta.array.M/receiveArrayConfig.meta.array.Mprime;
    sULA = phased.ULA('NumElements', numElemPerSub, ...
        'Element', phased.NRAntennaElement, ...
        'ElementSpacing', receiveArrayConfig.meta.array.dV_lambda*lambda, ...
        'ArrayAxis', 'z');

    receiveArray = phased.ReplicatedSubarray('Subarray', sULA, ...
        'GridSize', [receiveArrayConfig.meta.array.Mprime receiveArrayConfig.meta.array.Nprime], ...
        'SubarraySteering','Custom', ...
        'GridSpacing', [numElemPerSub*receiveArrayConfig.meta.array.dV_lambda receiveArrayConfig.meta.array.dH_lambda]*lambda);
    rxSV = phased.SteeringVector('SensorArray', receiveArray, 'IncludeElementResponse', true);
    rxPV = rxSV(fcHz, [aoaAz_loc; aoaEl_loc], reshape(receiveArrayConfig.beamformer.wElem, ...
        numElemPerSub,[])); % [Nrx x P]
else
    receiveArray = phased.URA([receiveArrayConfig.meta.array.M receiveArrayConfig.meta.array.N], ...
        'ElementSpacing',[lambda*receiveArrayConfig.meta.array.dV_lambda,lambda*receiveArrayConfig.meta.array.dH_lambda]);
    receiveArray.Element = phased.NRAntennaElement;
    rxSV = phased.SteeringVector('SensorArray', receiveArray, 'IncludeElementResponse', true);
    rxPV = rxSV(fcHz, [aoaAz_loc; aoaEl_loc]); % [Nrx x P]
end

%% Apply TX Pattern
% TX beamforming weigths
w_tx = transmitArrayConfig.beamformer.wElem;
txGain = (w_tx' * txPV).';   % [P x 1]

%% GET CIR
% Element-domain received CIR per path
% common.cirP: [P x T]
if isAnalog
    angleEstimate = estimateRxAngle(transmitArray,angleEstimation, fcHz, txPV, ...
        common.cirP, common.aoaAzLOS- yawDeg, common.aoaElLOS, ...
        'applyFrontBackMask', true, 'boresightAzEl', [0 0]);
    W_RF = sum(conj(rxSV(fcHz, angleEstimate')),2)/numel(common.aoaAzLOS);
    rxGain =  rxPV.' * W_RF;    % [P x 1]
    cirAntennaPortP = (common.cirP .* txGain) .* rxGain;  % [P x T]  (NRF=1)
    % [(P*T) x NRF]
else
    cirPT = common.cirP .* txGain;                 % [P x T]
    cirAntennaPortP = cirPT .* reshape(rxPV.', [size(cirPT,1), 1, size(rxPV,1)]);  % [P x T x NRF]

end

% Interpolate to sampled CIR if bandwidth is set
nSamples = common.nSamples;
nRealization = common.nRealization;
NRF = size(cirAntennaPortP,3);
outCirAntennaPort = zeros(nSamples, nRealization, NRF);

if ~isempty(common.bandwidth)
    % outCir = sincInterp(common.delays_sorted - common.delays_sorted(1), cirBeamP, common.timeSampling, common.bandwidth);
    for i = 1:NRF
        outCirAntennaPort(:,:,i) = nrRadar.dsp.sincInterp(common.delays_sorted - common.delays_sorted(1), ...
            cirAntennaPortP(:,:,i), common.timeSampling, common.bandwidth);
    end
else
    % If no interpolation requested, keep a simple mapping:
    % outCir is returned as [P x T] not matching [nSamples x T].
    % To keep output signature stable, we place impulses on nearest taps.
    % outCir = placeOnTaps(common.delays_sorted - common.delays_sorted(1), cirBeamP, common.timeSampling);
    for i = 1:NRF
        outCirAntennaPort(:,:,i) = placeOnTaps(common.delays_sorted - common.delays_sorted(1), ...
            cirAntennaPortP(:,:,i), common.timeSampling);
    end
end

if isFullDigital
    idxVec = nrRadar.array.getRfChainSortedIndex(receiveArray, receiveArrayConfig.beamformer.F_RF);
    outCirAntennaPort = outCirAntennaPort(:,:,idxVec);
end

end

function out = placeOnTaps(delays, cirP, timeSampling)
% Helper to return [nSamples x T] even when no interpolation is requested.
% Places each path on nearest tap.
%
% delays: [P x 1] seconds relative to first path
% cirP:   [P x T] complex
% timeSampling: [nSamples x 1]

nSamples = numel(timeSampling);
T = size(cirP,2);
out = zeros(nSamples, T);

Ts = timeSampling(2) - timeSampling(1);
tap = round(delays(:)/Ts) + 1;
tap = max(1, min(nSamples, tap));

for p = 1:numel(tap)
    out(tap(p), :) = out(tap(p), :) + cirP(p, :);
end
end
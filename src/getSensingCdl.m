function [out, outCir,tgtPg, syncOffset] = getSensingCdl(bsPos,tgtPos,tgtVel, fc, pri, varargin)
% GETSENSINGCDL Generate a channel delay line (CDL) model with a target
%
%   [OUT, OUTCIR] = GETSENSINGCDL(BSPOS, TGTPOS, TGTVEL, FC, PRI) computes
%   the output structure OUT and the time-varying channel impulse response
%   OUTCIR for a target moving with velocity TGTVEL and position TGTPOS
%   relative to a base station at position BSPOS.
%
%   BSPOS     - 3-element vector specifying the [x, y, z] coordinates of
%               the base station (in meters)
%   TGTPOS    - 3-element vector specifying the initial [x, y, z]
%               coordinates of the target (in meters)
%   TGTVEL    - 3-element vector specifying the [vx, vy, vz] velocity of
%               the target (in meters per second)
%   FC        - Carrier frequency (in Hz)
%   PRI       - Pulse repetition interval (in seconds)
%   TRANSMITARRAY     - phased array system object (phased.ULA, phased.URA..)
%   RECEIVEARRAY     - phased array system object (phased.ULA, phased.URA..)
%
%   Optional Name-Value Pair Arguments:
%   'bandwidth'      - Sampling rate of the channel (in Hz)
%                      Default: 100e6
%   'bistaticAngle'  - Initial bistatic angle of the target (in degrees)
%                      Default: 0
%   'aspectAngle'    - Initial aspect angle of the target (in degrees)
%                      Default: 0
%   'nRcsSamples'    - Number of RCS samples to generate
%                      Default: 1
%   'nRealization'   - Number of realizations for the channel model
%                      Default: 100
%   'angleEstimation' - Angle estimation method between ideal (pointing to
%   target), nearest (wrt codebook), scan (scan a codebook)
%
%   Outputs:
%   OUT       - Struct containing the sorted MPC descriptors:
%               .delay    : Sorted delays (in seconds)
%               .pathGain : Time-varying path gains (complex)
%               .aoaAz    : Azimuth angles of arrival (in degrees)
%               .aoaEl    : Elevation angles of arrival (in degrees)
%               .aodAz    : Azimuth angles of departure (in degrees)
%               .aodEl    : Elevation angles of departure (in degrees)
%
%   OUTCIR    - Time-varying channel impulse response (complex matrix)
%
%   Example:
%   bsPos = [0 0 10];
%   tgtPos = [100 0 0];
%   tgtVel = [0 0 0];
%   fc = 3.5e9; pri = 0.01;
%   [out, outCir] = getSensingCdl(bsPos, tgtPos, tgtVel, fc, pri);
%
%   See also nrCDLChannel, info, sincInterp

% Create input parser
p = inputParser;

validAngleEstimation = {'ideal', 'nearest', 'scan'};


% Add optional parameters with default values
addParameter(p, 'bandwidth', []);        % Default bistatic angle
addParameter(p, 'bistaticAngle', 0);        % Default bistatic angle
addParameter(p, 'aspectAngle', 0);          % Default aspect angle
addParameter(p, 'nRcsSamples', 1);          % Default number of RCS samples
addParameter(p, 'nRealization', 1);       % Default number of realizations
addParameter(p, 'transmitArray', phased.ULA(64));
addParameter(p, 'receiveArray', phased.ULA(64));
addParameter(p, 'firstScatterRangeMax', 30);
addParameter(p, 'scanvector', 0);
addParameter(p, 'scenario', 'UMi-AV');
addParameter(p, 'angleEstimation', 'ideal', @(x) any(validatestring(x, validAngleEstimation)));
addParameter(p, 'backgroundChannel', []);
addParameter(p, 'targetChannel', []);


% Parse inputs
parse(p, varargin{:});
nRealization = p.Results.nRealization;
bandwidth = p.Results.bandwidth;
transmitArray = p.Results.transmitArray;
receiveArray = p.Results.receiveArray;
scanvector = p.Results.scanvector;
scenario = p.Results.scenario;
angleEstimation = p.Results.angleEstimation;
backgroundChannel = p.Results.backgroundChannel;
targetChannel = p.Results.targetChannel;

%% Parameters
c = 299702547;
lambda = c/fc;
if ~isempty(bandwidth)
    sampleRate = bandwidth;
else
    sampleRate = 100e6;
end
timeVector = 0:pri:pri*(nRealization-1);

if isempty(targetChannel)
    tgtRcs =  getSigmaRCS('uav-small', 'returnLargeScale',1);
    % UAV delay and Doppler
    tgtRange = norm(tgtPos - bsPos);
    tgtDelay = 2*tgtRange/c;
    bs2tgtVec = (tgtPos - bsPos)/tgtRange;
    vRel = dot(tgtVel, bs2tgtVec);
    fD = 2*(vRel/lambda);
    [PG, HasLOSCluster] = computePathLoss(scenario, fc, bsPos, tgtPos, 0);
    tgtPg =  -(2*PG- tgtRcs + 10*log10(c^2/(4*pi*fc^2)));

    dopplerPhase = exp(1j * 2 * pi * fD * timeVector);  % Doppler phase shift

    % Extract components of the vector
    vx = bs2tgtVec(1);
    vy = bs2tgtVec(2);
    vz = bs2tgtVec(3);

    % Compute Azimuth AoA (in degrees)
    aoaAzTgt = atan2d(vy, vx);

    % Compute Elevation AoA (in degrees)
    aoaElTgt = atan2d(vz, sqrt(vx^2 + vy^2));

    % Since it's a monostatic system, AoD is the same as AoA
    aodAzTgt = aoaAzTgt;
    aodElTgt = aoaElTgt;
else

    tgtDelay = targetChannel.PathDelays;
    aoaAzTgt = wrapTo180(targetChannel.AnglesAoA);

    % Compute Elevation AoA (in degrees)
    aoaElTgt = 90-targetChannel.AnglesZoA;

    % Since it's a monostatic system, AoD is the same as AoA
    aodAzTgt = wrapTo180(targetChannel.AnglesAoD);
    aodElTgt = 90-targetChannel.AnglesZoD;

    dodVector = angle2vector(aodAzTgt, 90-aodElTgt, 1)./vecnorm(angle2vector(aodAzTgt, 90-aodElTgt, 1),2,2);
    doaVector = angle2vector(aoaAzTgt, 90-aoaElTgt, 1)./vecnorm(angle2vector(aoaAzTgt, 90-aoaElTgt, 1),2,2);
    vRel = dodVector * tgtVel'+ doaVector * tgtVel';

    fD = vRel/lambda;

    dopplerPhase = exp(1j * 2 * pi * fD * timeVector).* exp(1j*targetChannel.InitialPhases(:,1));  
    tgtPg = 10*log10(targetChannel.AveragePathGains);
    HasLOSCluster = targetChannel.HasLOSCluster;
end

%% CDL Channel Configuration
aoaAz =  backgroundChannel.AnglesAoA.';
aodAz =  backgroundChannel.AnglesAoD.';
aoaEl =  backgroundChannel.AnglesZoA.';
aodEl = backgroundChannel.AnglesZoD.';
delays = backgroundChannel.PathDelays.';
phases = backgroundChannel.InitialPhases.';
envPgLin = backgroundChannel.AveragePathGains.';

if isempty(delays)
    aoaAz =  0;
    aodAz = 0;
    aoaEl =  0;
    aodEl =0;
    delays = tgtDelay;
    phases = 0;
end

%% Combine MPC descriptors (UAV and CDL)
% Path gain for UAV and CDL
envPg = 10*log10(abs(envPgLin));
path_gain_combined = [tgtPg', envPg];

% Delays for UAV and CDL
delays_combined = [tgtDelay', delays];

% AoA (Azimuth and Elevation) for UAV and CDL
aoaAz_combined = [aoaAzTgt', aoaAz];
aoaEl_combined = [aoaElTgt', aoaEl];

% AoD (Azimuth and Elevation) for UAV and CDL
aodAz_combined = [aodAzTgt', aodAz];
aodEl_combined = [aodElTgt', aodEl];

% Phases for UAV and CDL (already computed)
phases_combined = [angle(dopplerPhase); repmat(phases, nRealization,1).'].';

%% Sort combined MPC descriptors by delay
[delays_sorted, sortIdx] = sort(delays_combined);  % Sort delays and get indices

% Apply sorting to all other descriptors
path_gain_sorted = path_gain_combined(sortIdx);
aoaAz_sorted = aoaAz_combined(sortIdx);
aoaEl_sorted = aoaEl_combined(sortIdx);
aodAz_sorted = aodAz_combined(sortIdx);
aodEl_sorted = aodEl_combined(sortIdx);
phases_sorted = phases_combined(:,sortIdx);

nSamples = ceil((delays_sorted(end)-delays_sorted(1))*sampleRate)+10;
timeSampling = (0:nSamples-1).'/sampleRate;
cir = (sqrt(10.^(path_gain_sorted/10)).*exp(1j*phases_sorted)).';
dropRayId = (aodEl_sorted>=90 | aoaEl_sorted>=90 | aodEl_sorted<=-90 | aoaEl_sorted<=-90 | path_gain_sorted<(max(tgtPg)-60));
cir(dropRayId,:) = [];
aodAz_sorted(dropRayId) = [];
aodEl_sorted(dropRayId) = [];
aoaAz_sorted(dropRayId) = [];
aoaEl_sorted(dropRayId) = [];
delays_sorted(dropRayId) = [];
% Beamforming
txSteeringVector = phased.SteeringVector('SensorArray', transmitArray);
rxSteeringVector = phased.SteeringVector('SensorArray', receiveArray);

txPV = txSteeringVector(fc, [aodAz_sorted; aodEl_sorted]);
rxPV = rxSteeringVector(fc, [aoaAz_sorted; aoaEl_sorted]);

switch angleEstimation

    case 'scan'
        rxBF = conj(rxSteeringVector(fc, scanvector.'));
        rxGain = rxBF.' * rxPV;
        rxGain = permute(rxGain, [3,2,1]);
        cirTg= cir-mean(cir,2);
        cirTg = repmat(cirTg.', [1 1 size(scanvector,1)]).*rxGain;
        rxPowerPerBeam = squeeze(sum(sum(abs(cirTg),1),2));
        if abs(aoaAzTgt - 90) < 1e-6
            rxPowerPerBeam(sind(scanvector(:,1))<0) = 0;
        elseif abs(aoaAzTgt + 90) < 1e-6
            rxPowerPerBeam(sind(azimuthRange)>0) = 0;
        else
            if cosd(aoaAzTgt) > eps
                rxPowerPerBeam(cosd(scanvector(:,1))<0) = 0;
            else
                rxPowerPerBeam(cosd(scanvector(:,1))>0) = 0;
            end
        end
        [~,idmax] = max(rxPowerPerBeam);
        angleEstimate = scanvector(idmax,:);

    case 'nearest'
        [~, idmax ] = min(sum(abs(scanvector - [aoaAzTgt(1), aoaElTgt(1)]),2));
        angleEstimate = scanvector(idmax,:);

    case 'ideal'
        angleEstimate =  [aoaAzTgt(1), aoaElTgt(1)];
end

txBF = conj(txSteeringVector(fc, angleEstimate'));
rxBF = conj(rxSteeringVector(fc, angleEstimate'));
txGain = txPV.' * txBF;
rxGain = rxBF.' * rxPV;
cir = (cir.'.*txGain.'.*rxGain).';

outCir = cir;
if ~isempty(bandwidth)
    cirInt = sincInterp(delays_sorted-delays_sorted(1), outCir, timeSampling, bandwidth);
    outCir = cirInt;
end

% figure, stem(delays, envPg, 'BaseValue', -250), hold on, stem(tgtDelay, tgtPg, 'BaseValue', -250)
% hold on, plot(timeSampling+delays_sorted(1),20*log10(abs(outCir(:,1))))
% grid on
% ylabel('Power (dB)')
% legend('Background', 'Target', 'CIR')
% xlabel('Delay (s)')
% figure, plot(aodAzTgt, aodElTgt,'o'), hold on, plot(aodAz, aodEl,'o'),
% plot(angleEstimate(1), angleEstimate(2), '*', 'MarkerSize',12)
% legend('Target', 'Background', 'Beamforming'), xlabel('Azimuth (deg)'), ylabel('Elevation (deg)')

%% Output
out.delay = delays_sorted;
out.pathGain = outCir;
out.aoaAz = aoaAz_sorted;
out.aoaEl = aoaEl_sorted;
out.aodAz = aodAz_sorted;
out.aodEl = aodEl_sorted;
out.aodAzTgt = aodAzTgt;
out.aoaAzTgt = aoaAzTgt;
out.aodElTgt = aodElTgt;
out.aoaElTgt = aoaElTgt;
out.angleEstimate = angleEstimate;
out.los = HasLOSCluster;
syncOffset = delays_sorted(1)*c/2;

end
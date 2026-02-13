function [simConfig, prsConfig, sensConfig, geometry] = makeTinyConfig(rxMode)
%MAKETINYCONFIG Create a small, deterministic configuration for unit tests.
%
%   [SIM, PRS, SENS, GEOM] = nrRadarTest.makeTinyConfig(RXMODE)
%   returns minimal structs that satisfy nrRadar.internal.precompute and the
%   downstream RX/Sensing pipeline.
%
%   RXMODE: 'full_digital' | 'hybrid'
%
% The intent is fast execution (seconds), determinism, and stable contracts
% for backwards-compatibility checks.

rxMode = validatestring(rxMode, {'full_digital','hybrid'});

%% Simulation config (minimal fields consumed by precompute and detection)

simConfig = struct();

% Carrier/system
simConfig.systemFc = 4e9;
simConfig.systemBw = 100e6;
simConfig.systemNF = 5;
simConfig.antennaCouplingEfficiency = 1;
simConfig.carrierSubcarrierSpacing = 30;      % kHz
simConfig.carrierNSizeGrid = 273;              % PRBs
simConfig.maxRangeInterest = 400;             % m (test)
simConfig.trpYawDeg = 30;
simConfig.txPower  = 52;
% RX antenna model (minimal schema used by precompute + rdmDetection)
M = 8; N = 8;
dV = 0.5; dH = 0.5;

rxAntenna = struct();
rxAntenna.meta = struct();
rxAntenna.meta.name = rxMode;
rxAntenna.meta.fc_Hz = simConfig.systemFc;
rxAntenna.meta.array = struct('M', M, 'N', N, 'dV_lambda', dV, 'dH_lambda', dH);

switch rxMode
    case 'full_digital'
        rxAntenna.meta.array.Mprime = M;
        rxAntenna.meta.array.Nprime = N;
        rxAntenna.beamformer = struct('wElem', ones(M*N,1));
    case 'hybrid'
        % 2x2 subarray grid, each subarray has M/Mprime elements along z.
        rxAntenna.meta.array.Mprime = 4;
        rxAntenna.meta.array.Nprime = 8;
        numElemPerSub = M / rxAntenna.meta.array.Mprime; % =2
        numSubarrays  = rxAntenna.meta.array.Mprime * rxAntenna.meta.array.Nprime; % =4
        % Uniform subarray weights (normalized)
        rxAntenna.beamformer = struct('wElem', ones(numElemPerSub,numSubarrays) / sqrt(numElemPerSub));
end

% TX antenna is only used in channel generation modules; keep a stub.
simConfig.txAntenna = struct('meta', struct('fc_Hz', simConfig.systemFc));
simConfig.rxAntenna = rxAntenna;

%% PRS config (minimal, valid NR PRS allocation)
prsConfig = struct();
prsConfig.PRSResourceSetPeriod  = [1 0];
prsConfig.PRSResourceOffset     = 0;
prsConfig.PRSResourceRepetition = 1;
prsConfig.PRSResourceTimeGap    = 1;
prsConfig.NumRB                 = simConfig.carrierNSizeGrid;
prsConfig.RBOffset              = 0;
prsConfig.CombSize              = 2;
prsConfig.REOffset              = 0;
prsConfig.NPRSID                = 1;
prsConfig.NumPRSSymbols         = 2;
prsConfig.SymbolStart           = 1;

%% Sensing config
sensConfig = struct();

% Keep CPI short to run quickly.
sensConfig.numberSensingSymbols = 8;
sensConfig.dopplerFftLen        = 8;

% Angle FFTs only matter for 'beamspaceFFT' suppression and DOA grid.
sensConfig.azFftLen = 32;
sensConfig.elFftLen = 32;

sensConfig.window        = 'blackmanharris';
sensConfig.windowLen     = sensConfig.dopplerFftLen;
sensConfig.windowOverlap = 0;

% CFAR (small windows for speed)
sensConfig.cfarGrdCellRange    = 1;
sensConfig.cfarGrdCellVelocity = 1;
sensConfig.cfarTrnCellRange    = 4;
sensConfig.cfarTrnCellVelocity = 2;
sensConfig.cfarThreshold       = 12;  % dB

% Pre-thresholding (keep permissive)
sensConfig.rdaThreshold = 0;          % dB above median

% Peak picking / clustering
sensConfig.nmsRadius    = [2 1 1 1];
sensConfig.nmsMaxPeaks  = 50;
sensConfig.dbscanMinMaxRatio = inf;

% DOA
switch rxMode
    case 'full_digital'
        sensConfig.doaEstimationMethod = 'beamspaceFFT';
    case 'hybrid'
        sensConfig.doaEstimationMethod = 'barlettScan';
end

%% Geometry
geometry = struct();
geometry.tx = [0 0 1.5];
geometry.rx = [0 0 1.5];

end

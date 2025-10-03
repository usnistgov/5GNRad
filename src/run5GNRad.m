function results = run5GNRad(simConfig, stConfig, prsConfig, geometry, sensConfig,backgroundChannel,targetChannel)
% RUN5GNRAD Run PRS-based radar simulation using 5G NR signals
%   RESULTS = RUN5GNRAD(SIMCONFIG, STCONFIG, PRSCONFIG, GEOMETRY,
%   SENSCONFIG, BACKGROUNDCHANNEL) simulates a monostatic
%   radar system based on 5G NR Positioning Reference Signals (PRS). It
%   computes a range-Doppler map and estimates the position and velocity of
%   the target across a number of sensing snapshots.
%
%   Inputs:
%     SIMCONFIG          - Structure with system-level configuration (e.g., fc, bandwidth, noise figure)
%     STCONFIG           - Structure with true target states:
%                          * position: [N x 3] target positions
%                          * velocity: [N x 3] target velocities
%     PRSCONFIG          - PRS configuration structure for nrPRSConfig
%     GEOMETRY           - Structure with TX and RX geometry (fields: tx, rx)
%     SENSCONFIG         - Sensing parameters (e.g., window type, Doppler FFT size)
%     BACKGROUNDCHANNEL  - Cell array or vector with background channel impulse responses per snapshot
%
%   Output:
%     RESULTS - Structure containing per-snapshot radar performance metrics:
%         * positionError        - Euclidean error in 3D position estimate
%         * rangeError           - Error in estimated range
%         * velocityError        - Error in estimated radial velocity
%         * azimuthError         - AOA azimuth error
%         * elevationError       - AOA elevation error
%         * dopplerPeakPower     - Peak value of Doppler map
%         * dopplerpeakToAverage - Peak-to-average power ratio in Doppler domain
%         * isLos                - Boolean LOS flag from channel model
%
%   The simulation includes:
%     - PRS generation and OFDM mapping
%     - CDL-based sensing channel modeling
%     - Range-Doppler processing and geometry estimation
%
%   Example:
%     results = RUN5GNRAD(simCfg, targetCfg, prsCfg, geometry, sensCfg, bgChannel);

%   2025 NIST/CTL Steve Blandino

%   This file is available under the terms of the NIST License.


%% CONSTANTS
c = 299702547;
subcarrierRB = 12;
maxPowerFcc_100MHz = 75; % https://www.ecfr.gov/current/title-47/chapter-I/subchapter-B/part-30/subpart-C
k = 1.380649*1e-23; 
T = 297;
%% SYSTEM PARAMS
fc = simConfig.systemFc; % Carrier frequency
numberSensingSymbols = sensConfig.numberSensingSymbols;
dopplerFftLen = sensConfig.dopplerFftLen;

% Carrier
carrier = nrCarrierConfig;
carrier.SubcarrierSpacing = simConfig.carrierSubcarrierSpacing;
carrier.NSizeGrid = simConfig.carrierNSizeGrid;

%% Configure PRS 
prs = nrPRSConfig;
prs.PRSResourceSetPeriod = prsConfig.PRSResourceSetPeriod;
prs.PRSResourceOffset = prsConfig.PRSResourceOffset;
prs.PRSResourceRepetition = prsConfig.PRSResourceRepetition;
prs.PRSResourceTimeGap = prsConfig.PRSResourceTimeGap;
prs.NumRB = prsConfig.NumRB;
prs.RBOffset = prsConfig.RBOffset;
prs.CombSize = prsConfig.CombSize;
prs.REOffset = prsConfig.REOffset;
prs.NPRSID = prsConfig.NPRSID;
prs.NumPRSSymbols = prsConfig.NumPRSSymbols;
prs.SymbolStart = prsConfig.SymbolStart;
prs.NumRB = carrier.NSizeGrid;

% Get the number of orthogonal frequency division multiplexing (OFDM) symbols per slot.
numSymPerSlot = carrier.SymbolsPerSlot;
numSlots = numberSensingSymbols*prs(1).PRSResourceSetPeriod(1);

%% OFDM PARAMS
ofdmInfo = nrOFDMInfo(carrier);
ofdmSymbolTime = ofdmInfo.SymbolLengths / ofdmInfo.SampleRate;  % Time in seconds
ofdmTs = mean(ofdmSymbolTime);
prsPeriodicity =  sum(ofdmSymbolTime)/carrier.SlotsPerSubframe*prs.PRSResourceSetPeriod(1); 

%% DEPENDENT PARAMS
numberSubcarriers = carrier.NSizeGrid * subcarrierRB; % Total subcarriers 

% Range resolution
prsRangeResolution = 1/(2*ofdmInfo.SampleRate)*c;
rangeFFTLen = ofdmInfo.Nfft;
rangeWindow = getDftWindow('hamming', numberSubcarriers);
rangeBinDestgrd = (0:rangeFFTLen*prs.CombSize-1) * prsRangeResolution * 1/prs.CombSize;

% Velocity Resolution
vosf = dopplerFftLen/numberSensingSymbols;
prsVelocityResolution = c / (2*numberSensingSymbols*prsPeriodicity*fc)/vosf;

% Velocity Bin
velocityBin = (-dopplerFftLen/2:dopplerFftLen/2-1)*prsVelocityResolution;

lambda = c/fc;

%% Antenna params
nAntH = simConfig.antennaNumH;
nAntV = simConfig.antennaNumV;
transmitArray = phased.URA([nAntH nAntV], 'ElementSpacing',lambda/2);
receiveArray = phased.URA([nAntH nAntV], 'ElementSpacing',lambda/2);
hpbwH = 0.886*2/nAntH*180/pi;
hpbwV = 0.886*2/nAntV*180/pi;
scanStepH = floor(hpbwH);
scanStepV = floor(hpbwV);
azimuthRange = -180:scanStepH:180;
elevationRange = -90:scanStepV:90;
[az, el] = meshgrid(azimuthRange, elevationRange);
scanVector = [az(:), el(:)];

%% Power
maxPowerFcc = maxPowerFcc_100MHz - 10*log10(100e6) + 10*log10(simConfig.systemBw);
maxPower = maxPowerFcc - (10*log10(nAntV*nAntH) + 10*log10(simConfig.antennaCouplingEfficiency));

%% SNR DEFINITION
NF = 10^(simConfig.systemNF/10);
N = k*T*simConfig.systemBw*NF;
P = 10^((maxPower)/10);
SNR = 10*log10(P/N);

%% CHANNEL MODEL PARAMETERS
% Generate random channel gains, delays, and Doppler shifts for L objects
txPos = geometry.tx;
rxPos = geometry.rx;
scenario = simConfig.channelScenario;
stPosition = stConfig.position; 
stPositionEstimate = nan(size(stPosition));

%% CHANNEL MODEL DEPENDENT PARAMETERS
velocityVectors = stConfig.velocity;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End Config
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Transmission: Generate OFDM grid
% Extract OFDM parameters
ofdmInfo    = nrOFDMInfo(carrier);
ofdmFftLen  = ofdmInfo.Nfft;
cpLengths   = ofdmInfo.CyclicPrefixLengths;  % length for each symbol in a slot

% Map the resource elements (RE) for both of the PRS resources on the carrier resource grid.
% Determine the start and end indices for the active subcarriers
startIdx = (ofdmFftLen - numberSubcarriers) / 2 + 1;
endIdx = startIdx + numberSubcarriers - 1;

maxElementsPerSlot = ofdmFftLen * numSymPerSlot;

totalMaxElements = numSlots * maxElementsPerSlot;

Iall = zeros(totalMaxElements, 1);
Jall = zeros(totalMaxElements, 1);
Vall = zeros(totalMaxElements, 1);
ptr = 1;

for slotIdx = 0:prs.PRSResourceSetPeriod(1):numSlots-1
    carrier.NSlot = slotIdx;
    indCell = nrPRSIndices(carrier,prs,'OutputResourceFormat','cell');
    symCell = nrPRS(carrier,prs,'OutputResourceFormat','cell');
    idx = indCell{1};   
    val = symCell{1};

    [rowIdx, colIdx] = ind2sub([numberSubcarriers, numSymPerSlot], idx);
    rowIdx = rowIdx+startIdx-1;
    colIdx = colIdx + slotIdx * numSymPerSlot;

    rangeIdx = ptr : ptr + numberSubcarriers - 1;
    Iall(rangeIdx) = rowIdx;
    Jall(rangeIdx) = colIdx;
    Vall(rangeIdx) = val;

    ptr = ptr + numel(val);
end

Iall(ptr:end) = [];
Jall(ptr:end) = [];
Vall(ptr:end) = [];

% sparse assembly
ofdmGrid = sparse(Iall, Jall, Vall, ofdmFftLen, carrier.SymbolsPerSlot*numSlots);


% Get the column indices of non-zero elements in ofdmGrid
[~, colIdx] = find(ofdmGrid);

% Determine the slot indices by dividing column indices by numSymPerSlot
symbolIndices = unique(colIdx);

%% Transmission: generate TX waveform
txWaveform = zeros((cpLengths(2)+ofdmFftLen)*prs.NumPRSSymbols*numberSensingSymbols, 1);
for symIdx = 0:prs.NumPRSSymbols*numberSensingSymbols-1

    % Frequency-domain samples for this symbol
    freqDomainSymbol = sqrt(prs.CombSize)*ofdmGrid(:, symbolIndices(symIdx+1));

    timeDomainNoCP = ifft(full(freqDomainSymbol), ofdmFftLen) * sqrt(ofdmFftLen);

    % Determine CP length for symbol symIdx within the slot:
    cpLen = cpLengths(mod(symbolIndices(symIdx+1)-1, length(cpLengths)) + 1);

    % Extract the CP samples from the end of the symbol
    cyclicPrefix = timeDomainNoCP(end - cpLen + 1 : end);

    % Concatenate CP + time-domain symbol
    timeDomainWithCP = [cyclicPrefix; timeDomainNoCP];

    % Append to the full TX waveform
    txWaveform((cpLengths(2)+ofdmFftLen)*symIdx+1: (cpLengths(2)+ofdmFftLen)*(symIdx+1) ) = timeDomainWithCP;
end

Nst = size(stPosition,1);
velocityError = zeros(1,Nst);
rangeError = zeros(1,Nst);
snrvarVec = zeros(1,Nst);
signalPower = zeros(1,Nst);
peakToAverage = zeros(1,Nst);
positionError = zeros(1,Nst);
azimuthError = zeros(1,Nst);
elevationError= zeros(1,Nst);
isLos= zeros(1,Nst);
tStart = tic;
for q = 1:Nst
    %% Realize Channel

    baseArgs = { ...
        txPos, stPosition(q,:), velocityVectors(q,:), fc, ofdmTs};

    opts = struct( ...
        'bandwidth', ofdmInfo.SampleRate, ...
        'nRealization', symbolIndices(end), ...
        'transmitArray', transmitArray, ...
        'receiveArray', receiveArray, ...
        'scanvector', scanVector, ...
        'scenario', scenario, ...
        'backgroundChannel', backgroundChannel(q),...
        'angleEstimation', sensConfig.angleEstimationMethod);

    if ~isempty(targetChannel)
        opts.targetChannel = targetChannel(q);
    end

    args = [baseArgs, namedargs2cell(opts)];

    % Call function with dynamically built args
    [info,H,tgtPwr,syncOffset] = getSensingCdl(args{:});

    %% Through the channel
    numSamplesTotal = length(txWaveform);

    snrvar = 10^(SNR/10);
    samplePointer = 0;
    rxWaveform = zeros(numSamplesTotal + size(H,1) - 1, 1);  

    for symIdx = 0:prs.NumPRSSymbols*numberSensingSymbols-1

        % Symbol block in TX waveform
        cpLen = cpLengths(mod(symbolIndices(symIdx+1)-1, length(cpLengths)) + 1);
        symbolLen = ofdmFftLen + cpLen;

        % Extract samples for this symbol from txWaveform
        thisSymbolTx = txWaveform(samplePointer + 1 : samplePointer + symbolLen);

        hSym = H(:,symbolIndices(symIdx+1));

        % Convolution
        thisSymbolRx = conv(thisSymbolTx, hSym);

        % Add the received symbol to the overall rxWaveform in the correct place
        rxWaveform(samplePointer + 1 : samplePointer + length(thisSymbolRx)) = ...
            rxWaveform(samplePointer + 1 : samplePointer + length(thisSymbolRx)) + thisSymbolRx;

        % Update pointer
        samplePointer = samplePointer + symbolLen;
    end
   
    noise = sqrt(1/(2*snrvar)) * (randn(size(rxWaveform)) + 1j *randn(size(rxWaveform)));

    rxWaveform  = rxWaveform+noise;
    snrvarVec(q) = 10*log10(sum(snrvar*10.^(tgtPwr/10)));

    %% Receiver Processing: retrieve OFDM grid

    rows = (1:ofdmFftLen).';
    samplePointer = 0;
    freqDomainSymbolRxStore =  zeros(prs.NumPRSSymbols*numberSensingSymbols*ofdmFftLen, 1);
    linearIdxStore = zeros(prs.NumPRSSymbols*numberSensingSymbols*ofdmFftLen, 1);
    for symIdx = 0:prs.NumPRSSymbols*numberSensingSymbols-1

        cpLen = cpLengths(mod(symbolIndices(symIdx+1)-1, length(cpLengths)) + 1);
        symbolLen = ofdmFftLen + cpLen;

        % Extract the samples for this symbol
        thisSymbolRx = rxWaveform(samplePointer + 1 : samplePointer + symbolLen);

        % Remove CP
        thisSymbolRxNoCP = thisSymbolRx(cpLen + 1 : end);

        % Take FFT (including the sqrt(NFFT) normalization)
        freqDomainSymbolRx = (1 / sqrt(ofdmFftLen)) * fft(thisSymbolRxNoCP, ofdmFftLen);

        % Store into the received OFDM grid
        linearIdx = rows + (symbolIndices(symIdx+1) - 1) * ofdmFftLen;
        freqDomainSymbolRxStore(symIdx*ofdmFftLen+1: (symIdx+1)*ofdmFftLen ) = freqDomainSymbolRx;
        linearIdxStore(symIdx*ofdmFftLen+1: (symIdx+1)*ofdmFftLen ) = linearIdx;

        % Advance pointer
        samplePointer = samplePointer + symbolLen;
    end

    [rowIdx, colIdx] = ind2sub([ofdmFftLen, carrier.SymbolsPerSlot * numSlots], linearIdxStore);
    rxOfdmGrid = sparse(rowIdx, colIdx, freqDomainSymbolRxStore, ofdmFftLen,  carrier.SymbolsPerSlot * numSlots);

    %% Receiver Processing: Channel estimate
    g_tilde_kn = rxOfdmGrid.*conj(ofdmGrid);

    %% Receiver Processing: Sensing processing
    % Range processing
    rangeFft = zeros(rangeFFTLen*prs.CombSize,numberSensingSymbols);

    % Loop over PRS resource symbols in each slot
    for slotIdx = 0:prs.PRSResourceSetPeriod(1):numSlots-1
        carrier.NSlot = slotIdx;
        indCell = nrPRSIndices(carrier,prs,'OutputResourceFormat','cell');
        slotGridRx =  g_tilde_kn(startIdx:endIdx,(1:numSymPerSlot)+numSymPerSlot*slotIdx);
        symCellRx = full(slotGridRx(indCell{1}));
        symIndSlot = symbolIndices(symbolIndices < (slotIdx+1)*ofdmInfo.SymbolsPerSlot &  symbolIndices > (slotIdx)*ofdmInfo.SymbolsPerSlot);

        prsSlotGridRx = reshape(symCellRx, prs.NumRB*12/prs.CombSize, []);
        prsSlotGridRxDestagrd = prsDestaggering(prsSlotGridRx, prs,mod(symIndSlot-1, ofdmInfo.SymbolsPerSlot));

        rangeFft(:, ceil(slotIdx/prs.PRSResourceSetPeriod(1))+1) = sqrt(rangeFFTLen*prs.CombSize)*ifft(prsSlotGridRxDestagrd.*rangeWindow, rangeFFTLen*prs.CombSize);
    end

    rangeFft = rangeFft - mean(rangeFft,2);

    rdMap  = 1/sqrt(dopplerFftLen) * stft(rangeFft, 'window', sensConfig.window,...
        'windowLen', sensConfig.windowLen,...
        'overlap', sensConfig.windowOverlap, ...
        'nfft', dopplerFftLen);

    %% Detection
    [signalPower(q), linearIndex] = max(rdMap(:));
    peakToAverage(q) = 10*log10(signalPower(q)/mean(rdMap(:)));

    % Converting linear index to subscript indices
    [row, col] = ind2sub(size(rdMap), linearIndex);
    d_hat_k = rangeBinDestgrd(row)+syncOffset; 
    v_hat_parallel_k = velocityBin(col); 
    gtVel = sum(((-txPos+stPosition(q,:))/norm(txPos-stPosition(q,:))).* (velocityVectors(q,:)));

    velocityError(q) = abs(gtVel - v_hat_parallel_k);
    rangeError(q) = abs(norm(txPos - stPosition(q,:)) - d_hat_k);
    delay = 2*d_hat_k(:)/c*1e9;

    [~, P1,P2] = estimateScatteringGeometry(txPos ,rxPos, delay,info.angleEstimate);
    [~, im] = min(vecnorm(([P1;P2]- stPosition(q,:))'));
    if im ==1
        stPositionEstimate(q,:) =  P1;
    else
        stPositionEstimate(q,:) =  P2;
    end
    positionError(q) = norm(stPositionEstimate(q,:)-stPosition(q,:));
    azimuthError(q) = info.aoaAzTgt(1) - info.angleEstimate(1);
    elevationError(q) = info.aoaElTgt(1) - (info.angleEstimate(2));
    isLos(q) = info.los;

    percent = q / Nst * 100;
    elapsed = toc(tStart);
    est_total = elapsed / q * Nst;
    eta = est_total - elapsed;

    % Format ETA as duration
    etaDur = duration(0, 0, eta);  % duration in HH:MM:SS
    fprintf('\rProgress: %5.1f%% | ETA: %s', percent, char(etaDur));
    % display(q)
end
% figure, surf(velocityBin,rangeBinDestgrd, 20*log10(rdMap)), xlabel('Velocity (m/s)'), ylabel('Range (m)'), zlabel('Doppler Power Spectrum (dB/Hz)')
% shading flat
results.positionError = positionError(:);
results.rangeError = rangeError(:);
results.velocityError = velocityError(:);
results.azimuthError = azimuthError(:);
results.elevationError = elevationError(:);
results.dopplerPeakPower = signalPower(:);
results.dopplerpeakToAverage = peakToAverage(:);
results.isLos = isLos(:);

end
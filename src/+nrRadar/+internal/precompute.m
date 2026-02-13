function [pre, txWaveform] = precompute(simConfig, prsConfig, sensConfig, geometry)
%RUN5GNRAD_PRECOMPUTE Precompute all drop-invariant quantities for run5GNRad.
%
%   PRE = RUN5GNRAD_PRECOMPUTE(SIMCONFIG, PRSCONFIG, SENSCONFIG, GEOMETRY)
%   builds and returns a struct PRE containing:
%     - Carrier/PRS objects and OFDM info
%     - Range/Doppler/angle grids and windows
%     - PRS OFDM sparse grid and TX waveform
%     - Precomputed PRS indices per slot (to avoid recomputing in each drop)
%     - Convenience constants used by the per-drop worker
%
%   Intended usage:
%     pre = run5GNRad_precompute(...);
%     parfor q = 1:Ndrop
%         drop = run5GNRad_dropWorker(q, pre, ...);
%     end

%% -------------------- Constants --------------------
pre.c  = 299702547;            % (your code's value)
pre.kB = 1.380649e-23;
pre.T0 = 297;
subcarrierRB = 12;

%% -------------------- System params --------------------
pre.fc                    = simConfig.systemFc;
pre.lambda                = pre.c / pre.fc;
pre.numberSensingSymbols  = sensConfig.numberSensingSymbols;
pre.dopplerFftLen         = sensConfig.dopplerFftLen;

% Geometry
pre.txPos = geometry.tx;
pre.rxPos = geometry.rx;

%% -------------------- Carrier --------------------
carrier = nrCarrierConfig;
carrier.SubcarrierSpacing = simConfig.carrierSubcarrierSpacing;
carrier.NSizeGrid         = simConfig.carrierNSizeGrid;
carrier.CyclicPrefix      = 'Normal';
pre.carrier = carrier;

%% -------------------- PRS --------------------
prs = nrPRSConfig;
prs.PRSResourceSetPeriod     = prsConfig.PRSResourceSetPeriod;
prs.PRSResourceOffset        = prsConfig.PRSResourceOffset;
prs.PRSResourceRepetition    = prsConfig.PRSResourceRepetition;
prs.PRSResourceTimeGap       = prsConfig.PRSResourceTimeGap;
prs.NumRB                    = prsConfig.NumRB;           %#ok<NASGU> (kept for compatibility)
prs.RBOffset                 = prsConfig.RBOffset;
prs.CombSize                 = prsConfig.CombSize;
prs.REOffset                 = prsConfig.REOffset;
prs.NPRSID                   = prsConfig.NPRSID;
prs.NumPRSSymbols            = prsConfig.NumPRSSymbols;
prs.SymbolStart              = prsConfig.SymbolStart;

% Your original code overwrote NumRB with carrier.NSizeGrid (keep same behavior)
prs.NumRB = carrier.NSizeGrid;

pre.prs = prs;

pre.numSymPerSlot = carrier.SymbolsPerSlot;
pre.numSlots      = pre.numberSensingSymbols * prs(1).PRSResourceSetPeriod(1);  % matches your code

%% -------------------- OFDM info --------------------
ofdmInfo          = nrOFDMInfo(carrier);
pre.ofdmInfo      = ofdmInfo;
pre.ofdmFftLen    = ofdmInfo.Nfft;
pre.cpLengths     = ofdmInfo.CyclicPrefixLengths;   % per symbol
pre.sampleRate    = ofdmInfo.SampleRate;

% Symbol time / Ts
ofdmSymbolTime    = ofdmInfo.SymbolLengths / ofdmInfo.SampleRate;  % seconds
pre.ofdmTs        = mean(ofdmSymbolTime);

% PRS periodicity as in your code
pre.prsPeriodicity = sum(ofdmSymbolTime)/carrier.SlotsPerSubframe * prs.PRSResourceSetPeriod(1);

% Total symbols on the timeline (used later in channel estimation packing)
pre.nSymTot = carrier.SymbolsPerSlot * pre.numSlots;

%% -------------------- Range/Doppler grids + windows --------------------
pre.numberSubcarriers    = carrier.NSizeGrid * subcarrierRB;
pre.rangeFFTLen          = pre.ofdmFftLen;                    % your code uses Nfft
pre.prsRangeResolution   = 1/(2*pre.sampleRate) * pre.c;
pre.rangeBinDestgrd      = (0:pre.rangeFFTLen-1) * pre.prsRangeResolution;

% Windows (your helper)
pre.rangeWindow          = nrRadar.dsp.getDftWindow('hamming', pre.numberSubcarriers);
dopplerWindow            = nrRadar.dsp.getDftWindow('hamming', pre.numberSensingSymbols);
pre.dopplerWindow        = dopplerWindow;
pre.dopplerWindow3D      = reshape(dopplerWindow(:), 1, [], 1);    % 1 x Nsense x 1

% Velocity bins
vosf = pre.dopplerFftLen / pre.numberSensingSymbols;  % oversampling factor
pre.prsVelocityResolution = pre.c / (2*pre.numberSensingSymbols*pre.prsPeriodicity*pre.fc) / vosf;
pre.velocityBin           = (-pre.dopplerFftLen/2 : pre.dopplerFftLen/2-1) * pre.prsVelocityResolution;

% Optional: range truncation length (used later right after range IFFT)
pre.rangeLength = [];
if isfield(simConfig, 'maxRangeInterest') && ~isempty(simConfig.maxRangeInterest)
    idx = find(pre.rangeBinDestgrd > simConfig.maxRangeInterest, 1, 'first');
    if isempty(idx)
        pre.rangeLength = pre.ofdmFftLen - 1;
    else
        pre.rangeLength = min(pre.ofdmFftLen - 1, idx);
    end
end

%% -------------------- RX antenna / angle FFT grids --------------------
pre.Nfft_h = sensConfig.azFftLen;
pre.Nfft_v = sensConfig.elFftLen;

Mprime = simConfig.rxAntenna.meta.array.Mprime;
Nprime = simConfig.rxAntenna.meta.array.Nprime;
pre.Nprime = Nprime;
pre.Mprime = Mprime;
pre.M = simConfig.rxAntenna.meta.array.M;
pre.N = simConfig.rxAntenna.meta.array.N;
pre.nChanRx = Mprime * Nprime;

% Taylor windows + separable 2D window over array dims (3,4)
wv = taylorwin(Mprime, 4, -35);
wh = taylorwin(Nprime, 4, -35);
pre.W =  wv * wh.'; %reshape(wv, 1,1,Mprime,1) .* reshape(wh, 1,1,1,Nprime);

% Angle grids (your helper call preserved)
dV = simConfig.rxAntenna.meta.array.M / simConfig.rxAntenna.meta.array.Mprime ...
     * pre.lambda * simConfig.rxAntenna.meta.array.dV_lambda;
dH = pre.lambda * simConfig.rxAntenna.meta.array.dH_lambda;

[pre.azGrid, pre.elGrid] = nrRadar.array.getArrayAngleGrid(pre.Nfft_v, pre.Nfft_h, dV, dH, pre.lambda);

%% SNR 
% (Your code defines a constant SNR based on system BW & NF and txPower, then uses tgtPwr separately.)
NF_lin   = 10^(simConfig.systemNF/10);
N_watt   = pre.kB * pre.T0 * simConfig.systemBw * NF_lin;
P_watt   = 10^((simConfig.txPower - 30)/10);     % dBm->dBW->W 
pre.SNR_dB  = 10*log10(P_watt / N_watt);
pre.snrvar  = 10^(pre.SNR_dB/10);                % linear

%% -------------------- Build OFDM sparse grid (PRS mapped) --------------------
% Active subcarrier placement inside Nfft grid
pre.startIdx = (pre.ofdmFftLen - pre.numberSubcarriers) / 2 + 1;
pre.endIdx   = pre.startIdx + pre.numberSubcarriers - 1;

% Slot indices where PRS is present (same stepping you use)
prsPeriodSlots = prs.PRSResourceSetPeriod(1);
slotIdxList    = 0:prsPeriodSlots:pre.numSlots-1;
pre.slotIdxList = slotIdxList;

% Precompute PRS indices per slot (so worker doesn't call nrPRSIndices repeatedly)
% Also build the global sparse grid in one pass.
Icell = cell(numel(slotIdxList),1);
Jcell = cell(numel(slotIdxList),1);
Vcell = cell(numel(slotIdxList),1);

pre.prsIndPerSlot      = cell(pre.numberSensingSymbols, 1);  % indexed by (thisSlot+1)
pre.symIndSlotPerSlot  = cell(pre.numberSensingSymbols, 1);  % indexed by (thisSlot+1)

for s = 1:numel(slotIdxList)
    slotIdx = slotIdxList(s);
    carrier.NSlot = slotIdx;

    indCell = nrPRSIndices(carrier, prs, 'OutputResourceFormat','cell');
    symCell = nrPRS(carrier, prs, 'OutputResourceFormat','cell');

    idx = indCell{1};
    val = symCell{1};

    % idx is in a [NumSubcarriers x NumSymPerSlot] slot grid (your original assumption)
    [rowIdx, colIdx] = ind2sub([pre.numberSubcarriers, pre.numSymPerSlot], idx);

    % Map into full Nfft grid + global timeline symbol indexing
    rowIdx = rowIdx + pre.startIdx - 1;
    colIdx = colIdx + slotIdx * pre.numSymPerSlot;

    Icell{s} = rowIdx(:);
    Jcell{s} = colIdx(:);
    Vcell{s} = val(:);

    % For the range-processing stage: store indCellSlot and symIndSlot per sensed slot
    thisSlot0 = ceil(slotIdx / prsPeriodSlots);   % 0-based, as in your code
    slotKey   = thisSlot0 + 1;                    % MATLAB 1-based

    pre.prsIndPerSlot{slotKey} = indCell{1};

    % symbolIndices are global columns; we can precompute the per-slot symbol columns later,
    % but we also store the raw set once we know symbolIndices (below). For now, store slotIdx.
    pre.symIndSlotPerSlot{slotKey} = slotIdx;     % placeholder; overwritten after symbolIndices computed
end

Iall = vertcat(Icell{:});
Jall = vertcat(Jcell{:});
Vall = vertcat(Vcell{:});

% Global sparse OFDM grid
pre.ofdmGrid = sparse(Iall, Jall, Vall, pre.ofdmFftLen, pre.numSymPerSlot*pre.numSlots);

% Determine which OFDM symbols contain PRS (global indices)
[~, colNZ] = find(pre.ofdmGrid);
pre.symbolIndices = unique(colNZ);

% Now overwrite symIndSlotPerSlot with the actual per-slot PRS symbol indices
for s = 1:numel(slotIdxList)
    slotIdx = slotIdxList(s);
    thisSlot0 = ceil(slotIdx / prsPeriodSlots);
    slotKey   = thisSlot0 + 1;

    symIndSlot = pre.symbolIndices( ...
        pre.symbolIndices < (slotIdx+1)*pre.numSymPerSlot & ...
        pre.symbolIndices > (slotIdx)*pre.numSymPerSlot );

    pre.symIndSlotPerSlot{slotKey} = symIndSlot;
end

%% -------------------- Precompute per-processed-symbol "comb parity" --------------------
% This replaces the repeated:
%   all(mod(find(ofdmGrid(:, symbolIndices(symIdx+1))),2)==0)
nProcSym = prs.NumPRSSymbols * pre.numberSensingSymbols;
pre.nProcSym = nProcSym;

pre.isEvenCombSym = false(nProcSym,1);
for k = 1:nProcSym
    symCol = pre.symbolIndices(k);
    rowsNZ = find(pre.ofdmGrid(:, symCol));
    pre.isEvenCombSym(k) = ~isempty(rowsNZ) && all(mod(rowsNZ,2)==0);
end

%% -------------------- TX waveform generation --------------------
% NOTE: Your original code uses a constant stride (cpLengths(2)+Nfft).
% We keep the same layout for compatibility with your samplePointer arithmetic.
pre.txSymbolStride = (pre.cpLengths(2) + pre.ofdmFftLen);

txWaveform = zeros(pre.txSymbolStride * nProcSym, 1);

samplePtr = 0; %#ok<NASGU> (kept for readability; we index directly like your code)
for symIdx0 = 0:nProcSym-1
    symCol = pre.symbolIndices(symIdx0+1);

    % Frequency-domain samples for this symbol
    freqDomainSymbol = sqrt(prs.CombSize) * pre.ofdmGrid(:, symCol);

    timeDomainNoCP = ifft(full(freqDomainSymbol), pre.ofdmFftLen) * sqrt(pre.ofdmFftLen);

    % CP length for the symbol index within the slot timeline:
    cpLen = pre.cpLengths(mod(symCol-1, numel(pre.cpLengths)) + 1);

    cyclicPrefix = timeDomainNoCP(end-cpLen+1:end);
    timeDomainWithCP = [cyclicPrefix; timeDomainNoCP];

    a = pre.txSymbolStride*symIdx0 + 1;
    b = pre.txSymbolStride*(symIdx0+1);

    % If cpLen varies, timeDomainWithCP might not exactly match stride.
    % To preserve your existing behavior, we zero-pad or truncate to stride.
    L = numel(timeDomainWithCP);
    if L < pre.txSymbolStride
        txWaveform(a:a+L-1) = timeDomainWithCP;
        % remaining samples (if any) stay zero
    else
        txWaveform(a:b) = timeDomainWithCP(1:pre.txSymbolStride);
    end
end

%% Convenience
pre.rows = (1:pre.ofdmFftLen).';

end

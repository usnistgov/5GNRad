function rdCube = getRangeDoppler(h, staticParams)
% GETRANGEDOPPLER Compute range–Doppler cube from PRS RX grids.
%   RDCUBE = GETRANGEDOPPLER(H, STATICPARAMS) performs
%   range processing (IFFT across PRS subcarriers) and Doppler processing
%   (FFT across sensing symbols) to form a range–Doppler cube from the
%   per-antenna received OFDM grids in H.
%
%   The function extracts PRS resource elements for each slot listed in
%   STATICPARAMS.SLOTIDXLIST using precomputed PRS index maps, applies PRS
%   destaggering, performs a windowed range IFFT, removes static clutter by
%   subtracting the mean across the slow-time dimension, optionally
%   truncates range bins, and finally applies a Doppler window and computes
%   a Doppler FFT with FFTSHIFT.
%
%   INPUTS
%     H     : 1×NRX cell array. Each cell contains a complex
%                     matrix of received OFDM grid samples with size
%                     [Nsubcarriers × NsymbolsTotal].
%     STATICPARAMS   : Struct of precomputed/static parameters used across
%                     drops
%
%   OUTPUTS
%     RDCUBE         : Complex range–Doppler cube with size
%                     [Nrng × Ndopp × NRX], where Nrng is RANGEFFTLEN (or
%                     STATICPARAMS.RANGELENGTH if truncation is applied),
%                     and Ndopp is STATICPARAMS.DOPPLERFFTLEN.
%
%
%   2026 NIST/CTL Steve Blandino
%
%   This file is available under the terms of the NIST License.

prs           = staticParams.prs;
rangeFFTLen   = staticParams.rangeFFTLen;
nChanRx       = staticParams.nChanRx;
numSymPerSlot = staticParams.numSymPerSlot;

%% Range Processing (Range FFT)
% Initialize grid: [RangeBins x SensingSymbols x Antennas]
rangeFft = zeros(rangeFFTLen, staticParams.numberSensingSymbols, nChanRx);

% Loop over the precomputed slot list
for s = 1:numel(staticParams.slotIdxList)
    slotIdx = staticParams.slotIdxList(s);

    % Calculate indices
    thisSlot0 = ceil(slotIdx / prs.PRSResourceSetPeriod(1));
    slotKey   = thisSlot0 + 1; % Index for retrieving precomputed maps

    indCellSlot = staticParams.prsIndPerSlot{slotKey};
    symIndSlot  = staticParams.symIndSlotPerSlot{slotKey};
    slotCols    = (1:numSymPerSlot) + numSymPerSlot * slotIdx;

    for nrx = 1:nChanRx
        % Extract slot and PRS symbols
        slotGridRx = h{nrx}(staticParams.startIdx:...
            staticParams.endIdx, slotCols);
        symCellRx  = slotGridRx(indCellSlot);

        % Reshape and Destagger
        numREPerSym = prs.NumRB * 12 / prs.CombSize;
        prsSlotGridRx = reshape(symCellRx, numREPerSym, []);

        symInSlot0 = mod(symIndSlot - 1, numSymPerSlot);
        prsSlotGridRxDestagrd = nrRadar.rx.prsDestaggering(prsSlotGridRx, ...
            prs, symInSlot0);

        % Range IFFT with Windowing
        rangeFft(:, slotKey, nrx) = ...
            sqrt(rangeFFTLen) * ifft(...
            full(prsSlotGridRxDestagrd) .* staticParams.rangeWindow,...
            rangeFFTLen);
    end
end

%% Clutter removal
% Static clutter removal (subtract mean across time/Doppler dim)
rangeFft = rangeFft - mean(rangeFft, 2);

% Truncate range bins
if ~isempty(staticParams.rangeLength)
    rangeFft(staticParams.rangeLength+1:end, :, :) = [];
end

%% Doppler Processing
% Apply Doppler window (broadcast across range and antennas)
rangeFftWindowed = rangeFft .* staticParams.dopplerWindow3D;

% Doppler FFT along dimension 2
rdCube = fftshift(1/sqrt(staticParams.dopplerFftLen) * ...
    fft(rangeFftWindowed, staticParams.dopplerFftLen, 2), 2);

end
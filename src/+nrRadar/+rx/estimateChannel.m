function h = estimateChannel(rxWaveform, staticParams, varargin)
% RUN5GNRAD_RXCHANNELESTIMATION PRS-based frequency-domain channel estimation.
%   H = ESTIMATECHANNEL(RXWAVEFORM, STATICPARAMS)
%   estimates the frequency-domain channel on the PRS resource grid for each
%   receive antenna. The function processes the received OFDM waveform
%   symbol-by-symbol: it selects an FFT window within each OFDM symbol
%   (optionally accounting for a maximum channel length), performs a comb-2
%   reconstruction/combining step (even/odd comb symbols), computes the FFT,
%   and then forms a sparse received grid. The channel estimate is obtained
%   by element-wise multiplication of the sparse received grid with the
%   conjugate of the known transmitted PRS grid.
%
%   H = ESTIMATECHANNEL(..., 'channelLength', L)
%   specifies the assumed channel length L (samples) used to place the FFT
%   window within each OFDM symbol. If omitted or empty, L defaults to the
%   minimum cyclic prefix length in STATICPARAMS.CPLENGTHS.
%
%   INPUTS
%     RXWAVEFORM     : NSAMP×NRX complex array of received time-domain OFDM
%                     samples, where NRX = STATICPARAMS.NCHANRX.
%     STATICPARAMS   : Struct of precomputed/static parameters
%
%   NAME-VALUE PAIRS
%     'channelLength': Scalar integer or [] specifying assumed channel
%                      impulse response length in samples for FFT window
%                      placement. Default: min(cpLengths).
%
%   OUTPUTS
%     H     : 1×NRX cell array. Each cell contains an
%                     ofdmFftLen×nSymTot complex matrix with the PRS-based
%                     frequency-domain channel estimate for that RX antenna.
%
%
%
%   2026 NIST/CTL Steve Blandino
%
%   This file is available under the terms of the NIST License.

    %% Input Parsing
    p = inputParser;
    addRequired(p, 'rxWaveform', @isnumeric);
    addRequired(p, 'staticParams', @isstruct);
    addParameter(p, 'channelLength', [], @(x) isempty(x) || ...
        (isnumeric(x) && isscalar(x)));
    
    parse(p, rxWaveform, staticParams, varargin{:});
    
    % Extract parsed results
    channelLen = p.Results.channelLength;

    %% Unpack Parameters
    prs             = staticParams.prs;
    ofdmFftLen      = staticParams.ofdmFftLen;
    cpLengths       = staticParams.cpLengths;
    symbolIndices   = staticParams.symbolIndices;
    ofdmGrid        = staticParams.ofdmGrid;
    nChanRx         = staticParams.nChanRx;
    carrier         = staticParams.carrier;
    numSlots        = staticParams.numSlots;
    isEvenCombSym   = staticParams.isEvenCombSym; 

    % Set default channel length if not provided by user
    if isempty(channelLen)
        channelLen = min(cpLengths); 
    end

    %% Initialize Buffers
    totalSyms = prs.NumPRSSymbols * staticParams.numberSensingSymbols;
    
    freqDomainSymbolRxStore = zeros(totalSyms * ofdmFftLen, 1);
    linearIdxStore          = zeros(totalSyms * ofdmFftLen, nChanRx);
    
    rows = (1:ofdmFftLen).';
    samplePointer = 0;

    %% Symbol Loop
    for symIdx = 0 : totalSyms - 1

        currentSymGlobal = symbolIndices(symIdx+1);
        cpLen = cpLengths(mod(currentSymGlobal - 1, length(cpLengths)) + 1);
        symbolLen = ofdmFftLen + cpLen;
        
        % Calculate FFT window placement based on channelLen
        symbolStart = min(channelLen, ofdmFftLen/2 + cpLen) + 1;
        symbolStart = max(symbolStart, cpLen + 1);
        tau         = symbolStart - (cpLen + 1);

        for nrx = 1:nChanRx
            % Extract samples
            thisSymbolRx = rxWaveform(samplePointer + 1 : samplePointer + symbolLen, nrx);
            thisSymbolRxNoCP = thisSymbolRx(symbolStart : symbolStart + ofdmFftLen/2 - 1);

            % Comb Combining (SNR Optimization)
            if isEvenCombSym(symIdx+1)
                % Even Comb Logic
                thisSymbolRxNoCP(ofdmFftLen/2-(symbolStart-1-cpLen)+1:end) = ...
                    -thisSymbolRxNoCP(ofdmFftLen/2-(symbolStart-1-cpLen)+1:end);
                thisSymbolRxNoCP(1:ofdmFftLen/2-(symbolStart-1-cpLen)) = ...
                    thisSymbolRxNoCP(1:ofdmFftLen/2-(symbolStart-1-cpLen))/2 - ...
                    thisSymbolRx(symbolStart+ofdmFftLen/2:end)/2;
                
                thisSymbolRxNoCP = circshift(thisSymbolRxNoCP, tau);
                thisSymbolRxReconstructed = [thisSymbolRxNoCP; -thisSymbolRxNoCP];    
            else
                % Odd Comb Logic
                thisSymbolRxNoCP(1:ofdmFftLen/2-(symbolStart-1-cpLen)) = ...
                    thisSymbolRxNoCP(1:ofdmFftLen/2-(symbolStart-1-cpLen))/2 + ...
                    thisSymbolRx(symbolStart+ofdmFftLen/2:end)/2;
                
                thisSymbolRxNoCP = circshift(thisSymbolRxNoCP, tau);
                thisSymbolRxReconstructed = [thisSymbolRxNoCP; thisSymbolRxNoCP];    
            end

            % FFT
            freqDomainSymbolRx = (1 / sqrt(ofdmFftLen)) * fft(thisSymbolRxReconstructed, ofdmFftLen);
            
            % Store indices
            linearIdx = rows + (currentSymGlobal - 1) * ofdmFftLen;
            storeIdxStart = symIdx * ofdmFftLen + 1;
            storeIdxEnd   = (symIdx+1) * ofdmFftLen;
            
            freqDomainSymbolRxStore(storeIdxStart:storeIdxEnd, nrx) = freqDomainSymbolRx;
            linearIdxStore(storeIdxStart:storeIdxEnd, nrx)          = linearIdx;
        end
        samplePointer = samplePointer + symbolLen;
    end

    %% Channel Estimation (Sparse Grid Construction)
    nSymTot  = carrier.SymbolsPerSlot * numSlots;
    N        = size(linearIdxStore, 1);

    lin = linearIdxStore(:);
    row = mod(lin-1, ofdmFftLen) + 1;
    col = floor((lin-1)/ofdmFftLen) + 1;
    
    rxOffset = reshape(repelem((0:nChanRx-1).', N), [], 1);
    colBig   = col + rxOffset * nSymTot;
    val      = freqDomainSymbolRxStore(:);

    rxBig = sparse(row, colBig, val, ofdmFftLen, nSymTot*nChanRx);
    txBig = repmat(conj(ofdmGrid), 1, nChanRx);
    gBig  = rxBig .* txBig;

    h = mat2cell(gBig, ofdmFftLen, repmat(nSymTot, 1, nChanRx));
end
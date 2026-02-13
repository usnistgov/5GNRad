function [rxWaveform] = getRxWaveform(txWaveform, H, staticParams, varargin)
% GETRXWAVEFORM Generate receive waveform 
%   RXWAVEFORM = GETRXWAVEFORM(TXWAVEFORM,H,STATICPARAMS) generates the
%   received time-domain waveform RXWAVEFORM by applying a (possibly
%   time-varying) impulse response H to the transmitted waveform TXWAVEFORM.
%   The function processes the waveform symbol-by-symbol using OFDM
%   parameters in STATICPARAMS, accumulates inter-symbol interference (ISI)
%   naturally via overlap-add, and adds complex AWGN based on the configured
%   noise variance.
%
%   Inputs:
%     TXWAVEFORM     - Transmitted time-domain waveform as a column vector
%                      of complex samples.
%     H              - Channel impulse response with size [L x NSYM x NRX],
%                      where L is the channel length (taps), NSYM is the
%                      number of processed OFDM symbols, and NRX is the
%                      number of receive channels (ports/antennas).
%     STATICPARAMS   - Structure containing OFDM and sensing parameters:
%                      * ofdmFftLen            : FFT length per OFDM symbol
%                      * cpLengths             : cyclic prefix lengths
%                      * symbolIndices         : indices of OFDM symbols in a slot/frame
%                      * nChanRx               : number of receive channels (NRX)
%                      * prs                   : structure with field NumPRSSymbols
%                      * numberSensingSymbols  : number of sensing symbol groups
%                      * snrvar                : SNR-related scale (linear), used
%                                                to derive noise power
%
%   Name-Value Pairs:
%     (none currently used) This input is reserved for future options.
%
%   Output:
%     RXWAVEFORM     - Received time-domain waveform of size
%                      [(length(TXWAVEFORM)+L-1) x NRX]. The waveform includes
%                      the convolution output and added complex Gaussian noise.
%
%   2026 NIST/CTL Steve Blandino
%
%   This file is available under the terms of the NIST License.


 
    %% Unpack Parameters
    ofdmFftLen    = staticParams.ofdmFftLen;
    cpLengths     = staticParams.cpLengths;
    symbolIndices = staticParams.symbolIndices;
    nChanRx       = staticParams.nChanRx;
    prs           = staticParams.prs;
    numSensSym    = staticParams.numberSensingSymbols;
    snrvar        = staticParams.snrvar;

    %% Channel Convolution
    numSamplesTotal = length(txWaveform);
    channelLen      = size(H, 1);
    
    waveformLen = numSamplesTotal + channelLen - 1;
    rxWaveform  = zeros(waveformLen, nChanRx);

    samplePointer = 0;
    totalSyms     = prs.NumPRSSymbols * numSensSym;

    for symIdx = 0 : totalSyms - 1

        % Determine CP length for this specific symbol to stride correctly
        currentSymGlobal = symbolIndices(symIdx+1);
        cpLen = cpLengths(mod(currentSymGlobal - 1, length(cpLengths)) + 1);
        symbolLen = ofdmFftLen + cpLen;

        % Extract the time-domain samples for this specific OFDM symbol
        % Note: Ensure indices don't exceed txWaveform length (safety check)
        idxEnd = min(samplePointer + symbolLen, numSamplesTotal);
        thisSymbolTx = txWaveform(samplePointer + 1 : idxEnd);
        
        if isempty(thisSymbolTx)
            break; 
        end

        for nrx = 1:nChanRx
            % Retrieve impulse response for this symbol/antenna
            hSym = H(:, symIdx+1, nrx);

            % Convolution: Symbol * Channel
            thisSymbolRx = conv(thisSymbolTx, hSym);

            % Superimpose onto the receive buffer (Overlapping adds up)
            % This correctly handles ISI between adjacent symbols
            outStart = samplePointer + 1;
            outEnd   = outStart + length(thisSymbolRx) - 1;
            
            rxWaveform(outStart:outEnd, nrx) = ...
                rxWaveform(outStart:outEnd, nrx) + thisSymbolRx;
        end
        
        % Advance pointer
        samplePointer = samplePointer + symbolLen;
    end

    %% Add AWGN
    noisePower = 1/snrvar;
    noiseScale = sqrt(noisePower/2); 
    
    rxWaveform = rxWaveform+ noiseScale * ...
        (randn(waveformLen, nChanRx) + 1j * randn(waveformLen, nChanRx));

end
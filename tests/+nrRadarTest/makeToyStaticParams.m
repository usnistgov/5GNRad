function staticParams = makeToyStaticParams(varargin)
%MAKETOYSTATICPARAMS Minimal staticParams struct for unit-testing getRxWaveform.
%
% This helper isolates nrRadar.rx.getRxWaveform from the full NR stack.
%
% Name/value:
%   'ofdmFftLen' : FFT length (default 64)
%   'cpLen'      : CP length (default 16)
%   'nChanRx'    : number of receive ports (default 2)
%   'numSensSym' : number of sensing slots (default 4)
%   'numPRSSym'  : number of PRS OFDM symbols per slot (default 1)
%
% Output:
%   staticParams contains only the fields consumed by getRxWaveform.

p = inputParser;
addParameter(p,'ofdmFftLen',64);
addParameter(p,'cpLen',16);
addParameter(p,'nChanRx',2);
addParameter(p,'numSensSym',4);
addParameter(p,'numPRSSym',1);
parse(p,varargin{:});
o = p.Results;

staticParams = struct();
staticParams.ofdmFftLen = o.ofdmFftLen;
staticParams.cpLengths  = repmat(o.cpLen, 14, 1);
staticParams.symbolIndices = (1:(o.numPRSSym*o.numSensSym)).';
staticParams.nChanRx = o.nChanRx;
staticParams.numberSensingSymbols = o.numSensSym;
staticParams.snrvar = Inf; % noise disabled

prs = struct();
prs.NumPRSSymbols = o.numPRSSym;
staticParams.prs = prs;
end

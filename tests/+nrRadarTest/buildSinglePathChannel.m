function [H, channelLen, truth] = buildSinglePathChannel(staticParams, simConfig, receiveArray, varargin)
%BUILDSINGLEPATHCHANNEL Synthetic time-varying single-tap channel for tests.
%
%   [H, channelLen, truth] = nrRadarTest.buildSinglePathChannel(...)
%   returns H with size [channelLen x totalSyms x nChanRx]. The channel has:
%     * a single delay tap (in samples)
%     * a single AoA (az/el in deg)
%     * a single Doppler bin (integer bin index of FFT, post-fftshift)
%
% Name/value:
%   'tapDelaySamples' : integer >=0
%   'dopplerBin'      : integer in [-Nd/2, Nd/2-1]
%   'azDeg','elDeg'   : degrees
%   'amplitude'       : linear amplitude
%
% truth outputs expected bins + physical values.

p = inputParser;
addParameter(p,'tapDelaySamples', 6, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'dopplerBin', 2, @(x)isnumeric(x)&&isscalar(x));
addParameter(p,'azDeg', 10, @(x)isnumeric(x)&&isscalar(x));
addParameter(p,'elDeg', 0, @(x)isnumeric(x)&&isscalar(x));
addParameter(p,'amplitude', 1, @(x)isnumeric(x)&&isscalar(x));
parse(p,varargin{:});
opt = p.Results;

fcHz = simConfig.systemFc;
lambda = 299792458/fcHz;

nChanRx = staticParams.nChanRx;
prs = staticParams.prs;
Nsense = staticParams.numberSensingSymbols;
Ndopp = staticParams.dopplerFftLen;

totalSyms = prs.NumPRSSymbols * Nsense;

% Doppler bin convention: after fftshift, zero is at center (idx = Ndopp/2+1).
% Here we accept dopplerBin in signed bin units (e.g., +2 means 2 bins above 0).
assert(mod(Ndopp,2)==0, 'Test helper expects even dopplerFftLen.');
dopplerBinSigned = opt.dopplerBin;
dopplerBinSigned = max(min(dopplerBinSigned, Ndopp/2-1), -Ndopp/2);

% Generate steering response across RF chains
rxSV = phased.SteeringVector('SensorArray', receiveArray, 'IncludeElementResponse', true);
if isa(receiveArray,'phased.URA')
    a = rxSV(fcHz, [opt.azDeg; opt.elDeg]);
else
    % For ReplicatedSubarray, pass custom subarray weights.
    % NOTE: for our receiveArray construction in tests, the subarray is a ULA
    % and wElem is [numElemPerSub x numSubarrays].
    numElemPerSub = simConfig.rxAntenna.meta.array.M / simConfig.rxAntenna.meta.array.Mprime;
    wElem = reshape(simConfig.rxAntenna.beamformer.wElem, numElemPerSub, []);
    a = rxSV(fcHz, [opt.azDeg; opt.elDeg], wElem);
end

% Ensure expected length
if numel(a) ~= nChanRx
    error('nrRadarTest:ChannelDimMismatch', ...
        'Steering vector length (%d) does not match staticParams.nChanRx (%d).', numel(a), nChanRx);
end

a = a(:).'; % 1 x nChanRx

% Allocate channel
channelLen = max(opt.tapDelaySamples+1, 1);
H = complex(zeros(channelLen, totalSyms, nChanRx, 'double'));

% Build a Doppler phase progression across sensing slots.
% Index over sensing slot (slow time) rather than PRS symbol within the slot.
% k = 0..Nsense-1
kSlow = 0:Nsense-1;
phaseSlow = exp(1j*2*pi*(dopplerBinSigned/ Ndopp) * kSlow);  % length Nsense

% Expand across PRS symbols within each sensing slot
for sym = 1:totalSyms
    k = floor((sym-1) / prs.NumPRSSymbols); % 0..Nsense-1
    H(opt.tapDelaySamples+1, sym, :) = opt.amplitude * phaseSlow(k+1) * a;
end

% Provide ground truth in terms of expected shifted-bin index
truth = struct();
truth.tapDelaySamples = opt.tapDelaySamples;
truth.rangeBin = opt.tapDelaySamples + 1;      % 1-based within IFFT output (approx)
truth.dopplerBinSigned = dopplerBinSigned;
truth.dopplerIndexShifted = (Ndopp/2 + 1) + dopplerBinSigned; % 1-based in fftshifted axis
truth.lambda = lambda;
truth.velocityExpected = dopplerBinSigned * staticParams.prsVelocityResolution;
truth.azDeg = opt.azDeg;
truth.elDeg = opt.elDeg;
end

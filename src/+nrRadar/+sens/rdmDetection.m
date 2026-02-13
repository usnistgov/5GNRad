function [stPositionEstimate_k, velParallel_k, peakRda,peakPostProcSnrEstimate] = ...
    rdmDetection(rdmCube, sensConfig, simConfig, receiveArray, syncOffset, staticParams,varargin) 

p = inputParser;
addParameter(p, 'angleEstimate', []);
parse(p, varargin{:});

angleEstimate      = p.Results.angleEstimate;
if ~isempty(angleEstimate)
    if any(isnan(angleEstimate))
        isAnalogBeamforming = false;
    else
        isAnalogBeamforming = true;
    end
else
    isAnalogBeamforming = false;
end

%% Detection
stPositionEstimate_k = double.empty(0,3);
velParallel_k = double.empty(0,1);
peakRda = double.empty(0,1);
doaEstimationMethod = sensConfig.doaEstimationMethod; %{'barlett', 'FFT', 'music'}


rdmPwrCube = abs(rdmCube).^2;
% figure, surf(20*log10(abs(rdmCube)))
[rangeLength,dopplerLength,numDigitalChain] = size(rdmPwrCube);
elevationLength = sensConfig.elFftLen;
azimuthLength = sensConfig.azFftLen;

%% 2-D projection
rdmPwr   = proj2D(rdmPwrCube, [1 2]);
% rdMapFloor = median(rdmPwr(:));
% rdMapThreshold = rdMapFloor*10^(sensConfig.rdaThreshold/10);
% rdmPwr(rdmPwr<rdMapThreshold) = 0;

%% 2-D CFAR
grdRD = [sensConfig.cfarGrdCellRange,    sensConfig.cfarGrdCellVelocity];
trnRD = [sensConfig.cfarTrnCellRange,    sensConfig.cfarTrnCellVelocity];
thrCFAR = sensConfig.cfarThreshold;
[rdmPwrMask, rdmPwrNoise]   = nrRadar.sens.cfar2D(rdmPwr,   grdRD,   trnRD,   thrCFAR);

%% Find Peaks
peaks = nrRadar.sens.pick_peaks_nms(rdmPwrMask, rdmPwrNoise,  thrCFAR, grdRD);

%% Process peaks
if isempty(peaks)
    return;
end



% Clamp + columnize indices (also consider rounding if peaks can be non-integer)
r = peaks(:,1);
d = peaks(:,2);
r = max(1, min(rangeLength, r(:)));
d = max(1, min(dopplerLength, d(:)));
nDet = numel(r);

linRD   = sub2ind([rangeLength dopplerLength], r, d);
rdmFlat = reshape(rdmCube, rangeLength*dopplerLength, numDigitalChain);     % view
rdmFlatDetected       = rdmFlat(linRD,:).';              % [nAnt x nDet]
peakPostProcSnrEstimate = 10*log10(rdmPwrMask(linRD)./(rdmPwrNoise(linRD)+eps));


if isAnalogBeamforming
    azIdx = angleEstimate(1)*ones(nDet,1);
    elIdx = angleEstimate(2)*ones(nDet,1);
    azHat = azIdx;
    elHat = elIdx;
else

    switch doaEstimationMethod
        case 'beamspaceFFT'
            rdmCubePreFft  = reshape(rdmFlatDetected, staticParams.M, staticParams.N, nDet) .* staticParams.W;
            angleMapPerDetection = fft2(rdmCubePreFft, staticParams.Nfft_v, staticParams.Nfft_h);
            angleMapPerDetection = fftshift(angleMapPerDetection,1);
            angleMapPerDetection = fftshift(angleMapPerDetection,2);
            Pwr  = abs(angleMapPerDetection).^2;                    % [Nfft_v x Nfft_h x nDet]
            [~, idx] = max(reshape(Pwr, [], nDet), [], 1);
            [elIdx, azIdx] = ind2sub([staticParams.Nfft_v, staticParams.Nfft_h], idx);
        case 'barlettScan'

            az_search = -60:1:60;
            el_search = 0:1:60;
            wElem = simConfig.rxAntenna.beamformer.wElem;
            respMatrix = nrRadar.sens.buildRespMatrixBartlett(receiveArray, staticParams.fc, az_search, el_search, wElem);
            % Dictionary A: [NRF x K], K=nAz*nEl
            [NRF, nAz, ~] = size(respMatrix);
            A = reshape(respMatrix, NRF, []);         % NRF x K

            % Bartlett scan for all angles and all detections:
            Y = A' * rdmFlatDetected;                                % K x nDet   (conj transpose = matched filter)
            S = abs(Y).^2;                             % K x nDet

            % Pick best angle per detection
            [bestVal, idx] = max(S, [], 1);            % 1 x nDet

            azIdx = mod(idx-1, nAz) + 1;
            elIdx = floor((idx-1)/nAz) + 1;

            azHat = az_search(azIdx).';                % nDet x 1
            elHat = el_search(elIdx).';                % nDet x 1
            bestVal = bestVal.';

        case 'MUSIC'
            error('MUSIC not implemented yet')
    end
end

pk = [peaks(:,1:2), elIdx(:), azIdx(:), peaks(:,end)];

dims = struct('R', rangeLength, 'D', dopplerLength, 'El', elevationLength, 'Az', azimuthLength, ...
    'wrapD', true, 'wrapAz', true);    % wrap Doppler/Azimuth

if isempty(pk)
    detection4D = [];
else
    opts = struct();
    opts.minVal    = pk(1,end)/sensConfig.dbscanMinMaxRatio;        % discard weak peaks (linear units)
    opts.norm      = [2 2 1 1];% normalize bin distances per axis
    opts.valWeight = 0.5;      % include value similarity (0 = ignore)
    opts.valScale  = [];       % leave empty to auto (robust MAD)
    opts.eps       = 6;      % clustering radius in normalized space
    opts.minPts    = 1;        % at least 1 peak to form a cluster

    cl = nrRadar.sens.cluster_peaks_4d(pk, dims, opts);

    detection4D = table2array(cl.clusters(:,[8:11 6]));
    detection4D = nrRadar.sens.suppressSidelobes(detection4D, inf, [rangeLength dopplerLength elevationLength azimuthLength]);
    detection4D(:,1:4) = round(detection4D(:,1:4));
    % detection4D(:,end) = [];
end

%% Parameter estimation

if isempty(detection4D)
    return
end

    nDetectedTarget = size(detection4D, 1);

    azEstimate_k  = zeros(nDetectedTarget, 1);
    elEstimate_k  = zeros(nDetectedTarget, 1);
    peakRda       = zeros(nDetectedTarget, 1);

    for det = 1:nDetectedTarget
        if isAnalogBeamforming
            azEstimate_k(det) = azHat(det)+ simConfig.trpYawDeg;
            elEstimate_k(det) = elHat(det);
            peakRda(det) = 10*log10(detection4D(det, 5));
        else
            switch doaEstimationMethod
                case 'barlettScan'
                    azEstimate_k(det) = azHat(det)+ simConfig.trpYawDeg;
                    elEstimate_k(det) = elHat(det);
                    peakRda(det) = bestVal(det);
                case 'beamspaceFFT'
                    azEstimate_k(det) = staticParams.azGrid(detection4D(det,3), detection4D(det,4)) + simConfig.trpYawDeg;
                    elEstimate_k(det) = nrRadar.array.converAliasAngleToReal(staticParams.elGrid(detection4D(det,3), detection4D(det,4)), round(receiveArray.ElementSpacing(2)/staticParams.lambda*100)/100);
                    peakRda(det)      = 10*log10(detection4D(det, 5));
            end
        end
    end

    rngEstimate_k = staticParams.rangeBinDestgrd(detection4D(:,1)) + syncOffset;
    velParallel_k = staticParams.velocityBin(detection4D(:,2));


    stPositionEstimate_k = nrRadar.util.estimateScatteringGeometry(staticParams.txPos, staticParams.rxPos, ...
        2*rngEstimate_k(:)/staticParams.c*1e9, [azEstimate_k, elEstimate_k], ...
        'keepSectorAz', [-60 60] + simConfig.trpYawDeg, ...
        'txHeading', [simConfig.trpYawDeg 0], ...
        'rxHeading', [simConfig.trpYawDeg 0]);

    validEstimateIndex = all(~isnan(stPositionEstimate_k), 2) & (stPositionEstimate_k(:,3) >= 0);

    stPositionEstimate_k = stPositionEstimate_k(validEstimateIndex, :);

    if ~isempty(velParallel_k) && ~isscalar(velParallel_k)
        velParallel_k = velParallel_k(validEstimateIndex);
    end

    if ~isempty(peakRda)
        peakRda = peakRda(validEstimateIndex);
    end






end


function M = proj2D(cube, keepDims)
% PROJ2D Project a 4-D cube onto a 2-D plane by summing across other dimensions
%   M = PROJ2D(CUBE, KEEPDIMS) reduces a 4-D array CUBE by summing along
%   all dimensions not specified in KEEPDIMS. The result M is a 2-D matrix
%   representing the projection of CUBE over the selected dimensions.
%
%   Inputs
%   ------
%   CUBE     : [R x D x El x Az] 4-D numeric array (e.g., range-Doppler-angle cube)
%   KEEPDIMS : [1 x 2] vector specifying which dimensions to preserve
%               (e.g., [1 2] for range–Doppler, [1 3] for range–elevation)
%
%   Outputs
%   -------
%   M : 2-D matrix obtained by summing CUBE along all dimensions not in KEEPDIMS.
%       The output is squeezed to remove singleton dimensions.
%
%   Example
%   -------
%       % Example with a random cube
%       cube = rand(32, 32, 8, 8);
%       M = proj2D(cube, [1 2]);   % Range–Doppler projection
%
%   2025 NIST/CTL Steve Blandino
%
%   This file is available under the terms of the NIST License.

allDims = 1:4;
dropDims = setdiff(allDims, keepDims);
M = cube;

for dd = dropDims
    M = sum(M, dd);
end

M = squeeze(M);
end
function respMatrix = buildRespMatrixBartlett(receiveArray, fcHz, az_search, el_search, wElem, varargin)
%BUILDRESPMATRIXBARTLETT Build manifold (respMatrix) for Bartlett scan.
%
% respMatrix: [NRF x nAz x nEl] where NRF = number of RF-chain ports.

p = inputParser;
addParameter(p, 'IncludeElementResponse', true);
addParameter(p, 'NormalizeColumns', true);   % normalize each a(az,el)
parse(p, varargin{:});
incElem = p.Results.IncludeElementResponse;
doNorm  = p.Results.NormalizeColumns;

rxSV = phased.SteeringVector('SensorArray', receiveArray, ...
    'IncludeElementResponse', incElem);

nAz = numel(az_search);
nEl = numel(el_search);

% Angle grid with az varying fastest (matches idx mapping: idx=(el-1)*nAz+az)
[AZ, EL] = ndgrid(az_search, el_search);
ang = [AZ(:).'; EL(:).'];  % 2 x (nAz*nEl)

% Custom subarray weights (for phased.ReplicatedSubarray with 'Custom' steering)
% wElem must be [numElemPerSub x numSubarrays]
if isa(receiveArray,'phased.URA')
    A = rxSV(fcHz, ang);
else
    A = rxSV(fcHz, ang, wElem);     % [NRF x (nAz*nEl)]
end

if doNorm
    A = A ./ max(vecnorm(A,2,1), eps(class(A)));
end

respMatrix = reshape(A, size(A,1), nAz, nEl); % [NRF x nAz x nEl]
end

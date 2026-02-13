function [receiveArray, wElem] = buildReceiveArray(simConfig)
%BUILDRECEIVEARRAY Construct a Phased Array object consistent with simConfig.rxAntenna.
%
% Returns:
%   receiveArray - phased.URA or phased.ReplicatedSubarray
%   wElem        - custom subarray weights matrix for ReplicatedSubarray ([]) for URA

fcHz = simConfig.systemFc;
lambda = 299792458/fcHz;

rxMeta = simConfig.rxAntenna.meta;
M = rxMeta.array.M;
N = rxMeta.array.N;
dV = rxMeta.array.dV_lambda * lambda;
dH = rxMeta.array.dH_lambda * lambda;

rxMode = rxMeta.name;

if strcmp(rxMode,'hybrid')
    Mprime = rxMeta.array.Mprime;
    Nprime = rxMeta.array.Nprime;
    numElemPerSub = M / Mprime;

    sULA = phased.ULA('NumElements', numElemPerSub, ...
        'Element', phased.NRAntennaElement, ...
        'ElementSpacing', dV, ...
        'ArrayAxis', 'z');

    receiveArray = phased.ReplicatedSubarray('Subarray', sULA, ...
        'GridSize', [Mprime Nprime], ...
        'SubarraySteering', 'Custom', ...
        'GridSpacing', [numElemPerSub*dV, dH]);

    wElem = reshape(simConfig.rxAntenna.beamformer.wElem, numElemPerSub, []);
else
    receiveArray = phased.URA([M N], 'ElementSpacing', [dV dH]);
    receiveArray.Element = phased.NRAntennaElement;
    wElem = [];
end
end

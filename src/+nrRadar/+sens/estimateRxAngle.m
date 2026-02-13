function angleEstimate = estimateRxAngle(tx, method, fc, txPV, cir, aoaAzLOS, aoaElLOS, varargin)
%ESTIMATERXANGLE Choose beam direction by 'ideal' | 'nearest' | 'scan'
%
% rxPV: [Nrx x P] steering vectors for each MPC path direction (AoA)
% cir : [P x T] complex path gains vs time (before RX digital beamforming)
% scanvector: [Nscan x 2] [az el] in degrees

p = inputParser;
addParameter(p,'applyFrontBackMask',false);
addParameter(p,'boresightAzEl',[0 0]);  % defines "front" direction if masking is used
parse(p,varargin{:});
applyMask = p.Results.applyFrontBackMask;
boresight = p.Results.boresightAzEl;

method = validatestring(method, {'ideal','nearest','scan'});

hpbwH = 0.886*2/tx.Size(2)*180/pi;
hpbwV = 0.886*2/tx.Size(1)*180/pi;
scanStepH = floor(hpbwH);
scanStepV = floor(hpbwV);
azimuthRange = -60:scanStepH:60;
elevationRange = -90:scanStepV:90;
[az, el] = meshgrid(azimuthRange, elevationRange);
scanvector = [az(:), el(:)];


switch method
    case 'ideal'
        angleEstimate = [aoaAzLOS', aoaElLOS'];

    case 'nearest' % Single target
        % nearest direction in the codebook to the LOS direction
        ref = [aoaAzLOS(1), aoaElLOS(1)];
        d = sum(abs(scanvector - ref).^2, 2);
        [~,id] = min(d);
        angleEstimate = scanvector(id,:);

    case 'barlettScan' % Single target
        % Build candidate beamformers (one per scan angle)
        % rxSV(fc, scanvector.') returns [Nrx x Nscan] if scanvector is 2xNscan

        txSV = phased.SteeringVector('SensorArray', tx);


        W = conj(txSV(fc, scanvector.'));    % [Nrx x Nscan]

        % Compute per-beam gains for each path:
        % W' * rxPV -> [Nscan x P]
        G = W.' * txPV;                       % [Nscan x P]

        % Compute received beam output per beam across time:
        % y_b(t) = sum_p G(b,p) * cir(p,t)
        % Y = G * cir -> [Nscan x T]
        Y = G * (cir-mean(cir,2));


        % Power metric
        pow = sum(abs(Y).^2, 2);             % [Nscan x 1]

        % Optional front/back masking using 3D dot product test
        if applyMask
            vScan = angle2vector(scanvector(:,1), 90 - scanvector(:,2), 1);   % [Nscan x 3]
            vBore = angle2vector(boresight(1), 90 - boresight(2), 1);         % [1 x 3]
            isFront = (vScan * vBore(:)) >= 0;  % keep same hemisphere as boresight
            pow(~isFront) = 0;
        end

        [~,id] = max(pow);
        angleEstimate = scanvector(id,:);
end
end

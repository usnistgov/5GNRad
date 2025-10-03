function [PL,isLos] = computePathLoss(scenario, fc, txPos, rxPos, isBackground)
% COMPUTEPATHLOSS Compute path loss and LOS condition for 5GNR channel models
%   [PL, ISLOS] = COMPUTEPATHLOSS(SCENARIO, FC, TXPOS, RXPOS, ISBACKGROUND)
%   computes the path loss (PL) in dB and determines whether the link is
%   line-of-sight (ISLOS), based on 3GPP TR 38.901 models for RMaAV, UMaAV,
%   and UMiAV scenarios.
%
%   Inputs:
%     scenario     - Character string: 'RMaAV', 'UMaAV', or 'UMiAV'
%     fc           - Carrier frequency in Hz
%     txPos        - 3D transmitter position [x, y, z] in meters
%     rxPos        - 3D receiver (user terminal) position [x, y, z] in meters
%     isBackground - Boolean flag to force NLOS for background channel modeling
%
%   Outputs:
%     PL     - Path loss in dB, including shadow fading
%     isLos  - Boolean indicating if the link is line-of-sight (1 = LOS, 0 = NLOS)
%
%   This function:
%     - Computes LOS probability based on scenario and geometry
%     - Computes deterministic path loss formulas from 3GPP models
%     - Adds scenario-specific shadow fading using computeShadowFading()
%
%
%   See also: COMPUTELOSPROBABILITY, COMPUTEUTEffectiveHeight, COMPUTESHADOWFADING

%   2025 NIST/CTL Steve Blandino
%   This file is available under the terms of the NIST License.

% Constants
c = 3e8; % Speed of light in m/s
fcGHz = fc / 1e9; % Convert frequency to GHz

% Compute distances
d3D = norm(txPos - rxPos); % 3D distance
d2D = norm(txPos(1:2) - rxPos(1:2)); % 2D distance in xy-plane
hUT = rxPos(3); % RX height (z component)
hBS = txPos(3);

losProbability = computeLOSProbability(scenario, d2D, hUT);
if isBackground
    isLos = 0;
else
isLos =rand<losProbability;
end
isAfterBP = [];

% Check valid height and range for each model
if hUT < 1.5 || hUT > 300
if hUT <1.5
hUT = 1.5;
else 
hUT = 300;
end
    warning('hUT out of range:. Value has been set to %.1f m.', hUT);
end

switch scenario
    case 'RMaAV'
        % if d2D<10 || d2D > 10e3
        %     error('d2D out of range for RMa-AV: 10m d2D  10km');
        % end

        if hUT>= 1.5 && hUT<=10
            h = 5; %avg. building height

            dBP = 2*pi*hBS *hUT*fc/c ;

            PL1 = @(d) 20 * log10((40 * pi * d * fcGHz) / 3) ...
                + min(0.03 * hUT^1.72, 10) * log10(d) ...
                - min(0.044 * hUT^1.72, 14.77) ...
                + 0.002 * log10(h) * d;
            isAfterBP = d2D >= dBP;
            if isAfterBP
                PL = PL1(dBP) + 40 * log10(d3D / dBP);
            else
                PL = PL1(d3D);
            end

            if ~isLos
                W = 10;
                PL_RMa_NLOS = 161.04 - 7.1 * log10(W) + 7.5 * log10(h) ...
                    - (24.37 - 3.7 * (h / hBS)^2) * log10(hBS) ...
                    + (43.42 - 3.11 * log10(hBS)) * (log10(d3D) - 3) ...
                    + 20 * log10(fc) - (3.2 * (log10(11.75 * hUT))^2 - 4.97);
                PL = max(PL, PL_RMa_NLOS);
            end


        elseif hUT> 10 && hUT<=300
            PL = max(23.9 - 1.8 * log10(hUT),20)...
                *log10(d3D)+ 20 * log10((40 * pi * fcGHz) / 3);
            if ~isLos
                PL_RMa_NLOS =   - 12 + (35 - 5.3*log10(hUT))*log10(d3D) + 20*log10((4*pi*fcGHz)/3);
                PL = max(PL, PL_RMa_NLOS);
            end

        else
            error('hUT out of range for RMa-AV: 1.5  hUT  300m');

        end

    case 'UMaAV'
        if hUT>= 1.5 && hUT<=22.5
            hBS1 = hBS-hUT;
            hE = computeUTEffectiveHeight(d2D, hUT,scenario);
            hUT1 = hUT-hE;
            dBP1 = 4 * (hBS1) * (hUT1) * fc / c;

            if d2D >= 10 &&  d2D<=dBP1
                PL = 28.0 + 22 * log10(d3D) + 20 * log10(fcGHz);
            elseif d2D > dBP1 &&  d2D<=5e3
                PL = 28.0 + 40 * log10(d3D) + 20 * log10(fcGHz) - 9 *log10(dBP1^2+(hBS-hUT)^2);
            else
                error('')
            end
            if ~isLos
                PL_UMa_NLOS = 13.54+39.08*log10(d3D) + 20*log10(fcGHz) - 0.6*(hUT-1.5);
                PL = max(PL_UMa_NLOS, PL);
            end


        elseif hUT>22.5 && hUT<300
            if isLos
                if d2D > 4e3
                    error('d2D out of range for UMa-AV: d2D  4km');
                end
                PL = 28.0 + 22 * log10(d3D) + 20 * log10(fcGHz);
            else
                PL = -17.5+(46-7*log10(hUT))*log10(d3D) + 20*log10((40*pi*fcGHz)/3);
            end

        else
            error('')
        end

    case 'UMiAV'
        if d2D > 4e3
            error('d2D out of range for UMi-AV: d2D  4km');
        end

        if hUT>= 1.5 && hUT<=22.5
            hBS1 = hBS-hUT;
            hE = computeUTEffectiveHeight(d2D, hUT,scenario);
            hUT1 = hUT-hE;
            dBP1 = 4 * (hBS1) * (hUT1) * fc / c;

            % Choose PL based on d2D
            if d2D>=10 && d2D <= dBP1
                PL = 32.4 + 21 * log10(d3D) + 20 * log10(fcGHz);
            elseif d2D>dBP1 && d2D <= 5e3
                PL = 32.4 + 40 * log10(d3D) + 20 * log10(fcGHz) ...
                    - 9.5 * log10((dBP1)^2 + (hBS - hUT)^2);
            else
                error('')
            end
            if ~isLos
                PL_UMi_NLOS = 35.3*log10(d3D)+22.4 + 21.3*log10(fcGHz) - 0.3*(hUT-1.5);
                PL = max(PL_UMi_NLOS, PL);
            end

        elseif hUT> 22.5 && hUT<300
            PL1 = 92.4 + 20*log10(fcGHz) + 20*log10(d3D*1e-3);
            PL = max(PL1, ...
                30.9 + (22.25 - 0.5 * log10(hUT)) * log10(d3D) + 20 * log10(fcGHz));

            if ~isLos
                PL_UMi_NLOS = 32.4 + (43.2-7.6*log10(hUT))+20*log10(fc);
                PL = max(PL, PL_UMi_NLOS);
            end
        else
            error('')
        end

    otherwise
        error('Invalid scenario. Choose "RMa-AV", "UMa-AV", or "UMi-AV".');
end
sigma_SF = computeShadowFading(scenario, isLos, hUT, isAfterBP);
PL = PL + randn*sigma_SF;

end

function losProbability = computeLOSProbability(scenario, d2D, hUT)
% COMPUTELOSPROBABILITY Compute LOS probability based on 3GPP 38.901/36.873
%   LOSPROB = COMPUTELOSPROBABILITY(SCENARIO, D2D, HUT) returns the line-of-sight (LOS)
%   probability for a given scenario, horizontal distance, and user height.
%
%   Inputs:
%     scenario - Scenario string: 'RMaAV', 'UMaAV', or 'UMiAV'
%     d2D      - 2D ground distance between TX and RX [meters]
%     hUT      - User terminal (UE) height [meters], must be within 1.5–300
%
%   Output:
%     losProbability - Line-of-sight probability ∈ [0, 1]
%
%
%   Example:
%     p = computeLOSProbability('UMiAV', 50, 1.5);
%
%   See also: 3GPP TR 38.901, 3GPP TR 36.873

%   2025 NIST/CTL Steve Blandino
%   This file is available under the terms of the NIST License.



if hUT < 1.5
    warning('hUT should be above 1.5m. Setting hUT = 1.5')
    hUT = 1.5;
end
if hUT > 300
    warning('hUT should be below 300m. Setting hUT = 300')
    hUT = 300;
end
        
switch scenario
    case 'RMaAV'
        if hUT >= 1.5 && hUT <= 10
            if d2D <=10
                losProbability = 1;
            else
                losProbability = exp(-((d2D-10)/100));
            end
        elseif hUT > 10 && hUT <= 40
            % Compute parameters p1 and d1
            p1 = max(15021 * log10(hUT) - 16053, 1000);
            d1 = max(1350.8 * log10(hUT) - 1602, 18);  %m

            % Compute P_LOS based on the given conditions
            if d2D <= d1
                losProbability = 1;
            else
                losProbability = (d1 / d2D)+ exp(-d2D / p1) * (1 - d1/d2D);
            end
        elseif hUT> 40 && hUT <= 300
            losProbability = 1;
        else
            error('hUT should be between 1.5m and 300')
        end
    case 'UMaAV'
        if hUT >= 1.5 && hUT <= 22.5
            % 38.901
            if d2D <=18
                losProbability = 1;
            else
                if hUT <= 13
                    C_prime = 0;
                elseif hUT > 13 && hUT <= 23
                    C_prime = ((hUT - 13) / 10)^1.5;
                else
                    error('hUT must be in the range hUT  23m');
                end
                d2Dout = d2D; % assume target is outdoor
                losProbability = ((18 / d2Dout) +exp(-(d2Dout/63))*(1 - 18/d2Dout))* ...
                    (1 + C_prime * (5 / 4) * (d2Dout / 100)^3 * exp(-d2Dout / 150));

            end


        elseif hUT > 22.5 && hUT <= 100
            % Compute parameters p1 and d1
            p1 = 4300*log10(hUT)-3800;
            d1 = max(460*log10(hUT)-700,18);  %m

            % Compute P_LOS based on the given conditions
            if d2D <= d1
                losProbability = 1;
            else
                losProbability = (d1 / d2D)+ exp(-d2D / p1) * (1 - d1/d2D);
            end

        elseif hUT > 100 && hUT <= 300
            losProbability = 1;

        else
            error('hUT should be between 1.5m and 300')
        end

    case 'UMiAV'
        if hUT >= 1.5 && hUT <= 22.5
            d2Dout = d2D; % assume target is outdoor
            losProbability = ((18 / d2Dout) +exp(-(d2Dout/36))*(1 - 18/d2Dout));
        elseif  hUT > 22.5 && hUT <= 300
            % Compute parameters p1 and d1
            p1 = 233.98*log10(hUT)-0.95;
            d1 = max(294.05*log10(hUT)-432.94,18);  %m
            losProbability = (d1 / d2D)+ exp(-d2D / p1) * (1 - d1/d2D);
        end

end
end
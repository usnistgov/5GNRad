function hE = computeUTEffectiveHeight(d2D, hUT, scenario)
% COMPUTEUTEFFECTIVEHEIGHT Compute effective environment height (hE)
%   HE = COMPUTEUTEFFECTIVEHEIGHT(D2D, HUT, SCENARIO) returns the effective 
%   environment height hE based on the user terminal height (hUT), 2D ground 
%   distance (d2D), and propagation scenario.
%
%   Inputs:
%     d2D     - 2D horizontal distance between TX and RX [meters]
%     hUT     - User terminal (UE) height [meters], valid range: 1 ≤ hUT ≤ 300
%     scenario - Scenario string: 'UMiAV' or 'UMaAV'
%
%   Output:
%     hE      - Effective environment height [meters], used in breakpoint
%               distance and path loss calculations
%
%   Behavior:
%     - For 'UMiAV', hE is fixed to 1.0 m.
%     - For 'UMaAV', a probabilistic model (from 3GPP TR 36.777) is used:
%         * A function g(d2D) controls the randomness of hE.
%         * With high probability, hE = 1.0; otherwise hE is randomly selected
%           from discrete values between 12 and floor(hUT - 1.5).
%
%   Example:
%     hE = computeUTEffectiveHeight(50, 20, 'UMaAV');
%
%   See also: COMPUTEPATHLOSS, COMPUTELOSPROBABILITY

%   2025 NIST/CTL Steve Blandino
%   This file is available under the terms of the NIST License.

switch scenario
    case 'UMiAV'
         hE = 1.0;
    case 'UMaAV'
        % Check input validity
        if hUT < 1 || hUT > 300
            error('hUT out of range: 1m  hUT  300m');
        end
        if d2D < 0
            error('d2D must be non-negative');
        end

        % Compute g(d2D) as per the formula
        if d2D <= 18
            g_d2D = 0;
        else
            g_d2D = (5/4) * (d2D / 100)^3 * exp(-d2D / 150);
        end

        % Compute C(d2D, hUT)
        if hUT < 13
            C_d2D_hUT = 0;
        elseif hUT <= 23
            C_d2D_hUT = ((hUT - 13) / 10)^1.5 * g_d2D;
        else
            C_d2D_hUT = g_d2D;
        end

        % Compute hE for UMa and UMi
        probability = 1 / (1 + C_d2D_hUT);

        if rand() < probability
            hE = 1.0; % UMi effective environment height
        else
            if floor(hUT - 1.5)<12
                hE = 1.0; % Default to 1.0m if the range is invalid
            else
                hE = randi([12, floor(hUT - 1.5)]); % Random integer selection
            end
        end
end


end
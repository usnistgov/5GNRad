function sigma_SF = computeShadowFading(scenario, isLOS, hUT, varargin)
% COMPUTESHADOWFADING Shadow fading standard deviation in dB
%   SIGMA_SF = COMPUTESHADOWFADING(SCENARIO, ISLOS, HUT) returns the standard
%   deviation of shadow fading based on the scenario, LOS condition, and UE height.
%
%   SIGMA_SF = COMPUTESHADOWFADING(..., ISAFTERBP) includes a flag for whether
%   the UE is beyond the breakpoint distance (RMaAV only, hUT ≤ 10).
%
%   Inputs:
%     scenario   - Scenario string: 'RMaAV', 'UMaAV', or 'UMiAV'
%     isLOS      - Boolean flag indicating line-of-sight condition
%     hUT        - User terminal height in meters (1.5 ≤ hUT ≤ 300)
%     varargin   - Optional:
%                   * isAfterBP (boolean): Only for RMaAV when hUT ≤ 10
%
%   Output:
%     sigma_SF   - Shadow fading standard deviation in dB
%
%   See also: COMPUTEPATHLOSS

%   2025 NIST/CTL Steve Blandino
%   This file is available under the terms of the NIST License.

if hUT < 1.5 || hUT > 300
    error('hUT must be within the range 1.5m < hUT < 300m');
end

% Compute shadow fading std based on scenario and condition
switch scenario
    case 'RMaAV'
        if isLOS
            if hUT <= 10
                isAfterBP = varargin{1};
                if isAfterBP
                    sigma_SF = 6;
                else
                    sigma_SF = 4;
                end
            else
                sigma_SF = 4.2 * exp(-0.0046 * hUT);
            end
        else
            if hUT <= 10
                sigma_SF = 8;%6
            else
                sigma_SF = 6;
            end
        end

    case 'UMaAV'
        if isLOS
            if hUT <= 22.5
                sigma_SF = 4;
            else
                sigma_SF = 4.64 * exp(-0.0066 * hUT);
            end
        else
            if hUT <= 22.5
                sigma_SF = 4;
            else
                sigma_SF = 6;
            end
        end

    case 'UMiAV'
        if isLOS
            if hUT <= 22.5
                sigma_SF = 4;
            else
                sigma_SF = max(5 * exp(-0.01 * hUT), 2);
            end
        else
            if hUT <= 22.5
                sigma_SF = 7.82;
            else
                sigma_SF = 8;
            end
        end

    otherwise
        error('Invalid scenario. Choose RMa-AV, UMa-AV, or UMi-AV');
end
end
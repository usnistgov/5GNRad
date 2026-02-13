function hE = computeEffectiveHeight(d2D, hUT)
% COMPUTEEFFECTIVEHEIGHT Compute effective environment height hE (3GPP TR 36.777)
%   HE = COMPUTEEFFECTIVEHEIGHT(D2D, HUT) returns the effective environment
%   height based on the 3GPP TR 36.777 model for urban macro (UMa) and 
%   urban micro (UMi) scenarios.
%
%   Inputs:
%     d2D - 2D ground distance between transmitter and receiver [m]
%     hUT - User terminal (UE) height [m], valid range: 1 ≤ hUT ≤ 300
%
%   Output:
%     hE  - Effective environment height [m]
%
%   The method models a probabilistic selection between:
%     - UMi: fixed effective height of 1.0 m
%     - UMa: random selection from set {12, 15, ..., floor(hUT - 1.5)} based on probability
%
%   Reference:
%     3GPP TR 36.777, Annex B.1.2 - Effective Environment Height Model
%
%   Example:
%     hE = computeEffectiveHeight(50, 20);
%
%   See also: TR 36.777 Table B.1.2.1-1

%   2025 NIST/CTL Steve Blandino
%   This file is available under the terms of the NIST License.

% Check input validity
if hUT < 1 || hUT > 300
    error('hUT out of range: 1m ≤ hUT ≤ 300m');
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

% UMi: hE = 1.0m
% UMa: hE chosen randomly from {12, 15, ..., (hUT - 1.5)}
if rand() < probability
    hE = 1.0; % UMi effective environment height
else
    if isempty(hE_values)
        hE = 1.0; % Default to 1.0m if the range is invalid
    else
        hE = randi([12, floor(hUT - 1.5)]); % Random integer selection
    end
end
end
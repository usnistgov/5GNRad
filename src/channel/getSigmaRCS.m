function [sigmaRCSdBsm, sigmaMdBsm, sigmaDdB, sigmaSdB] = getSigmaRCS(st, varargin)
%GETSIGMARCS 38.901 Radar Cross Section
% Return RCS components (total, mean, directional, random) in dBsm for a
% given target type (st) and scenario (single/multi scatterers, mono/bistatic).
%
%   This function returns the total RCS in dBsm along with its components:
%     - RCS
%     - Mean component (sigmaM in dBsm)
%     - Directional gain component (sigmaD in dB)
%     - Log-normal fluctuation component (sigmaS in dB)
%
%   The function supports angular-dependent and angular-independent targets
%   based on 3GPP RCS tables (TR 38.901 Section 7.9.2.1-x).
%
%   USAGE:
%     [sigmaRCS, sigmaM, sigmaD, sigmaS] = getSigmaRCS('uav-small');
%     [sigmaRCS, sigmaM, sigmaD, sigmaS] = getSigmaRCS('vehicle', 'spst', 'multi', 'stDirection', [60, 100], 'N', 1000);
%
%   INPUTS:
%     st              - Target type (string):
%                      'uav-small', 'uav-large', 'human', 'agv', or 'vehicle'
%
%   OPTIONAL PARAMETERS (name-value pairs):
%     'spst'          - Scattering point type: 'single' (default) or 'multi'
%     'stDirection'   - [theta, phi] direction of incidence in degrees (used for angular-dependent targets)
%     'anglesRCS'- [theta_i, phi_i, theta_s, phi_s] (used for bistatic sensing, optional)
%     'N'             - Number of RCS samples to generate (only for random models)
%     'returnLargeScale'  - sigmaRCSdBsm is sigmaM.
%
%
%   OUTPUTS:
%     sigmaRCSdBsm    - Total RCS value(s) in dBsm (size 1xN if N is specified, else scalar)
%     sigmaMdBsm      - Mean RCS value (scalar, in dBsm)
%     sigmaDdB        - Directional gain component (scalar, in dB)
%     sigmaSdB        - Random fluctuation component (scalar, in dB)
%
%   NOTES:
%     - Angular-dependent targets include: 'uav-large', 'vehicle', 'agv'
%     - 'uav-small' and 'human' are angular-independent and ignore direction inputs
%     - For log-normal fluctuation (sigmaSdB), the mean in dB is adjusted
%       using the formula: μ = -ln(10)/20 * σ² to ensure mean = 1 in linear scale
%     - If bistatic is detected based on 'bistaticAngles', the function enters bistatic mode
%
%   EXAMPLES:
%     % Generate 1000 RCS samples for a vehicle with single scatterer at [60,90]
%     [rcs, sigmaM, sigmaD, sigmaS] = getSigmaRCS('vehicle', 'spst', 'single', 'stDirection', [60,90], 'N', 1000);
%
%
%     Generate RCS for large scale calibration for UAV STs
%     rcs = getSigmaRCS('uav-small', 'returnLargeScale',1)
%
%     Generate RCS for full calibration for UAV STs
%     rcs = getSigmaRCS('uav-small');
%
%     Generate RCS for large scale calibration for automotive STs
%     rcs = getSigmaRCS('vehicle', 'returnLargeScale',1)
%
%     Plot pattern of car RCS pattern
%     for i = 1:180
%       radiationCar(i) = getSigmaRCS('vehicle', 'spst', 'single', 'stDirection', [i,90], 'N', 1);
%     end
%     figure, plot(radiationCar)
%
%     Reference: 38.901 7.9.2.1 RCS of a sensing target

%   2025 NIST/CTL Steve Blandino

%   This file is available under the terms of the NIST License.

% Define valid options
validST = {'uav-small', 'uav-large', 'human-model-1', 'agv', 'vehicle'};
validSPST = {'single', 'multi'};

% Input parser setup
p = inputParser;
addRequired(p, 'st', @(x) any(validatestring(x, validST)));
addParameter(p, 'spst', 'single', @(x) any(validatestring(x, validSPST)));
addParameter(p, 'stDirection', [], @(x) isnumeric(x) && numel(x) == 2 && isreal(x));
addParameter(p, 'anglesRCS', zeros(1,4), @(x) (isnumeric(x) && isreal(x)));
addParameter(p, 'N', 1, @(x) isnumeric(x) && mod(x,1) == 0);
addParameter(p, 'returnLargeScale', 0, @(x) isnumeric(x) );


parse(p, st, varargin{:});

st = p.Results.st;
spst = p.Results.spst;
stDirection = p.Results.stDirection;
anglesRCS = p.Results.anglesRCS/180*pi;
N =  p.Results.N;
returnLargeScale =  p.Results.returnLargeScale;

theta_i = anglesRCS(:,1);  % Incident elevation
phi_i   = anglesRCS(:,2);  % Incident azimuth
theta_s = anglesRCS(:,3);  % Scattered elevation
phi_s   = anglesRCS(:,4);  % Scattered azimuth

if all(theta_i==theta_s) && all(phi_i==phi_s)
    sensingMode = 'monostatic';
    beta = 0;
else
    sensingMode = 'bistatic';

    u_i = [sin(theta_i).*cos(phi_i), sin(theta_i).*sin(phi_i), cos(theta_i)];  % incident
    u_s = [sin(theta_s).*cos(phi_s), sin(theta_s).*sin(phi_s), cos(theta_s)];  % scattered

    u_bis = u_i + u_s;
    u_bis = u_bis ./ vecnorm(u_bis, 2,2);  % normalize

    theta = acos(u_bis(:,3));              % Zenith angle
    phi   = atan2(u_bis(:,2), u_bis(:,1));   % Azimuth angle

    % Compute dot product
    dotProd = dot(u_i', u_s');

    % Compute bistatic angle (in radians)
    beta = acos(dotProd);
end


% Check for invalid configurations
invalidConfigs = (strcmp(st, 'uav-small') && strcmp(spst, 'multi')) || ...
    (strcmp(st, 'human')     && strcmp(spst, 'multi'));

if invalidConfigs
    error('Invalid configuration: ''%s'' cannot be used with spst = ''%s''.', st, spst);
end

% Determine angular dependency
isAngularDependent = any(strcmp(st, {'uav-large', 'vehicle', 'agv'}));

if ~isAngularDependent && (~isempty(stDirection))
    warning('Direction provided, but ST is angularly independent. Inputs will be ignored.')
end
if ~isAngularDependent && (any(anglesRCS~=0))
    warning('Bistatic angle provided, but ST is angularly independent. Inputs will be ignored.')
end

if isAngularDependent &&  (isempty(stDirection))
    stDirection = [0 90];
    warning('Angular-dependent target ''%s'' is using default direction [0, 90]. Consider specifying ''direction'' explicitly.', st);
end

T = getRCSTable(st, spst);

switch sensingMode
    case 'monostatic'
        % Define example sigma values for each target type
        switch st
            case 'uav-small'
                sigmaMdBsm = -12.81;
                sigmaDdB = 0;
                sigmaSigmaSdB = 3.74;

                sigmaSdB = getAngularIndependentRCS(sigmaSigmaSdB,[1 N]);

            case 'uav-large'
                sigmaMdBsm = NaN;
                sigmaDdB = NaN;
                sigmaSdB = NaN;

            case 'human-model-1'
                sigmaMdBsm = -1.37;
                sigmaDdB = 0;
                sigmaSigmaSdB = 3.94;

                sigmaSdB = getAngularIndependentRCS(sigmaSigmaSdB,[1 N]);

            case 'agv'
                sigmaMdBsm = NaN;
                sigmaDdB = NaN;
                sigmaSdB = NaN;

            case 'vehicle'
                sigmaMdBsm = -20;
                sigmaDdB = getAngularDependentRCS(stDirection(1), stDirection(2), T);
                sigmaSigmaSdB = 3.41;
                sigmaSdB = getAngularIndependentRCS(sigmaSigmaSdB,[size(sigmaDdB,1), N]);
        end
    case 'bistatic'

        switch st
            case 'uav-small'
                sigmaMdBsm = max(-12.81 - 3*sin(beta/2), -inf);
                sigmaDdB = 0;
                sigmaSigmaSdB = 3.74;
                sigmaSdB = getAngularIndependentRCS(sigmaSigmaSdB,[1 N]);

            case 'uav-large'
                k1 = 6.05;
                k2 = 1.33;
                sigmaMdBsm = getSigmaMD(T, theta, phi, k1,k2);
                sigmaDdB = 1;
                sigmaSigmaSdB = T.sigmaSdB;
                sigmaSdB = getAngularIndependentRCS(sigmaSigmaSdB,[size(sigmaDdB,1), N]);

            case 'vehicle'
                k1 = 6.05;
                k2 = 1.33;
                sigmaMdBsm = getSigmaMD(T, theta, phi, k1,k2);
                sigmaDdB = 1;
                sigmaSigmaSdB = T.sigmaSdB;
                sigmaSdB = getAngularIndependentRCS(sigmaSigmaSdB,[size(sigmaDdB,1), N]);

            case 'agv'
                k1 = 12;
                k2 = 1.45;
                sigmaMdBsm = getSigmaMD(T, theta, phi, k1,k2);
                sigmaDdB = 1;
                sigmaSigmaSdB = T.sigmaSdB;
                sigmaSdB = getAngularIndependentRCS(sigmaSigmaSdB,[size(sigmaDdB,1), N]);

            case 'human-model-1'
                sigmaMdBsm =  max(-1.37 - 3*sin(beta/2), -inf) ;
                sigmaDdB = 0;
                sigmaSigmaSdB = 3.94;
                sigmaSdB = getAngularIndependentRCS(sigmaSigmaSdB,[1 N]);
        end
end

if returnLargeScale
    sigmaRCSdBsm = sigmaMdBsm;
else
    sigmaRCSdBsm = sigmaMdBsm+sigmaDdB+sigmaSdB;
end

end

function sigmaMDdB = getAngularDependentRCS(theta, phi, T)
%GETANGULARDEPENDENTRCS Computes the directional RCS gain component σ_MD (in dB)
%   based on incident direction (theta, phi) and a 3GPP-style RCS table.
%
%   This implements the model defined in 3GPP TR 38.901, Table 7.9.2.1-x, for
%   angular-dependent monostatic RCS. The output represents the direction-dependent
%   gain due to both elevation (θ) and azimuth (φ) deviation from sector center.
%
%   INPUTS:
%     theta - Elevation angle (degrees)
%     phi   - Azimuth angle (degrees)
%     T     - Table containing the directional RCS model parameters
%
%   OUTPUT:
%     sigma_MD_dB - Directional RCS gain in dB

if phi>=315, phi = 360-phi; end
indexSPT = (T.phi_range(:,1)<=phi & T.phi_range(:,2)>=phi) & (T.theta_range(:,1)<=theta & T.theta_range(:,2)>=theta);


% Extract parameters
thetaCenter = T.theta_center(indexSPT);
theta3dB    = T.theta_3dB(indexSPT);
phiCenter   = T.phi_center(indexSPT);
phi3dB      = T.phi_3dB(indexSPT);
gMax        = T.G_max(indexSPT);
sigmaMax    = T.sigma_max(indexSPT);

sigmaVdB = 12 * ((theta - thetaCenter)./theta3dB).^2;
sigmaVdB = -min(sigmaVdB, sigmaMax);

sigmaHdB = 12 * ((phi -phiCenter) ./ phi3dB).^2;
sigmaHdB = -min(sigmaHdB, sigmaMax);
sigmaMDdB = gMax - min(-(sigmaVdB + sigmaHdB), sigmaMax);

end

function sigmaRCSdB = getAngularIndependentRCS(sigmaSdB,N)
%GETANGULARINDEPENDENTRCS Generates N samples of angular-independent RCS in dBsm
%   based on a log-normal fluctuation model with a given mean (in dB) and standard deviation.
%
%   INPUTS:
%     sigmaMdBsm - Mean RCS value in dBsm
%     sigmaSdB   - Standard deviation of the RCS fluctuation in dB
%     N          - Number of samples to generate
%
%   OUTPUT:
%     sigmaRCSdBsm - 1×N vector of RCS samples in dBsm, where the linear-scale mean is 10^(sigmaMdBsm/10)


sigmaRCSdB =  (-log(10)/20)*sigmaSdB^2 +  sigmaSdB*randn(N);

end

function rcsTable = getRCSTable(st, spst)
%GETRCSTABLE Returns the angular-dependent RCS parameter table for a given target configuration.
%
%   rcsTable = GETRCSTABLE(st, spst, sensingMode)
%
%   INPUTS:
%     st          - Target type (e.g., 'vehicle', 'agv', 'uav-large')
%     spst        - Scattering point type: 'single' or 'multi'
%     sensingMode - Either 'monostatic' or 'bistatic' (currently only monostatic implemented)
%
%   OUTPUT:
%     rcsTable    - A MATLAB table containing angular RCS model parameters from 3GPP
%                   (e.g., Table 7.9.2.1-4 to 7.9.2.1-7 of TR 38.901)
rcsTable = [];

tol = 1e-10;
% Table 7.9.2.1-2
if strcmp(st, 'uav-large')
    rcsTable = table( ...
        ["Left"; "Back"; "Right"; "Front"; "Bottom"; "Roof"], ...     % Direction labels
        [90; 180; 270; 0; NaN; NaN], ...                                       % phi_center
        [7.13; 10.09; 7.13; 14.19; NaN; NaN], ...                                       % phi_3dB
        [90; 90; 90; 90; 180; 0], ...                                       % theta_center
        [8.68; 11.43; 8.68; 16.53; 4.93; 4.93], ...                                       % theta_3dB
        [7.43; 3.99; 7.43; 1.02; 13.55; 13.55], ...                                       % G_max
        [14.30; 10.86; 14.30; 7.89; 20.42; 20.42], ...                                       % sigma_max
        [45 135; 45 135; 45 135; 45 135; 135 180; 0 45], ...  % theta_range
        [45 135; 135 225; 225 315; -45 45; 0 360; 0 360], ...  % phi_range
        'VariableNames', { ...
        'Direction', ...
        'phi_center', 'phi_3dB', ...
        'theta_center', 'theta_3dB', ...
        'G_max', 'sigma_max', ...
        'theta_range', 'phi_range' ...
        } ...
        );

    rcsTable.sigmaMdBsm = repmat(-5.85, height(rcsTable), 1);  % fill later
    rcsTable.sigmaSdB = repmat(2.5, height(rcsTable), 1);
end

% Table 7.9.2.1-4
if strcmp(st, 'vehicle') & strcmp(spst, 'single')
    rcsTable = table( ...
        ["Left"; "Back"; "Right"; "Front"; "Roof"], ...
        [90; 180; 270; 0; NaN], ...                    % phi_center
        [26.90; 36.32; 26.90; 40.54; NaN], ...         % phi_3dB
        [79.70; 79.65; 79.70; 71.75; 0], ...         % theta_center
        [44.42; 36.73; 44.42; 29.13; 18.13], ...         % theta_3dB
        [20.75; 14.56; 20.75; 15.52; 21.26], ...         % G_max
        [13.68; 7.50; 13.68; 8.45; 14.19], ...           % sigma_max
        [30 180; 30 180; 30 180;  30 180; 0 30-tol], ...             % theta_range (numeric vector)
        [45+tol 135; 135+tol 225; 225+tol 315;  -45+tol 45; 0 360-tol], ...  % phi_range (numeric, with wrap-around for 'Front')
        'VariableNames', { ...
        'Direction', ...
        'phi_center', 'phi_3dB', ...
        'theta_center', 'theta_3dB', ...
        'G_max', 'sigma_max', ...
        'theta_range', 'phi_range' ...
        } ...
        );

    rcsTable.sigmaMdBsm =  repmat(11.25, height(rcsTable), 1);
    rcsTable.sigmaSdB = repmat(3.41, height(rcsTable), 1);

end

% Table 7.9.2.1-5
if strcmp(st, 'vehicle') & strcmp(spst, 'multi')
    rcsTable = table( ...
        ["Left"; "Back"; "Right"; "Front"; "Roof"], ...
        [90; 180; 270; 0; NaN], ...                        % phi_center
        [26.90; 36.32; 26.90; 40.54; NaN], ...             % phi_3dB
        [79.70; 79.65; 79.70; 71.75; 0.00], ...            % theta_center
        [44.42; 36.73; 44.42; 29.13; 18.13], ...           % theta_3dB
        [20.60; 13.90; 20.60; 14.99; 21.12], ...           % G_max
        [20.52; 13.82; 20.52; 14.91; 21.05], ...           % sigma_max
        [0 180; 0 180; 0 180; 0 180; 0 180], ...           % theta_range (numeric matrix: 5x2)
        [0 360; 0 360; 0 360; 0 360; 0 360], ...          % phi_range (as 5x2 matrix, wrap-around for front)
        'VariableNames', { ...
        'Direction', ...
        'phi_center', 'phi_3dB', ...
        'theta_center', 'theta_3dB', ...
        'G_max', 'sigma_max', ...
        'theta_range', 'phi_range' ...
        } ...
        );

    % Add sigmaM and sigmaS columns
    rcsTable.sigmaMdBsm = repmat(11.25, height(rcsTable), 1);        % Mean RCS (empty for now)
    rcsTable.sigmaSdB = repmat(3.41, height(rcsTable), 1); % Standard deviation of fluctuation

end

% Table 7.9.2.1-6
if strcmp(st, 'agv') & strcmp(spst, 'single')
    rcsTable = table( ...
        ["Left"; "Back"; "Right"; "Front"; "Roof"], ...     % Direction labels
        [0; 90; 180; 270; NaN], ...                                       % phi_center
        [13.68; 15.53; 12.49; 15.53; NaN], ...                                       % phi_3dB
        [90; 75; 90; 75; 0], ...                                       % theta_center
        [13.68; 20.03; 11.89; 20.03; 11.44], ...                                       % theta_3dB
        [13.02; 7.33; 11.01; 7.33; 11.79], ...                                       % G_max
        [23.29; 17.60; 21.28; 17.60; 22.06], ...                                       % sigma_max
        [30 180; 30 180; 30 180; 30 180; 0 30], ...  % theta_range
        [-45 45; 45 135; 135 225; 225 315; 0 360], ...  % phi_range
        'VariableNames', { ...
        'Direction', ...
        'phi_center', 'phi_3dB', ...
        'theta_center', 'theta_3dB', ...
        'G_max', 'sigma_max', ...
        'theta_range', 'phi_range' ...
        } ...
        );

    rcsTable.sigmaMdBsm = repmat(-4.25, height(rcsTable), 1);
    rcsTable.sigmaSdB = repmat(2.51, height(rcsTable), 1);
end

% Table 7.9.2.1-7
if strcmp(st, 'agv') & strcmp(spst, 'multi')
    rcsTable = table( ...
        ["Left"; "Back"; "Right"; "Front"; "Roof"], ...     % Direction labels
        [0; 90; 180; 270; NaN], ...                                       % phi_center
        [13.68; 15.53; 12.49; 15.53; NaN], ...                                       % phi_3dB
        [90; 75; 90; 75; 0], ...                                       % theta_center
        [13.68; 20.03; 11.89; 20.03; 11.44], ...                                       % theta_3dB
        [13; 7.27; 10.98; 7.27; 11.77], ...                                       % G_max
        [30.26; 24.53; 28.24; 24.53; 29.03], ...                                       % sigma_max
        [0 180; 0 180; 0 180; 0 180; 0 180], ...  % theta_range
        [0 360; 0 360; 0 360; 0 360; 0 360], ...  % phi_range
        'VariableNames', { ...
        'Direction', ...
        'phi_center', 'phi_3dB', ...
        'theta_center', 'theta_3dB', ...
        'G_max', 'sigma_max', ...
        'theta_range', 'phi_range' ...
        } ...
        );

end


end

function sigmaMDdBsm = getSigmaMD(T, theta, phi, k1, k2)

sigma_v_dB = -min([12 * ((theta - T.theta_center) ./ T.theta_3dB).^2, T.sigma_max],[],2);
sigma_h_dB = -min([12 * ((phi - T.phi_center) ./ T.phi_3dB).^2, T.sigma_max],[],2);
term1 = T.G_max - min([-(sigma_v_dB + sigma_h_dB), T.sigma_max], [],2);
term2 = -k1 * sin(k2 * beta / 2) + 5 * log10(cos(beta / 2));
sigmaMDdBsm = max([term1 + term2, T.G_max - T.sigma_max, -inf*ones(5,1)], [], 2);

end
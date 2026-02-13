% HYBRIDBEAMFORMER Create, visualize, and export hybrid (subarray) beamforming weights
%   This script configures a hybrid beamforming architecture for an 8-by-8
%   array partitioned into Mp-by-Np replicated subarrays using Phased Array
%   System Toolbox. A vertical ULA subarray is replicated to form a
%   ReplicatedSubarray with custom subarray steering. The script then builds
%   element-level (analog) weights by combining (i) a vertical steering phase
%   progression within each subarray and (ii) a quadratic phase taper to
%   widen the beam. Digital (RF-chain) weights are set to all ones for
%   visualization. Radiation patterns are plotted for the full 2D scan and
%   for azimuth/elevation cuts. Finally, the array and beamformer parameters
%   are packaged into a structure and exported to JSON.
%
%   Requirements:
%     - Phased Array System Toolbox
%
%   See also phased.ReplicatedSubarray, phased.ULA, phased.URA, pattern
%
%   2026 NIST/CTL Jian Wang
%
%   This file is available under the terms of the NIST License.

clear;
close all;
clc;

%% Configuration
c      = 3e8;
fc     = 4e9;
lambda = c / fc;

M  = 8;
N  = 8;
Mp = 4;
Np = 8;

numElemPerSub = M / Mp;

d_H = 0.5;
d_V = 0.8;

%% Arrays
sULA = phased.ULA( ...
    'NumElements',      numElemPerSub, ...
    'Element',          phased.NRAntennaElement, ...
    'ElementSpacing',   d_V*lambda, ...
    'ArrayAxis',        'z');

aURA = phased.ReplicatedSubarray( ...
    'Subarray',           sULA, ...
    'GridSize',           [Mp Np], ...
    'SubarraySteering',   'Custom', ...
    'GridSpacing',        [M*d_V/Mp d_H]*lambda);

dURA = phased.URA( ...
    'Size',             [8 8], ...
    'Element',          phased.NRAntennaElement, ...
    'ElementSpacing',   [d_V*lambda, d_H*lambda]);

pos = getElementPosition(dURA) / lambda;

%% Element ordering adjustment
idx = 1:64;

% Reshape to 2 rows, swap the rows, and flatten back
idx_swapped = reshape(idx, 2, []).';   % Step 1: [1 2; 3 4; ...]
idx_swapped = fliplr(idx_swapped).';   % Step 2: [2 1; 4 3; ...]
idx_swapped = idx_swapped(:);

pos = pos(:, idx_swapped);

% Define steering directions
steerAz = 0;    % Horizontal steering
steerEl = 9;    % Vertical steering  % 5 %25 %42 %-80

%% Hybrid weights
% A) Analog weights: steering the elements inside each subarray
local_pos_v = (0:M-1)' * d_V * lambda;
v_steering  = exp(1i * 2 * pi * local_pos_v * sind(steerEl) / lambda);

% B) Wide-beam shaping (quadratic phase)
alpha_h = 0.6;  % Tuning parameter for width  % 0.5
alpha_v = 0.5;

[grid_v, grid_h] = meshgrid(1:M, 1:N);
centered_v       = grid_v - (M+1)/2;
centered_h       = grid_h - (N+1)/2;
quad_phases      = alpha_v * centered_v.^2 + alpha_h * centered_h.^2;
widen_weights    = exp(1i * quad_phases(:));

% C) Combine into analog weight matrix
v_steer_m      = repmat(v_steering, 1, N);
analogWeights  = v_steer_m(:) .* widen_weights;

% Digital weights for the 32 RF chains (all ones), for plot purpose
w_digital = ones(32, 1);

% Analog weights reshaped per RF chain (2 elements per chain)
w_analog = reshape(analogWeights, 2, 32);   % size 2x32

% analogWeights = combined_steer_vec;
% w_analog      = reshape(analogWeights, 2, 32); % size 2x32
w_analog_plot = w_analog([2, 1], :);

%% Pattern plots
% Plot the pattern seen by the analog weights
figure;
pattern(aURA, fc, -180:180, -90:90, ...
    'Type',              'directivity', ...
    'Weights',           w_digital, ...
    'ElementWeights',    w_analog_plot, ...
    'PropagationSpeed',  c);

% 2D cut: azimuth (el = 0)
figure;
pattern(aURA, fc, -180:180, 0, ...
    'Type',              'directivity', ...
    'Weights',           w_digital, ...
    'ElementWeights',    w_analog);
grid on;
title('Horizontal Digital Beamwidth (RF Chain Level)');
ylabel('Directivity (dBi)');
xlabel('Azimuth Angle (deg)');

% 2D cut: elevation (az = 0)
figure;
pattern(aURA, fc, 0, -90:90, ...
    'Type',              'directivity', ...
    'Weights',           w_digital, ...
    'ElementWeights',    w_analog);
grid on;
title('Vertical Digital Beamwidth (RF Chain Level)');
ylabel('Directivity (dBi)');
xlabel('Elevation Angle (deg)');



%% Export configuration
cfg      = struct();
cfg.meta = struct("name", "hybrid", "fc_Hz", fc);
cfg.meta.array = struct( ...
    "type",       "URA", ...
    "M",          M, ...
    "Mprime",     Mp, ...
    "N",          N, ...
    "Nprime",     Np, ...
    "dV_lambda",  d_V, ...
    "dH_lambda",  d_H);

cfg.beamformer = struct();
cfg.beamformer.type   = "hybrid";
cfg.beamformer.wElem  = w_analog;   % Analog weights (2 x 32)
cfg.beamformer.Mprime = Mp;
cfg.beamformer.Nprime = Np;

% nrRadar.io.exportBFtoJSON("bf_hybrid2.json", cfg);

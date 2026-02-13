% FULLDIGITALBEAMFORMER Create, visualize, and export full-digital beamforming weights
%   This script configures a single 3GPP NR antenna element and an 8-by-8
%   uniform rectangular array (URA) using Phased Array System Toolbox. It
%   constructs full-digital beamforming weights by averaging steering vectors
%   over a specified elevation scan and a specified azimuth scan. 
%   The resulting separable 2-D weight vector is formed via a Kronecker 
%   product, normalized to unit power, and used to plot array patterns 
%   (2-D pattern and principal cuts). Finally, the array and beamformer 
%   parameters are packaged into a structure suitable for JSON export.
%
%   Requirements:
%     - Phased Array System Toolbox
%
%   See also phased.NRAntennaElement, phased.URA, phased.ULA, phased.SteeringVector, pattern
%
%   2026 NIST/CTL Jian Wang
%
%   This file is available under the terms of the NIST License.

clear;
close all;
clc;

%% Basic parameters
c      = 3e8;
fc     = 4e9;
lambda = c / fc;

% Single 3GPP NR antenna element
nRAntennaElement = phased.NRAntennaElement;

%% Element pattern sanity checks
figure;
pattern(nRAntennaElement, fc, 0, -90:90, 'CoordinateSystem', 'polar');

figure;
pattern(nRAntennaElement, fc, -180:180, 0, 'CoordinateSystem', 'polar');
% pattern(nRAntennaElement, fc);

%% Array geometry
M  = 8;
N  = 8;
Mp = 8;
Np = 8;

d_H = 0.5;  % Horizontal element spacing (in lambda)
d_V = 0.5;  % Vertical element spacing (in lambda)

%% Construct the MxN URA array
aURA = phased.URA( ...
    'Size',            [M N], ...
    'ElementSpacing',  [d_V*lambda, d_H*lambda]);
aURA.Element = nRAntennaElement;

%% Build elevation steering (average over elevation scan)
% Define 1-D vertical array (along z)
sULA = phased.ULA( ...
    'NumElements',     M, ...
    'ElementSpacing',  d_V*lambda, ...
    'ArrayAxis',       'z');
sULA.Element = nRAntennaElement;

el_angles = -60:1:60;

% Because the array is on the Z-axis, azimuth does not change the phase
% delay relative to the elements. Fix azimuth to 0.
az_angles   = zeros(size(el_angles));
scan_angles = [az_angles; el_angles];  % [Azimuth; Elevation]

sv = phased.SteeringVector( ...
    'SensorArray',        sULA, ...
    'PropagationSpeed',   c);

% Output is an M x numAngles matrix
steer_vecs = sv(fc, scan_angles);
steerEl    = mean(steer_vecs, 2);

%% Build azimuth steering (average over azimuth scan)
% Define 1-D horizontal array (along y)
sULA2 = phased.ULA( ...
    'NumElements',     N, ...
    'ElementSpacing',  d_H*lambda, ...
    'ArrayAxis',       'y');
sULA2.Element = nRAntennaElement;

az_angles = -60:1:60;

% Because the array is on the Y-axis, elevation does not change the phase
% delay relative to the elements. Fix elevation to 0.
el_angles   = zeros(size(az_angles));
scan_angles = [az_angles; el_angles];  % [Azimuth; Elevation]

sv = phased.SteeringVector( ...
    'SensorArray',        sULA2, ...
    'PropagationSpeed',   c);

% Output is an N x numAngles matrix
steer_vecs2 = sv(fc, scan_angles);
steerAz     = mean(steer_vecs2, 2);

%% Combine elevation/azimuth into a separable 2D weight vector
combined_steer_vec = kron(steerEl, steerAz);
combined_steer_vec = combined_steer_vec / norm(combined_steer_vec);

%% Pattern plots
figure(100);
pattern(aURA, fc, ...
    'PropagationSpeed',  c, ...
    'Weights',           combined_steer_vec, ...
    'Type',              'powerdb', ...
    'CoordinateSystem',  'polar', ...
    'Normalize',         false);

figure(110);
pattern(aURA, fc, -90:90, 0, ...
    'PropagationSpeed',  c, ...
    'Weights',           combined_steer_vec, ...
    'Type',              'powerdb', ...
    'CoordinateSystem',  'polar', ...
    'Normalize',         false);

figure(120);
pattern(aURA, fc, 0, -90:90, ...
    'PropagationSpeed',  c, ...
    'Weights',           combined_steer_vec, ...
    'Type',              'powerdb', ...
    'CoordinateSystem',  'polar', ...
    'Normalize',         false);

%% Export configuration
cfg      = struct();
cfg.meta = struct("name", "full_digital", "fc_Hz", fc);
cfg.meta.array = struct( ...
    "type",       "URA", ...
    "M",          M, ...
    "Mprime",     M, ...
    "N",          N, ...
    "Nprime",     N, ...
    "dV_lambda",  d_V, ...
    "dH_lambda",  d_H);

cfg.beamformer = struct();
cfg.beamformer.type   = "full_digital";
cfg.beamformer.wElem  = combined_steer_vec(:);  % 64x1 complex (canonical)
cfg.beamformer.Mprime = M;
cfg.beamformer.Nprime = N;

% nrRadar.io.exportBFtoJSON("bf_full_digital.json", cfg);

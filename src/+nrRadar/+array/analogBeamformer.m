% ANALOGBEAMFORMER Create, visualize, and export analog beamforming weights
%   This script configures a single 3GPP NR antenna element and an 8-by-8
%   uniform rectangular array (URA) using Phased Array System Toolbox. It
%   then constructs uniform (all-ones) analog beamforming weights, applies
%   the element ordering expected by phased.URA, and plots the resulting
%   radiation patterns. Finally, it packages the array and beamformer
%   configuration into a structure suitable for JSON export.
%
%
%   Requirements:
%     - Phased Array System Toolbox
%
%   See also phased.NRAntennaElement, phased.URA, pattern
%
%   2026 NIST/CTL Steve Blandino
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

%% Array geometry
M = 8;   % Vertical elements
N = 8;   % Horizontal elements

d_V = 0.5;   % Vertical spacing (in lambda)
d_H = 0.5;   % Horizontal spacing (in lambda)

% URA (8 x 8)
aURA = phased.URA( ...
    'Size',            [M N], ...
    'ElementSpacing',  [d_V*lambda, d_H*lambda]);
aURA.Element = nRAntennaElement;

%% Beamforming weights (analog/uniform)
F_RF    = ones(N*M, 1);
wBB     = 1;
wHybrid = F_RF * wBB;
wHybrid = wHybrid / norm(wHybrid);   % Normalize total power

% Element re-ordering for URA column-major ordering expected by phased.URA
wHybrid_rm = wHybrid;
W          = reshape(wHybrid_rm, [N, M]);
wHybrid    = W(:);

%% ------------------------------------------------------------------------
% Check elevation/azimuth patterns
%% ------------------------------------------------------------------------

% 2D pattern
figure(101);
pattern(aURA, fc, ...
    'PropagationSpeed',  c, ...
    'Weights',           wHybrid, ...
    'Type',              'powerdb', ...
    'CoordinateSystem',  'polar', ...
    'Normalize',         false);
title('Hybrid 8\times8 with 2-element vertical subarrays (2D)');

% Azimuth cut (elevation = 0 deg)
figure(110);
pattern(aURA, fc, -90:90, 0, ...
    'PropagationSpeed',  c, ...
    'Weights',           wHybrid, ...
    'Type',              'powerdb', ...
    'CoordinateSystem',  'polar', ...
    'Normalize',         false);
title('Hybrid: Azimuth cut (el = 0^\circ)');

% Elevation cut (azimuth = 0 deg)
figure(120);
pattern(aURA, fc, 0, -90:90, ...
    'PropagationSpeed',  c, ...
    'Weights',           wHybrid, ...
    'Type',              'powerdb', ...
    'CoordinateSystem',  'polar', ...
    'Normalize',         false);
title('Hybrid: Elevation cut (az = 0^\circ)');

%% Export configuration
cfg      = struct();
cfg.meta = struct("name", "analog", "fc_Hz", fc);
cfg.meta.array = struct( ...
    "type",       "URA", ...
    "M",          M, ...
    "Mprime",     1, ...
    "N",          N, ...
    "Nprime",     1, ...
    "dV_lambda",  d_V, ...
    "dH_lambda",  d_H);

cfg.beamformer = struct();
cfg.beamformer.type   = "analog";
cfg.beamformer.wElem  = wHybrid(:);   % 64x1 complex (canonical)
cfg.beamformer.Mprime = 1;
cfg.beamformer.Nprime = 1;

% nrRadar.io.exportBFtoJSON("bf_analog.json", cfg);

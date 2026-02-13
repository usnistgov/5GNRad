function R = evalTP_PosVel(gtPos, estPos, gtVel, estVel, sensorPos)
% EVALTP_POSVEL Evaluate position and velocity errors for target estimates
%   R = EVALTP_POSVEL(GTPOS, ESTPOS, GTVEL, ESTVEL, SENSORPOS) computes 
%   position, range, angular, and velocity errors between ground-truth and 
%   estimated target states relative to a reference sensor position.
%
%   Inputs
%   ------
%   GTPOS     : [N x 3] Ground-truth Cartesian positions [x y z] of targets
%   ESTPOS    : [N x 3] Estimated Cartesian positions [x y z] of targets
%   GTVEL     : [N x 1] Ground-truth radial velocity (or scalar velocity)
%   ESTVEL    : [N x 1] Estimated radial velocity (or scalar velocity)
%   SENSORPOS : [1 x 3] Sensor position [x y z] used as reference origin
%
%   Outputs
%   -------
%   R : struct containing position and velocity error metrics
%       R.pos.errXYZ     : [N x 3] Cartesian position error components
%       R.pos.absXYZ     : [N x 3] Absolute value of position errors
%       R.pos.errMag     : [N x 1] Magnitude of position error
%       R.pos.range_gt   : [N x 1] Ground-truth range from sensor
%       R.pos.range_est  : [N x 1] Estimated range from sensor
%       R.pos.range_err  : [N x 1] Range error (est - gt)
%       R.pos.az_err_deg : [N x 1] Azimuth angle error in degrees
%       R.pos.el_err_deg : [N x 1] Elevation angle error in degrees
%       R.vel.vr_err     : [N x 1] Velocity error (est - gt)
%
%   * This function assumes consistent ordering between ground-truth and 
%     estimated entries.
%
%   Example
%       gtPos = [0 0 0; 1 1 0];
%       estPos = [0 0 0.1; 1.1 1 0];
%       gtVel = [1.0; 2.0];
%       estVel = [1.2; 1.9];
%       sensorPos = [0 0 0];
%       R = evalTP_PosVel(gtPos, estPos, gtVel, estVel, sensorPos);
%
%   2025 NIST/CTL Steve Blandino
%
%   This file is available under the terms of the NIST License.


    % -------- Position errors --------
    dPos        = estPos - gtPos;                 % [N x 3]
    R.pos.errXYZ = dPos;
    R.pos.absXYZ = abs(dPos);
    R.pos.errMag = vecnorm(dPos, 2, 2);

    % Ranges from sensor
    r_gt_vec     = gtPos  - sensorPos;            % [N x 3]
    r_est_vec    = estPos - sensorPos;            % [N x 3]

    vx = r_gt_vec(:,1);
    vy = r_gt_vec(:,2);
    vz = r_gt_vec(:,3);

    % Compute Azimuth AoA (in degrees)
    az_gt = atan2d(vy, vx);

    % Compute Elevation AoA (in degrees)
    el_gt = atan2d(vz, sqrt(vx.^2 + vy.^2));

    vx = r_est_vec(:,1);
    vy = r_est_vec(:,2);
    vz = r_est_vec(:,3);

    % Compute Azimuth AoA (in degrees)
    az_est = atan2d(vy, vx);

    % Compute Elevation AoA (in degrees)
    el_est = atan2d(vz, sqrt(vx.^2 + vy.^2));

    R.pos.range_gt  = vecnorm(r_gt_vec,  2, 2);
    R.pos.range_est = vecnorm(r_est_vec, 2, 2);
    R.pos.range_err = R.pos.range_est - R.pos.range_gt;

    az_err = wrapTo180(az_est - az_gt);
    el_err = wrapTo180(el_est - el_gt);

    R.pos.az_err_deg = az_err;
    R.pos.el_err_deg = el_err;

    dVel         = estVel(:) - gtVel(:);

    R.vel.vr_err = dVel;
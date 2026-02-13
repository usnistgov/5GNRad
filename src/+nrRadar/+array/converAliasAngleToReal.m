function theta_real = converAliasAngleToReal(theta_alias, d_lambda)
%CONVERALIASANGLETORREAL Convert aliased angle to real physical angle.
%   THETA_REAL = CONVERALIASANGLETORREAL(THETA_ALIAS, D_LAMBDA) converts an
%   observed (aliased) angle THETA_ALIAS (in degrees) into a "real" physical
%   angle THETA_REAL (in degrees), based on the normalized inter-element
%   spacing D_LAMBDA = d/lambda.
%
%   The function maps the aliased angle to direction sine u = sin(theta),
%   applies an unaliasing shift u_real = u_alias + 1/D_LAMBDA, and then maps
%   back to angle. If the resulting direction sine is not physically valid
%   (|u_real| > 1), the function returns THETA_ALIAS unchanged.
%
%   Inputs
%     THETA_ALIAS : Aliased angle in degrees (scalar or array).
%     D_LAMBDA    : Normalized spacing d/lambda (positive scalar).
%
%   Output
%     THETA_REAL  : Unaliased (real) angle in degrees (same size as THETA_ALIAS).
%
%   Example
%     % For d/lambda = 0.8, shift = 1/0.8 = 1.25
%     thetaReal = converAliasAngleToReal(-30, 0.8);
%
%   2026 NIST/CTL Jian Wang
%
%   This file is available under the terms of the NIST License.


    % Convert observed angle to Direction Sin (u)
    u_alias = sin(deg2rad(theta_alias));
    
    % Calculate the shift factor (L/d)
    % For d = 0.8, shift is 1.25
    shift = 1 / d_lambda;
    
    % Shift the alias back to the positive domain
    % We add the shift to move from the negative side to the positive side
    u_real = u_alias + shift;
    
    % Check for physical validity
    if abs(u_real) > 1
        theta_real = theta_alias;
    else
        theta_real = rad2deg(asin(u_real));
    end
    
end
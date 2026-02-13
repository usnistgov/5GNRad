function assumeToolboxes(testCase)
%ASSUMETOOLBOXES Skip tests if required MathWorks products are unavailable.
%
% The 5GNRad pipeline uses:
%   * 5G Toolbox: nrCarrierConfig, nrPRSConfig, nrOFDMInfo, nrOFDMModulate
%   * Phased Array System Toolbox: phased.URA, phased.SteeringVector
%
% We check for key symbols instead of license feature names.

has5G = (exist('nrCarrierConfig','class')==8) && ...
        (exist('nrPRSConfig','class')==8) && ...
        (exist('nrOFDMInfo','file')==2) && ...
        (exist('nrOFDMModulate','file')==2);

hasPhased = (exist('phased.URA','class')==8) && ...
            (exist('phased.SteeringVector','class')==8);

testCase.assumeTrue(has5G && hasPhased, ...
    'Required toolbox missing. These tests need 5G Toolbox and Phased Array System Toolbox.');
end

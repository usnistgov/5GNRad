function pwr = getTargetPower(targetChannel)
% GETTARGETPOWER Compute summed target-channel power from path gains
%   PWR = GETTARGETPOWER(TARGETCHANNEL) returns the summed power for each
%   entry of the struct array TARGETCHANNEL.
%
%   Each element TARGETCHANNEL(k) is expected to include the field
%   AveragePathGains, containing per-path gains in linear units (not dB).
%   The output PWR is a K-by-1 vector, where K = numel(TARGETCHANNEL), and
%   PWR(k) is computed from all elements of TARGETCHANNEL(k).AveragePathGains
%   as:
%       PWR(k) = sum(AveragePathGains(:))
%
%   If AveragePathGains is empty for an entry, the corresponding power is
%   set to 0.
%
%   2026 NIST/CTL Steve Blandino
%
%   This file is available under the terms of the NIST License.

K = numel(targetChannel);

pwr = zeros(K,1);

for k = 1:K
    g = targetChannel(k).AveragePathGains;
    if isempty(g)
        pwr(k) = 0;
    else
        pwr(k) = sum(g(:));
    end
end


end
function y1 = prsDestaggering(y, prs,l, varargin)
% PRSDESTAGGERING Reorders staggered PRS subcarriers across symbols
%   Y1 = PRSDESTAGGERING(Y, PRS, L) destaggers the received PRS subcarriers
%   for the specified PRS symbols based on the comb pattern and
%   resource element (RE) offset defined in the PRS structure.
%
%   Inputs:
%     Y   - [N x S] matrix of received PRS signals in the frequency domain,
%           where N is the number of PRS subcarriers per symbol and S is
%           the number of PRS symbols.
%     PRS - structure with fields:
%           * SymbolStart     - Index of first PRS symbol
%           * CombSize        - Comb spacing (e.g., 2, 4, 6, 12)
%           * REOffset        - Subcarrier offset within each PRS symbol
%           * SymbolDuration  - (Unused) Duration of one OFDM symbol [s]
%     L   - Vector of PRS symbol indices within the slot (1-based)
%
%   Output:
%     Y1  - Column vector of destaggered PRS samples (subcarriers reordered)
%
%

%   2025 NIST/CTL Steve Blandino

%   This file is available under the terms of the NIST License.


% Varargin
lPRSStart = prs.SymbolStart;
kPRSComb = prs.CombSize;
kPRSOffset = prs.REOffset;

% Constants
kPrimeValsTab = [2   0   1   0   1   0   1   0   1   0   1   0   1;
    4   0   2   1   3   0   2   1   3   0   2   1   3;
    6   0   3   1   4   2   5   0   3   1   4   2   5;
    12  0   6   3   9   1   7   4   10  2   8   5   11];

% Find kPrime

rowIndex  = (kPrimeValsTab(:,1) == kPRSComb);
thisRow   = kPrimeValsTab(rowIndex, :);
symbolIdx = l - lPRSStart + 1;  % which symbol offset from the start of PRS
kPrimeVal = thisRow(symbolIdx + 1);  % +1 because col 1 is comb size
foffset  = mod(kPRSOffset + kPrimeVal, kPRSComb);
[~,i] = sort(foffset);
y1 = reshape(y(:,i).', [],1);
end








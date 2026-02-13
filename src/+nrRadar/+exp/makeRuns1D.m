function runs = makeRuns1D(path, values, varargin)
% MAKERUNS1D Convenience helper to build a 1-D sweep run list.
%
%   RUNS = MAKERUNS1D(PATH, VALUES) creates a struct array RUNS where each
%   element corresponds to one value in VALUES. Each run has fields:
%     - tag : string used for folder names / identification
%     - set : N-by-2 cell array of {PATH, VALUE} pairs
%
%   NAME-VALUE
%     'TagPrefix' - prefix used to build tags (default: 'run')
%     'TagFormat' - sprintf format applied to the run index (default: '%03d')
%
%   Example
%     thr_dB = 0:30;
%     thr = 10.^ (thr_dB/10);
%     runs = nrRadar.exp.makeRuns1D('sens.cfarThreshold', thr, 'TagPrefix','cfar');
%
%   2026 NIST/CTL Steve Blandino
%   This file is available under the terms of the NIST License.

p = inputParser;
addParameter(p, 'TagPrefix', 'run');
addParameter(p, 'TagFormat', '%03d');
parse(p, varargin{:});
opt = p.Results;

values = values(:);
n = numel(values);
runs = repmat(struct('tag', "", 'set', []), n, 1);

for k = 1:n
    runs(k).tag = string(opt.TagPrefix) + "_" + string(sprintf(opt.TagFormat, k));
    runs(k).set = {path, values(k)};
end
end

function exportResults(scenarioNameStr, results, detStats, detectionOutput, varargin)
% EXPORTRESULTS Write simulation outputs to disk.
%   EXPORTRESULTS(SCENARIONAMESTR, RESULTS, DETSTATS, DETECTIONOUTPUT)
%   writes CSV/JSON files to <scenario>/Output, replacing any existing
%   Output folder.
%
%   Name-Value Pair Arguments
%     'OutputPath' - Destination folder (default: fullfile(scenarioNameStr,'Output'))
%     'Overwrite'  - If true, deletes existing OutputPath before writing (default: true)
%
%   2025-2026 NIST/CTL Steve Blandino
%   This file is available under the terms of the NIST License.

p = inputParser;
addParameter(p, 'OutputPath', "", @(x) ischar(x) || isstring(x));
addParameter(p, 'Overwrite', true, @(x) islogical(x) || isnumeric(x));
parse(p, varargin{:});
opt = p.Results;

fprintf('Storing Results\n');

resultsTab    = struct2table(results);
detStatsTab   = struct2table(detStats);
detectionJson = jsonencode(detectionOutput, 'PrettyPrint', true);

if strlength(string(opt.OutputPath)) == 0
    outputPath = fullfile(scenarioNameStr, 'Output');
else
    outputPath = char(opt.OutputPath);
end

if isfolder(outputPath)
    if logical(opt.Overwrite)
        rmdir(outputPath, 's');
        mkdir(outputPath);
    end
else
    mkdir(outputPath);
end

writetable(resultsTab,  fullfile(outputPath, 'error.csv'));
writetable(detStatsTab, fullfile(outputPath, 'detStats.csv'));

fid = fopen(fullfile(outputPath, 'detection.json'), 'w');
fwrite(fid, detectionJson, 'char');
fclose(fid);

end

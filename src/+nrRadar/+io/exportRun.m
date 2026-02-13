function exportRun(outputPath, results, detStats, detectionOutput, varargin)
% EXPORTRUN Export results for a single run to disk.
%
%   EXPORTRUN(OUTPUTPATH, RESULTS, DETSTATS, DETECTIONOUTPUT) writes:
%     - error.csv      (RESULTS struct2table)
%     - detStats.csv   (DETSTATS struct2table)
%     - detection.json (JSON of DETECTIONOUTPUT)
%
%   NAME-VALUE
%     'Overwrite' - true/false. If true and OUTPUTPATH exists, it is removed
%                   before exporting. (default: true)
%
%   2026 NIST/CTL Steve Blandino
%   This file is available under the terms of the NIST License.

p = inputParser;
addParameter(p, 'Overwrite', true);
parse(p, varargin{:});
opt = p.Results;

if isstring(outputPath)
    outputPath = char(outputPath);
end

if ~isfolder(outputPath)
    mkdir(outputPath);
else
    if opt.Overwrite
        rmdir(outputPath, 's');
        mkdir(outputPath);
    end
end

resultsTab = struct2table(results);
detStatsTab = struct2table(detStats);

writetable(resultsTab, fullfile(outputPath, 'error.csv'));
writetable(detStatsTab, fullfile(outputPath, 'detStats.csv'));

try
    detectionOutputJsonStr = jsonencode(detectionOutput, 'PrettyPrint', true);
catch
    % Fallback if PrettyPrint is not supported or data contains unsupported types
    detectionOutputJsonStr = jsonencode(detectionOutput);
end

fid = fopen(fullfile(outputPath, 'detection.json'), 'w');
fwrite(fid, detectionOutputJsonStr, 'char');
fclose(fid);

end

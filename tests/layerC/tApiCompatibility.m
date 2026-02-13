classdef tApiCompatibility < nrRadarTest.BaseTestCase
    % Layer C compatibility gate: prevent accidental public API signature changes.
    %
    % Policy:
    % - Public API is defined as functions under src/+nrRadar excluding +internal.
    % - A change in nargin/nargout profile is considered breaking unless the
    %   baseline is intentionally updated.
    %
    % To update baseline intentionally:
    %   addpath("tests"); generatePublicApiBaseline
    %
    methods (Test, TestTags=["compat","layerC"])
        function publicApiSignaturesMatchBaseline(tc)
            repoRoot = fileparts(tc.RepoRoot);
            baselineFile = fullfile(repoRoot, "tests", "compat", "publicApiSignatures.json");
            tc.assertTrue(isfile(baselineFile), ...
                "Missing API baseline. Run tests/compat/generatePublicApiBaseline.m and commit the JSON.");

            txt = fileread(baselineFile);
            baseline = jsondecode(txt);

            srcRoot = fullfile(repoRoot, "src", "+nrRadar");
            currentApi = nrRadarTest.scanPublicApiSignatures(srcRoot);

            [ok, report] = localCompareApis(baseline.api, currentApi);

            tc.verifyTrue(ok, report);
        end
    end
end

function [ok, report] = localCompareApis(base, cur)
ok = true;
report = "";

baseNames = string({base.name});
curNames  = string({cur.name});

missing = setdiff(baseNames, curNames);
added   = setdiff(curNames, baseNames);

breaks = strings(0);
common = intersect(baseNames, curNames);

for i = 1:numel(common)
    nm = common(i);
    b = base(baseNames==nm);
    c = cur(curNames==nm);

    if ~isequaln([b.narginMin b.narginMax b.nargoutMin b.nargoutMax], ...
                [c.narginMin c.narginMax c.nargoutMin c.nargoutMax])
        breaks(end+1) = sprintf("%s: baseline in/out=(%s..%s,%s..%s) current in/out=(%s..%s,%s..%s)", ...
            nm, num2str(b.narginMin), num2str(b.narginMax), num2str(b.nargoutMin), num2str(b.nargoutMax), ...
                num2str(c.narginMin), num2str(c.narginMax), num2str(c.nargoutMin), num2str(c.nargoutMax)); %#ok<AGROW>
    end
end

if ~isempty(missing) || ~isempty(breaks)
    ok = false;
    report = "Public API compatibility check failed.\n";
    if ~isempty(missing)
        report = report + "Removed public functions:\n  - " + join(missing, "\n  - ") + "\n";
    end
    if ~isempty(breaks)
        report = report + "Signature changes:\n  - " + join(breaks, "\n  - ") + "\n";
    end
    if ~isempty(added)
        report = report + "New public functions (non-breaking, informational):\n  - " + join(added, "\n  - ") + "\n";
    end
    report = report + "\nIf changes are intentional, update baseline with generatePublicApiBaseline.m and commit.";
else
    % allow additions without failing
    report = "";
end
end

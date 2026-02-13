classdef BaseTestCase < matlab.unittest.TestCase
    %BASETESTCASE Common setup utilities for nrRadar tests.

    properties
        RepoRoot (1,1) string
    end

    methods (TestMethodSetup)
        function setRepoRoot(tc)
            tc.RepoRoot = string(fileparts(fileparts(mfilename('fullpath')))); % tests/.. -> repo
        end
    end

    methods
        function assumeToolbox(tc, licenseName, friendlyName)
            if nargin < 3 || strlength(friendlyName)==0
                friendlyName = licenseName;
            end
            tc.assumeTrue(license('test', licenseName), "Missing required toolbox/license: " + friendlyName);
        end
    end
end

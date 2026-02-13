function results = runTests(varargin)
%RUNTESTS Project test entrypoint.
%
%   results = runTests() runs all tests under ./tests.
%   results = runTests('Tag','layerB') runs only tests with that tag.
%
% Usage from repo root:
%   setup; results = runTests;
%
% CI usage:
%   matlab -batch "cd('repoRoot'); setup; r = runTests; assert(all([r.Passed]));"

p = inputParser;
addParameter(p,'Tag','');
parse(p,varargin{:});
opt = p.Results;

import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.plugins.DiagnosticsValidationPlugin
import matlab.unittest.plugins.FailOnWarningsPlugin
import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.XMLPlugin

repoRoot = fileparts(mfilename('fullpath'));
repoRoot = fileparts(repoRoot); % tests/.. -> repo root

if strlength(opt.Tag) > 0
    suite = TestSuite.fromFolder(fullfile(repoRoot,'tests'), 'IncludingSubfolders', true, 'Tag', char(opt.Tag));
else
    suite = TestSuite.fromFolder(fullfile(repoRoot,'tests'), 'IncludingSubfolders', true);
end

runner = TestRunner.withTextOutput('Verbosity', 2);
runner.addPlugin(DiagnosticsValidationPlugin);
runner.addPlugin(FailOnWarningsPlugin);

% Code coverage (folder-based, avoids needing per-file list)
runner.addPlugin(CodeCoveragePlugin.forFolder(fullfile(repoRoot,'src')));

% JUnit XML output if running in CI (optional)
xmlOut = fullfile(repoRoot,'testResults.xml');
runner.addPlugin(XMLPlugin.producingJUnitFormat(xmlOut));

results = runner.run(suite);
end

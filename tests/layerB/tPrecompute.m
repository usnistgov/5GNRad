classdef tPrecompute < matlab.unittest.TestCase
    % Layer B contract tests for nrRadar.internal.precompute

    methods(TestClassSetup)
        function addProjectToPath(testCase)
            % setup();
            nrRadarTest.assumeToolboxes(testCase);
        end
    end

    methods(Test, TestTags={'layerB'})
        function returnsRequiredFields_fullDigital(testCase)
            [sim, prs, sens, geom] = nrRadarTest.makeTinyConfig('full_digital');
            [st, tx] = nrRadar.internal.precompute(sim, prs, sens, geom);

            % Basic type/shape contracts
            testCase.verifyClass(st, 'struct');
            testCase.verifyClass(tx, 'double');
            testCase.verifyEqual(size(tx,2), 1, 'txWaveform must be a column vector.');
            testCase.verifyTrue(~isempty(tx) && any(tx~=0), 'txWaveform should be non-zero.');

            % Required fields used downstream
            mustHave = {'ofdmFftLen','cpLengths','symbolIndices','prs','numberSensingSymbols', ...
                'dopplerFftLen','rangeFFTLen','rangeBinDestgrd','velocityBin','snrvar','nChanRx'};
            for k = 1:numel(mustHave)
                testCase.verifyTrue(isfield(st, mustHave{k}), "Missing field: " + mustHave{k});
            end

            % Dimensions consistent
            testCase.verifyEqual(st.nChanRx, sim.rxAntenna.meta.array.M * sim.rxAntenna.meta.array.N);
            testCase.verifyEqual(st.dopplerFftLen, sens.dopplerFftLen);
            testCase.verifyEqual(st.numberSensingSymbols, sens.numberSensingSymbols);
        end

        function returnsRequiredFields_hybrid(testCase)
            [sim, prs, sens, geom] = nrRadarTest.makeTinyConfig('hybrid');
            [st, tx] = nrRadar.internal.precompute(sim, prs, sens, geom);

            testCase.verifyClass(st, 'struct');
            testCase.verifyEqual(size(tx,2), 1);

            Mprime = sim.rxAntenna.meta.array.Mprime;
            Nprime = sim.rxAntenna.meta.array.Nprime;
            testCase.verifyEqual(st.nChanRx, Mprime*Nprime);
        end
    end
end

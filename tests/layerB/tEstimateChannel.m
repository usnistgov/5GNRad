classdef tEstimateChannel < matlab.unittest.TestCase
    % Layer B contract tests for nrRadar.rx.estimateChannel

    methods(TestClassSetup)
        function addProjectToPath(testCase)
            % setup();
            nrRadarTest.assumeToolboxes(testCase);
        end
    end

    methods(Test, TestTags={'layerB'})
        function channelEstimate_hasExpectedSchema_fullDigital(testCase)
            [sim, prs, sens, geom] = nrRadarTest.makeTinyConfig('full_digital');
            [st, tx] = nrRadar.internal.precompute(sim, prs, sens, geom);
            st.snrvar = Inf; % deterministic (no noise)

            [rxArray, ~] = nrRadarTest.buildReceiveArray(sim);
            [H, channelLen, ~] = nrRadarTest.buildSinglePathChannel(st, sim, rxArray, ...
                'tapDelaySamples', 4, 'dopplerBin', 2, 'azDeg', 15, 'elDeg', 0, 'amplitude', 1);

            rx = nrRadar.rx.getRxWaveform(tx, H, st);

            hEst = nrRadar.rx.estimateChannel(rx, st, 'channelLength', channelLen);

            testCase.verifyClass(hEst, 'cell');
            testCase.verifyEqual(numel(hEst), st.nChanRx);

            nSymTot = st.carrier.SymbolsPerSlot * st.numSlots;
            for r = 1:st.nChanRx
                testCase.verifyTrue(issparse(hEst{r}), 'Expected sparse output per channel.');
                testCase.verifySize(hEst{r}, [st.ofdmFftLen nSymTot]);
                testCase.verifyFalse(any(isnan(nonzeros(hEst{r}))), 'Channel estimate contains NaNs.');
            end

            % Nonzero columns should only occur at PRS symbols
            cols = find(any(hEst{1} ~= 0, 1));
            testCase.verifyTrue(all(ismember(cols, st.symbolIndices)));
        end
    end
end

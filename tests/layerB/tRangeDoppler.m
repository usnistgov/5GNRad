classdef tRangeDoppler < matlab.unittest.TestCase
    % Layer B functional tests for nrRadar.sens.getRangeDoppler

    methods(TestClassSetup)
        function addProjectToPath(testCase)
            % setup();
            nrRadarTest.assumeToolboxes(testCase);
        end
    end

    methods(Test, TestTags={'layerB'})
        function peakAppearsAtExpectedDoppler_fullDigital(testCase)
            [sim, prs, sens, geom] = nrRadarTest.makeTinyConfig('full_digital');
            [st, tx] = nrRadar.internal.precompute(sim, prs, sens, geom);
            st.snrvar = Inf;

            [rxArray, ~] = nrRadarTest.buildReceiveArray(sim);
            [H, channelLen, truth] = nrRadarTest.buildSinglePathChannel(st, sim, rxArray, ...
                'tapDelaySamples', 6, 'dopplerBin', 2, 'azDeg', 10, 'elDeg', 0, 'amplitude', 1);

            rx = nrRadar.rx.getRxWaveform(tx, H, st);
            hEst = nrRadar.rx.estimateChannel(rx, st, 'channelLength', channelLen);
            rdCube = nrRadar.sens.getRangeDoppler(hEst, st);

            % 2D power projection across channels
            rdPwr = sum(abs(rdCube).^2, 3);
            [~, idx] = max(rdPwr(:));
            [rMax, dMax] = ind2sub(size(rdPwr), idx);

            testCase.verifyEqual(dMax, truth.dopplerIndexShifted, ...
                'Doppler peak moved (possible regression in slow-time indexing or windowing).');

            % Range location can smear due to limited PRS bandwidth; keep tolerant.
            testCase.verifyLessThanOrEqual(abs(rMax - truth.rangeBin), 3);
        end
    end
end

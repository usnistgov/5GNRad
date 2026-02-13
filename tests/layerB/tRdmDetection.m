classdef tRdmDetection < matlab.unittest.TestCase
    % Layer B end-to-end test for nrRadar.sens.rdmDetection.

    methods(TestClassSetup)
        function addProjectToPath(testCase)
            % setup();
            nrRadarTest.assumeToolboxes(testCase);
        end
    end

    methods(Test, TestTags={'layerB'})
        function detectsSingleTarget_fullDigital(testCase)
            [sim, prs, sens, geom] = nrRadarTest.makeTinyConfig('full_digital');
            [st, tx] = nrRadar.internal.precompute(sim, prs, sens, geom);
            st.snrvar = Inf;

            [rxArray, ~] = nrRadarTest.buildReceiveArray(sim);
            [H, channelLen, truth] = nrRadarTest.buildSinglePathChannel(st, sim, rxArray, ...
                'tapDelaySamples', 6, 'dopplerBin', 2, 'azDeg', 10, 'elDeg', 0, 'amplitude', 5);

            rx = nrRadar.rx.getRxWaveform(tx, H, st);
            hEst = nrRadar.rx.estimateChannel(rx, st, 'channelLength', channelLen);
            rdCube = nrRadar.sens.getRangeDoppler(hEst, st);

            [stPos, velPar, peakRda, peakSNR] = nrRadar.sens.rdmDetection( ...
                rdCube, sens, sim, rxArray, 0, st, 'angleEstimate', [NaN NaN]);

            testCase.verifyNotEmpty(stPos, 'No detections returned.');
            testCase.verifySize(stPos, [size(stPos,1) 3]);
            testCase.verifyNotEmpty(velPar);

            % Closest detected velocity should match the injected Doppler bin.
            [~,i] = min(abs(velPar - truth.velocityExpected));
            testCase.verifyLessThanOrEqual(abs(velPar(i) - truth.velocityExpected), st.prsVelocityResolution);

            testCase.verifyGreaterThan(max(peakRda), 0);
            testCase.verifyGreaterThan(max(peakSNR), 0);
        end

        function detectsSingleTarget_hybridBartlett(testCase)
            [sim, prs, sens, geom] = nrRadarTest.makeTinyConfig('hybrid');
            [st, tx] = nrRadar.internal.precompute(sim, prs, sens, geom);
            st.snrvar = Inf;

            [rxArray, ~] = nrRadarTest.buildReceiveArray(sim);
            [H, channelLen, truth] = nrRadarTest.buildSinglePathChannel(st, sim, rxArray, ...
                'tapDelaySamples', 4, 'dopplerBin', 1, 'azDeg', 20, 'elDeg', 0, 'amplitude', 5);

            rx = nrRadar.rx.getRxWaveform(tx, H, st);
            hEst = nrRadar.rx.estimateChannel(rx, st, 'channelLength', channelLen);
            rdCube = nrRadar.sens.getRangeDoppler(hEst, st);

            [stPos, velPar] = nrRadar.sens.rdmDetection( ...
                rdCube, sens, sim, rxArray, 0, st, 'angleEstimate', [NaN NaN]);

            testCase.verifyNotEmpty(stPos, 'No detections returned (hybrid/Bartlett).');
            [~,i] = min(abs(velPar - truth.velocityExpected));
            testCase.verifyLessThanOrEqual(abs(velPar(i) - truth.velocityExpected), st.prsVelocityResolution);
        end
    end
end

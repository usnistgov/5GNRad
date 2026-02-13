classdef tGetRxWaveform < matlab.unittest.TestCase
    % Layer B unit tests for nrRadar.rx.getRxWaveform

    methods(TestClassSetup)
        function addProjectToPath(testCase)
            %#ok<MANU>
            % setup();
            nrRadarTest.assumeToolboxes(testCase); %#ok<*MCNPN>
        end
    end

    methods(Test, TestTags={'layerB'})
        function matchesTimeInvariantConvolution(testCase)
            % Build a toy staticParams with constant CP to make symbol striding exact.
            ofdmFftLen = 64;
            cpLen      = 16;
            numSensSym = 4;
            numPRSSym  = 1;
            nChanRx    = 3;

            st = nrRadarTest.makeToyStaticParams( ...
                'ofdmFftLen', ofdmFftLen, ...
                'cpLen',      cpLen, ...
                'nChanRx',    nChanRx, ...
                'numSensSym', numSensSym, ...
                'numPRSSym',  numPRSSym);

            totalSyms = numSensSym * numPRSSym;
            symbolLen = ofdmFftLen + cpLen;

            % Deterministic tx waveform
            rng(42);
            txWaveform = complex(randn(totalSyms*symbolLen,1), randn(totalSyms*symbolLen,1));

            % Time-invariant channel impulse response
            h = [0.8+0.1j; 0.3-0.2j; 0.1+0.05j];
            channelLen = numel(h);

            H = complex(zeros(channelLen, totalSyms, nChanRx));
            for s = 1:totalSyms
                for r = 1:nChanRx
                    H(:,s,r) = h;
                end
            end

            rx = nrRadar.rx.getRxWaveform(txWaveform, H, st);

            % Expected: per-RF-chain identical LTI convolution
            expRx = zeros(numel(txWaveform)+channelLen-1, nChanRx);
            y = conv(txWaveform, h);
            for r = 1:nChanRx
                expRx(:,r) = y;
            end

            testCase.verifySize(rx, size(expRx));
            testCase.verifyEqual(rx, expRx, 'AbsTol', 1e-12);
        end
    end
end

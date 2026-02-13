function [results, detStats, detectionOutput, info] = run(simConfig, ...
    stConfig, prsConfig, geometry, sensConfig,backgroundChannel,targetChannel,desiredWorkers, varargin)
% RUN5GNRAD Run PRS-based radar simulation using 5G NR PRS waveforms.
%   [RESULTS, DETSTATS, DETECTIONOUTPUT, INFO] = RUN(SIMCONFIG, ...
%   STCONFIG, PRSCONFIG, GEOMETRY, SENSCONFIG, BACKGROUNDCHANNEL, ...
%   TARGETCHANNEL, DESIREDWORKERS) simulates a monostatic PRS-based radar
%   sensing pipeline over multiple sensing drops. The function precomputes
%   static PRS/OFDM parameters and a transmit waveform, partitions the 
%   target into drops, and runs one worker per drop using
%   either a serial FOR loop or a PARFOR loop (when enabled). Each drop is
%   processed by RUN5GNRAD_DROPWORKER, and results are aggregated into
%   outputs.
%
%   Inputs
%     SIMCONFIG          - Structure with simulation configuration. Must
%                          include:
%                          * nStDrop : number of state samples per drop
%     STCONFIG           - Structure with true target state time series:
%                          * position : [N x 3] positions (x,y,z)
%                          * velocity : [N x 3] velocities (vx,vy,vz)
%     PRSCONFIG          - Structure with PRS configuration used by the
%                          waveform generator and receiver processing.
%     GEOMETRY           - Structure describing scenario geometry (e.g.,
%                          TX/RX poses, boresight, etc.) as expected by
%                          RUN5GNRAD_PRECOMPUTE / RUN5GNRAD_DROPWORKER.
%     SENSCONFIG         - Structure with sensing configuration (e.g.,
%                          FFT sizes, thresholds, DOA settings) as expected
%                          by RUN5GNRAD_PRECOMPUTE / RUN5GNRAD_DROPWORKER.
%     BACKGROUNDCHANNEL  - Background channel model per drop. If empty, the
%                          simulation runs without background channel.
%                          If provided, BACKGROUNDCHANNEL(q) is passed to
%                          drop q.
%     TARGETCHANNEL      - Target channel model as a struct array aligned
%                          with STCONFIG samples. If empty, target channel
%                          is omitted. If provided, entries i1:i2 are
%                          passed to drop q.
%     DESIREDWORKERS     - Desired number of parallel workers. If empty or
%                          not provided, the function auto-selects the
%                          number of workers when parallel execution is
%                          enabled.
%
%   Name-Value Pair Arguments
%     'Parallel'     - Parallel execution mode:
%                      * 'auto' (default) : use PARFOR if DESIREDWORKERS>1,
%                        or if DESIREDWORKERS is empty and number of drops
%                        is large enough.
%                      * 'on'             : force PARFOR (if supported).
%                      * 'off'            : force serial execution.
%     'PoolProfile'  - Parallel pool profile name used by PARPOOL
%                      (default 'Processes').
%     'ReusePool'    - Logical flag to reuse an existing pool when possible
%                      (default true). If false, an existing pool is
%                      deleted and recreated.
%
%   Outputs
%     RESULTS            - Structure with aggregated per-snapshot error and
%                          SNR metrics (concatenated across all drops), for
%                          example:
%                          * timeIndex
%                          * positionErrorX, positionErrorY, positionErrorZ
%                          * rangeError, velocityError
%                          * azimuthError, elevationError
%                          * positionErrorH, positionErrorV
%                          * snr
%     DETSTATS           - Structure with aggregated detection statistics:
%                          * truePositive  (TP)
%                          * falseNegative (FN)
%                          * falsePositve  (FP)
%                          * falseAlarmProb (FPR)
%     DETECTIONOUTPUT    - Struct array concatenating each drop's detection
%                          payload (DROPRESULTS(k).detectionPayload).
%     INFO               - Structure with run-level info:
%                          * K   : number of drops with at least one
%                                  detected object
%                          * cpi : coherent processing interval in PRS
%                                  symbols (derived from precompute)
%
%   Notes
%   - This function falls back to serial execution if Parallel Computing
%     Toolbox or licensing is not available, even when 'Parallel' is 'on'
%     or 'auto' selects parallel.
%
%   See also RUN5GNRAD_PRECOMPUTE, RUN5GNRAD_DROPWORKER, PARFOR, PARPOOL.
%
%   2026 NIST/CTL Steve Blandino
%
%   This file is available under the terms of the NIST License.


p = inputParser;
addParameter(p, 'Parallel', 'auto');        % 'auto' | 'on' | 'off'
addParameter(p, 'PoolProfile', 'Processes');% 'Processes' or 'local'
addParameter(p, 'ReusePool', true);         % reuse existing pool if possible
parse(p, varargin{:});
opt = p.Results;

if nargin < 8 || isempty(desiredWorkers)
    desiredWorkers = [];   % let it auto-pick workers if parallel is used
end

[staticParams, txWaveform] = nrRadar.internal.precompute(simConfig, prsConfig, sensConfig, geometry);
nStDrop = simConfig.nStDrop;
NDrop = floor(size(stConfig.position,1)/simConfig.nStDrop);
NDrop = min(NDrop, simConfig.nMaxDrop);

switch lower(opt.Parallel)
    case 'off'
        useParallel = false;
    case 'on'
        useParallel = true;
    otherwise % 'auto'
        % Only parallelize if user asked for >1 worker,
        % or if they left it empty and NDrop is "big enough".
        useParallel = ( ~isempty(desiredWorkers) && desiredWorkers > 1 ) || ...
                      ( isempty(desiredWorkers) && NDrop >= 4 );
end

% If toolbox/license not available, fall back to serial
if useParallel
    hasPCT = (exist('parpool','file') == 2) && ...
        license('test','Distrib_Computing_Toolbox');
    if ~hasPCT
        useParallel = false;
    end
end

% ---- Pack per-q inputs ----
posCell = cell(NDrop,1);
velCell = cell(NDrop,1);
tgtCell = cell(NDrop,1);

for q = 1:NDrop
    i1 = (q-1)*nStDrop + 1;
    i2 = q*nStDrop;

    posCell{q} = stConfig.position(i1:i2, :);
    velCell{q} = stConfig.velocity(i1:i2, :);

    if isempty(targetChannel)
        tgtCell{q} = [];
    else
        tgtCell{q} = targetChannel(i1:i2);   % slice of struct array
    end
end


    bgCell = cell(NDrop,1);
    for q = 1:NDrop
        if  isempty(backgroundChannel)
            bgCell{q} = [];
        else

            bgCell{q} = backgroundChannel(q);    % struct/array element
        end
    end

% ---- Progress tracking ----
t0     = tic;
pN     = NDrop;
pCount = 0;

% Store results in cell to avoid struct preallocation headaches in parfor
dropResultsCell = cell(NDrop,1);

if useParallel
    % Create queue only if parallel is actually used
    dq = parallel.pool.DataQueue;
    afterEach(dq, @onProgress);

    % Get/reuse pool
    pool = gcp('nocreate');
    if isempty(pool) || ~opt.ReusePool
        % If a pool exists but user doesn't want reuse, delete it
        if ~isempty(pool) && ~opt.ReusePool
            delete(pool);
        end
        if isempty(desiredWorkers)
            pool = parpool(opt.PoolProfile);          % default size
        else
            pool = parpool(opt.PoolProfile, desiredWorkers);
        end
    end

    w = pool.NumWorkers;
    if ~isempty(desiredWorkers)
        w = min(w, desiredWorkers);
    end
    w = min(w, NDrop);

    fprintf('Running PARFOR: %d drops using %d workers (pool has %d)...\n', ...
        NDrop, w, pool.NumWorkers);

    % Ship big stuff once per worker
    preC = parallel.pool.Constant(staticParams);
    txWC = parallel.pool.Constant(txWaveform);

    parfor (q = 1:NDrop, w)
        preLoc = preC.Value;
        txWLoc = txWC.Value;

        dropResultsCell{q} = nrRadar.internal.dropWorker( ...
            q, txWLoc, preLoc, simConfig, sensConfig, ...
            posCell{q}, velCell{q}, bgCell{q}, tgtCell{q});

        send(dq, q);
    end

else
    fprintf('Running FOR (serial): %d drops...\n', NDrop);

    for q = 1:NDrop
        dropResultsCell{q} = nrRadar.internal.dropWorker( ...
            q, txWaveform, staticParams, simConfig, sensConfig, ...
            posCell{q}, velCell{q}, bgCell{q}, tgtCell{q});

        onProgress(q);
    end
end

dropResults = vertcat(dropResultsCell{:});

results.timeIndex = vertcat(dropResults.timeIndex);
results.positionErrorX = vertcat(dropResults.positionErrorX);
results.positionErrorY = vertcat(dropResults.positionErrorY);
results.positionErrorZ = vertcat(dropResults.positionErrorZ);

results.rangeError = vertcat(dropResults.rangeError);
results.velocityError = vertcat(dropResults.velocityError);
results.azimuthError = vertcat(dropResults.azimuthError);
results.elevationError = vertcat(dropResults.elevationError);
results.positionErrorH = vertcat(dropResults.positionErrorH);
results.positionErrorV = vertcat(dropResults.positionErrorV);
detStats.truePositive = [dropResults.TP].';
detStats.falseNegative = [dropResults.FN].';
detStats.falsePositve = [dropResults.FP].';
detStats.falseAlarmProb = [dropResults.FPR].';
results.snr = vertcat(dropResults.snrSlice);

detectionOutput = vertcat(dropResults.detectionPayload);
info.K = sum(detStats.truePositive>=1); %number of drops with at least one detected object


info.cpi = staticParams.prsPeriodicity*staticParams.numberSensingSymbols;

    function onProgress(~)
        pCount = pCount + 1;

        if pCount == 1 || mod(pCount, 10) == 0 || pCount == pN
            elapsed = toc(t0);
            eta = elapsed * (pN - pCount) / max(pCount, 1);
            fprintf('[%s] %d/%d (%.0f%%)  | elapsed %.1fs | ETA %.1fs\n', ...
                datestr(now,'HH:MM:SS'), pCount, pN, 100*pCount/pN, elapsed, eta); %#ok<TNOW1,DATST>
        end
    end
end
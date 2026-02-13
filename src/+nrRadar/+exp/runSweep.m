function [summary, runOut] = runSweep(scenarioPath, runs, varargin)
% RUNSWEEP Run a parameter sweep without creating multiple scenario folders.
%
%   [SUMMARY, RUNOUT] = RUNSWEEP(SCENARIOPATH, RUNS) executes each run
%   definition in RUNS. Each element in RUNS must contain:
%     - tag : string/name for this run
%     - set : cell array of overrides (N-by-2 of {PATH, VALUE})
%
%   The sweep engine supports an optimization for detection-only sweeps:
%   if runs differ only by detection parameters in sensConfig (e.g., CFAR threshold),
%   the range-Doppler cube is computed once per drop and reused across the run list.
%
%   NAME-VALUE OPTIONS
%     'Parallel'        - 'auto' | 'on' | 'off' (default: 'auto')
%     'PoolProfile'     - Parallel pool profile name (default: 'Processes')
%     'ReusePool'       - true/false reuse existing pool (default: true)
%     'desiredWorkers'  - requested workers (default: [])
%     'SaveRoot'        - folder for outputs (default: '') (no export)
%     'Overwrite'       - overwrite existing run folders (default: true)
%
%   Outputs
%     SUMMARY - table with one row per run (basic detection stats)
%     RUNOUT  - struct array with fields results, detStats, detectionOutput, info
%
%   Example (CFAR threshold sweep in dB)
%     thr_dB = 0:30;
%     thr = 10.^(thr_dB/10);
%     runs = nrRadar.exp.makeRuns1D('sens.cfarThreshold', thr, 'TagPrefix','cfar');
%     [summary, runOut] = nrRadar.exp.runSweep('examples/uma_trp1_3gpp', runs, ...
%         'Parallel','on','desiredWorkers',20, 'SaveRoot','results/cfar');
%
%   2026 NIST/CTL Steve Blandino
%   This file is available under the terms of the NIST License.

p = inputParser;
addParameter(p, 'Parallel', 'auto');
addParameter(p, 'PoolProfile', 'Processes');
addParameter(p, 'ReusePool', true);
addParameter(p, 'desiredWorkers', []);
addParameter(p, 'SaveRoot', '');
addParameter(p, 'Overwrite', true);
parse(p, varargin{:});
opt = p.Results;

% Normalize runs input
if ~isstruct(runs)
    error('runSweep:InvalidRuns', 'RUNS must be a struct array with fields tag and set.');
end
if ~isfield(runs,'tag')
    error('runSweep:InvalidRuns', 'Each run must contain a "tag" field.');
end
if ~isfield(runs,'set')
    % allow alternate field name
    if isfield(runs,'Overrides')
        [runs.set] = runs.Overrides;
    else
        error('runSweep:InvalidRuns', 'Each run must contain a "set" field (N-by-2 overrides).');
    end
end

nRunsTotal = numel(runs);
if nRunsTotal == 0
    summary = table();
    runOut = struct([]);
    return;
end

% Load base scenario once
% [simulation0, target0, prs0, geometry0, sens0, backgroundChannel0, targetChannel0] = ...
%     nrRadar.cfg.configureScenario(scenarioPath);

% Split overrides into:
%   - rdmOverrides : anything that changes channel / OFDM / RDM computation
%   - detOverrides : detection-only tweaks (sensConfig)
detOnlySensFields = { ...
    'doaEstimationMethod', ...
    'cfarGrdCellRange','cfarGrdCellVelocity','cfarGrdCellAzimuth','cfarGrdCellElevation', ...
    'cfarTrnCellRange','cfarTrnCellVelocity','cfarTrnCellAzimuth','cfarTrnCellElevation', ...
    'cfarThreshold', ...
    'nmsMaxPeaks','nmsRadius', ...
    'dbscanMinMaxRatio', ...
    'rdaThreshold' ...
    };

rdmKey = strings(nRunsTotal,1);
rdmOverridesAll = cell(nRunsTotal,1);
detOverridesAll = cell(nRunsTotal,1);

for i = 1:nRunsTotal
    setCell = runs(i).set;
    if isempty(setCell)
        rdmOverridesAll{i} = {};
        detOverridesAll{i} = {};
        rdmKey(i) = "";
        continue;
    end
    if isvector(setCell) && mod(numel(setCell),2)==0
        setCell = reshape(setCell,2,[]).';
    end
    if size(setCell,2) ~= 2
        error('runSweep:InvalidOverrideFormat', 'runs(%d).set must be N-by-2 {PATH,VALUE}.', i);
    end

    rdmSet = {};
    detSet = {};
    for k = 1:size(setCell,1)
        path = setCell{k,1};
        value = setCell{k,2};
        [isDetOnly, normPath] = iIsDetectionOnly(path, detOnlySensFields);
        if isDetOnly
            detSet(end+1, :) = {normPath, value}; %#ok<AGROW>
        else
            rdmSet(end+1, :) = {normPath, value}; %#ok<AGROW>
        end
    end
    rdmOverridesAll{i} = rdmSet;
    detOverridesAll{i} = detSet;
    rdmKey(i) = iOverridesKey(rdmSet);
end

% Group runs by identical rdmOverrides
[keysUnique, ~, keyIdx] = unique(rdmKey);
nGroups = numel(keysUnique);

% Preallocate outputs
runOut = repmat(struct('results',[],'detStats',[],'detectionOutput',[],'info',[],'tag',""), nRunsTotal, 1);

% Loop groups
for g = 1:nGroups
    runIdx = find(keyIdx == g);
    if isempty(runIdx)
        continue;
    end

    % Apply the common RDM-affecting overrides for this group
    rdmOverrides = rdmOverridesAll{runIdx(1)};

    [simulation, target, prs, geometry, sens, backgroundChannel, targetChannel] = ...
        nrRadar.cfg.configureScenario(scenarioPath, 'Overrides', rdmOverrides);

    % Precompute drop-invariant quantities for this group
    [staticParams, txWaveform] = nrRadar.internal.precompute(simulation, prs, sens, geometry);

    nStDrop = simulation.nStDrop;
    NDrop = floor(size(target.position,1)/nStDrop);
    NDrop = min(NDrop, simulation.nMaxDrop);

    % Per-run sensing configs (only detection-only tweaks)
    nRunsGroup = numel(runIdx);
    sensPerRun = cell(nRunsGroup,1);
    for r = 1:nRunsGroup
        idxRun = runIdx(r);
        sensPerRun{r} = sens;
        if ~isempty(detOverridesAll{idxRun})
            [~, ~, ~, ~, sensPerRun{r}] = nrRadar.util.applyOverrides(...
                simulation, target, prs, geometry, sensPerRun{r}, detOverridesAll{idxRun});
        end
    end

    % Slice target and channel per drop (as in nrRadar.run)
    posCell = cell(NDrop,1);
    velCell = cell(NDrop,1);
    tgtCell = cell(NDrop,1);
    bgCell  = cell(NDrop,1);
    for q = 1:NDrop
        targetIdx = (q-1)*nStDrop + (1:nStDrop);
        posCell{q} = target.position(targetIdx,:);
        velCell{q} = target.velocity(targetIdx,:);
        if isempty(targetChannel)
            tgtCell{q} = [];
        else
            tgtCell{q} = targetChannel(targetIdx);
        end
        if isempty(backgroundChannel)
            bgCell{q} = [];
        else
            bgCell{q} = backgroundChannel(q);
        end
    end

    % Decide parallel usage (same behavior as nrRadar.run)
    switch lower(opt.Parallel)
        case 'off'
            useParallel = false;
        case 'on'
            useParallel = true;
        otherwise
            useParallel = isempty(getCurrentTask());
    end

    if useParallel
        desiredWorkers = opt.desiredWorkers;
        if isempty(desiredWorkers)
            desiredWorkers = round(feature('numcores')/3);
        end
        pool = gcp('nocreate');
        if isempty(pool)
            parpool(opt.PoolProfile, desiredWorkers);
        else
            if pool.NumWorkers ~= desiredWorkers
                if opt.ReusePool
                    % keep existing pool; user controls pool size
                else
                    delete(pool);
                    parpool(opt.PoolProfile, desiredWorkers);
                end
            end
        end
    end

    dropResultsCell = cell(NDrop,1);

    if useParallel
        parfor q = 1:NDrop
            cache = nrRadar.internal.dropWorker_rdmOnly(q, txWaveform, staticParams, simulation, ...
                sens, posCell{q}, velCell{q}, bgCell{q}, tgtCell{q}); 

            dropRow = cell(nRunsGroup,1);
            for r = 1:nRunsGroup
                [stPositionEstimate, velParallel, peakRda] = nrRadar.sens.rdmDetection( ...
                    cache.rdmCube, sensPerRun{r}, simulation, cache.receiveArray, ...
                    cache.syncOffset, staticParams, 'angleEstimate', cache.analogAngle); %#ok<PFBNS>

                metrics = nrRadar.util.scoreAssociationsPos(cache.stPosition, stPositionEstimate, ...
                    cache.velocityVectors, velParallel, staticParams.txPos);

                dropRow{r} = nrRadar.internal.buildDropOutput(cache.q, metrics, cache.stSNR, ...
                    cache.stPosition, stPositionEstimate, velParallel, peakRda, simulation.nStDrop);
            end

            dropResultsCell{q} = dropRow;
        end
    else
        for q = 1:NDrop
            cache = nrRadar.internal.dropWorker_rdmOnly(q, txWaveform, staticParams, simulation, ...
                sens, posCell{q}, velCell{q}, bgCell{q}, tgtCell{q});

            dropRow = repmat(struct(), 1, nRunsGroup);
            for r = 1:nRunsGroup
                [stPositionEstimate, velParallel, peakRda] = nrRadar.sens.rdmDetection( ...
                    cache.rdmCube, sensPerRun{r}, simulation, cache.receiveArray, ...
                    cache.syncOffset, staticParams, 'angleEstimate', cache.analogAngle);

                metrics = nrRadar.util.scoreAssociationsPos(cache.stPosition, stPositionEstimate, ...
                    cache.velocityVectors, velParallel, staticParams.txPos);

                dropRow(r) = nrRadar.internal.buildDropOutput(cache.q, metrics, cache.stSNR, ...
                    cache.stPosition, stPositionEstimate, velParallel, peakRda, simulation.nStDrop);
            end

            dropResultsCell{q} = dropRow;
        end
    end

    rows = cellfun(@(c) reshape([c{:}], 1, []), dropResultsCell, 'UniformOutput', false);
    dropMat = vertcat(rows{:});% NDrop x nRunsGroup struct

    % dropMat = vertcat(dropResultsCell{:}); % NDrop x nRunsGroup struct

    % Build outputs per run in group
    for r = 1:nRunsGroup
        idxRun = runIdx(r);
        dropResults = dropMat(:,r);

        results = struct();
        results.rangeError = vertcat(dropResults.rangeError);
        results.velocityError = vertcat(dropResults.velocityError);
        results.azimuthError = vertcat(dropResults.azimuthError);
        results.elevationError = vertcat(dropResults.elevationError);
        results.positionErrorX = vertcat(dropResults.positionErrorX);
        results.positionErrorY = vertcat(dropResults.positionErrorY);
        results.positionErrorZ = vertcat(dropResults.positionErrorZ);
        results.positionErrorH = vertcat(dropResults.positionErrorH);
        results.positionErrorV = vertcat(dropResults.positionErrorV);

        detStats = struct();
        detStats.TP = sum([dropResults.TP]);
        detStats.FN = sum([dropResults.FN]);
        detStats.FP = sum([dropResults.FP]);
        detStats.FPR = mean([dropResults.FPR]);
        detStats.nDetectedTarget = sum([dropResults.nDetectedTarget]);

        detectionOutput = [dropResults.detectionPayload];

        info = struct();
        info.SNR_dB = staticParams.SNR_dB;
        info.snrvar = staticParams.snrvar;
        info.NDrop = NDrop;
        info.nStDrop = simulation.nStDrop;
        info.rdmOverrides = rdmOverrides;
        info.detOverrides = detOverridesAll{idxRun};

        runOut(idxRun).results = results;
        runOut(idxRun).detStats = detStats;
        runOut(idxRun).detectionOutput = detectionOutput;
        runOut(idxRun).info = info;
        runOut(idxRun).tag = runs(idxRun).tag;

        % Optional export
        if ~isempty(opt.SaveRoot)
            outFolder = fullfile(opt.SaveRoot, char(runs(idxRun).tag));
            nrRadar.io.exportRun(outFolder, results, detStats, detectionOutput, ...
                'Overwrite', opt.Overwrite);

            % Save metadata (MAT) to preserve complex arrays in configs if needed
            metaFile = fullfile(outFolder, 'runMeta.mat');
            tag = runs(idxRun).tag;
            detOverrides = detOverridesAll{idxRun}; 
            save(metaFile, 'scenarioPath', 'tag', 'rdmOverrides', 'detOverrides', 'simulation', 'prs', 'sens');
        end
    end
end

% Build summary table (lightweight)
tag = strings(nRunsTotal,1);
TP = zeros(nRunsTotal,1);
FN = zeros(nRunsTotal,1);
FP = zeros(nRunsTotal,1);
FPR = zeros(nRunsTotal,1);
for i = 1:nRunsTotal
    tag(i) = string(runs(i).tag);
    if ~isempty(runOut(i).detStats)
        TP(i) = runOut(i).detStats.TP;
        FN(i) = runOut(i).detStats.FN;
        FP(i) = runOut(i).detStats.FP;
        FPR(i) = runOut(i).detStats.FPR;
    end
end
summary = table(tag, TP, FN, FP, FPR);

end

function [isDetOnly, normPath] = iIsDetectionOnly(path, detOnlySensFields)
% Normalize prefix aliases and decide if this override is detection-only.
normPath = string(path);
normPath = strrep(normPath, 'sensConfig.', 'sens.');
normPath = strrep(normPath, 'simConfig.', 'simulation.');
normPath = strrep(normPath, 'stConfig.', 'target.');
normPath = strrep(normPath, 'prsConfig.', 'prs.');

parts = strsplit(normPath, '.');
if numel(parts) < 2
    isDetOnly = false;
    return;
end

prefix = parts{1};
field1 = parts{2};

isDetOnly = (prefix == "sens") && any(strcmp(field1, detOnlySensFields));
end

function key = iOverridesKey(setCell)
% Build a deterministic string key for a set of overrides.
if isempty(setCell)
    key = "";
    return;
end
if isvector(setCell) && mod(numel(setCell),2)==0
    setCell = reshape(setCell,2,[]).';
end

n = size(setCell,1);
pairs = strings(n,1);
for k = 1:n
    p = string(setCell{k,1});
    v = setCell{k,2};
    if ischar(v) || isstring(v)
        vs = string(v);
    elseif islogical(v)
        vs = string(v);
    else
        try
            vs = string(mat2str(v));
        catch
            vs = "val";
        end
    end
    pairs(k) = p + "=" + vs;
end
pairs = sort(pairs);
key = strjoin(cellstr(pairs), "|");
end

# 5GNRad Test Suite

This folder contains the MATLAB Unit Test suite used to validate backward compatibility as new features are added to **5GNRad**.

## Test layers

### Component 
Deterministic tests around core pipeline building blocks:
- `nrRadar.internal.precompute`
- `nrRadar.rx.getRxWaveform`
- `nrRadar.rx.estimateChannel`
- `nrRadar.sens.getRangeDoppler`
- `nrRadar.sens.rdmDetection`

**Purpose:** ensure key **function contracts** (I/O schema and invariants) do not regress.

### Integration / regression
Scenario-based tests (subset of `examples/`) comparing output summaries against stored baselines under `tests/baselines/`.

Includes an **API compatibility gate** that detects changes to public function signatures (unless intentionally updated).

**Purpose:** catch subtle numerical or behavioral changes that pass component tests but alter end-to-end results.

### Performance 
Performance regression checks (runtime + output-size proxy for memory) against stored baseline budgets.

**Purpose:** detect accidental slowdowns or memory blow-ups without brittle timing checks.


---

## Requirements

Tests require toolboxes:
- 5G Toolbox
- Phased Array System Toolbox


## Running tests

From the repo root:

```matlab
setup;
addpath("tests");

results = runTests;                      % all tests discovered
results = runTests('Tag','layerB');      % only Layer B
results = runTests('Tag','compat');      % API compatibility gate (Layer C)
results = runTests('Tag','regression');  % Layer C baseline comparisons
results = runTests('Tag','performance'); % Layer D performance budget

assert(all([results.Passed]));
```

---

## Generating/updating baselines

### Layer C baselines (scenario regression)

```matlab
setup;
addpath("tests");
makeBaselines;
```

This writes baseline files under `tests/baselines/`.

### Public API signature baseline (compatibility gate)

```matlab
setup;
addpath("tests");
generatePublicApiBaseline;
```

This writes `tests/compat/publicApiSignatures.json`.



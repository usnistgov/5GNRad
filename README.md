# 5GNRad: 5G New Radio Radar

> MATLAB end to end radar processing chain for 5G NR.

5GNRad is a MATLAB-based simulation platform to evaluate **sensing with cellular infrastructure** using **standard-compliant 5G NR Positioning Reference Signals (PRS)** and a radar-style processing chain that produces **range–Doppler–angle detections** and **target position/velocity estimates**.




---

## Key features

- **Standard-compliant PRS waveform generation** (NR numerology + OFDM mapping)
- **Two operating modes**
  1. **Uploaded MPC channels (3GPP Release 19, TR 38.901 sensing extensions)** via `backgroundChannel.json` / `targetChannel.json`
  2. **Simplified single-point target** (single scatterer) with **RCS-based amplitude scaling** aligned with TR 38.901 sensing formulations
- **Radar processing chain**
  - Range processing (IFFT / delay profile)
  - Doppler processing (slow-time FFT)
  - Angle estimation (beamspace FFT or Bartlett scan, depending on config)
  - Detection: 2D CFAR + peak picking + clustering 
- **Scenario-driven**: each example is self-contained under `examples/<scenarioName>/Input/`
- **Parallel execution** over drops (`parfor`) when enabled


## Requirements

Validated on **MATLAB R2025a**.

### Required MathWorks products
- **MATLAB**
- **5G Toolbox** (NR numerology, PRS generation, OFDM grid mapping)
- **Phased Array System Toolbox** (arrays, steering vectors, angle utilities)
- **Signal Processing Toolbox** (windowing and FFT helpers)

### Optional
- **Parallel Computing Toolbox** (parallel execution across drops)

---

## Repository structure

```
<repo-root>/
  main.m
  setup.m
  src/
    +nrRadar/
      +cfg/        (scenario configuration readers/validators)
      +rx/         (RX synthesis, OFDM demod, RD/RDA formation)
      +sens/       (CFAR, NMS, DOA, clustering, sidelobe suppression, geometry)
      +array/      (array utilities and coordinate conventions)
      +io/         (export utilities)
      +internal/   (precompute + per-drop worker entry points)
  examples/
    <scenarioName>/
      Input/       (configs + optional channel JSONs)
      Output/      (generated results)
```

---

## Quick start

From the repository root:

```matlab
setup();
main
```

`main.m` runs a batch of scenarios listed inside the script (config → precompute → per-drop processing → export).


## Run a single scenario

```matlab
setup();

scenarioPath = "examples/uma_trp1_3gpp";

[simConfig, stConfig, prsConfig, geometry, sensConfig, backgroundChannel, targetChannel] = ...
    nrRadar.cfg.configureScenario(scenarioPath);

desiredWorkers = 8;          % optional
parallelMode   = "off";      % "on" | "off" 

[results, detStats, detectionOutput] = nrRadar.run( ...
    simConfig, stConfig, prsConfig, geometry, sensConfig, ...
    backgroundChannel, targetChannel, desiredWorkers, ...
    "Parallel", parallelMode);

nrRadar.io.exportResults(scenarioPath, results, detStats, detectionOutput);
```

---

## Scenario inputs

Each scenario under `examples/<scenarioName>/Input/`  includes:

- `simulationConfig.txt` — carrier, bandwidth, noise figure, scenario options, etc.
- `prsConfig.txt` — PRS resource definition (periodicity, comb size, symbol mapping, etc.)
- `sensingConfig.txt` — processing and detection settings (FFT sizes, CFAR, DOA method, thresholds, clustering)
- `txAntennaModel.json`, `rxAntennaModel.json` — array geometry + beamforming/combining configuration
- `backgroundChannel.json` *(optional)* — environment/clutter MPCs (uploaded-channel mode)
- `targetChannel.json` *(optional)* — target MPCs (uploaded-channel mode)
- `bsConfiguration.txt`  — TRP/BS configuration
- `targetConfiguration.txt`  — target truth states (position/velocity timeline)


### Operating mode selection
- If `targetChannel.json` / `backgroundChannel.json` exist, the simulator uses **uploaded MPC mode**.
- If target MPCs are omitted (or the scenario is configured accordingly), the simulator can run in **simplified single-point target mode** (single scatterer with RCS-based amplitude scaling).



## Outputs

After a run, results are written under:

```
examples/<scenarioName>/Output/
```

Outputs include:

- `detection.json` — per-drop detection payloads (estimated positions/velocities and related metrics)
- `detStats.csv` — per-drop detection statistics (TP / FN / FP / false-alarm probability)
- `error.csv` — per-target/per-drop error table (position/range/angle/velocity errors; NaN if unassigned)

---

## Performance tips

- Use  `"Parallel","on"` with a sensible `desiredWorkers` to speed up batch runs.
- Runtime is dominated by:
  - number of drops / snapshots
  - range/Doppler/angle FFT sizes
  - number of RF chains (digital channels)
  - MPC count and channel length
- If you hit out-of-memory errors, reduce worker count.


---

## Documentation

Documentation is available in the [docs/5G_NR_Radar_doc.pdf](docs/5G_NR_Radar_doc.pdf) file.

---

## Publications and 3GPP Contributions

- S. Blandino et al., “[Detecting Airborne Objects with 5G NR Radars](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=11310754),” IEEE MILCOM 2025.

- National Institute of Standards and Technology (NIST), “[R1-2505684 Discussion on ISAC Performance Evaluation](https://www.3gpp.org/ftp/tsg_ran/WG1_RL1/TSGR1_122/Docs/R1-2505684.zip),” 3GPP TSG-RAN WG1 Meeting #122

- National Institute of Standards and Technology (NIST), “[R1-2508883 Discussion on ISAC Performance Evaluation](https://www.3gpp.org/ftp/tsg_ran/WG1_RL1/TSGR1_123/Docs/R1-2508883.zip),” 3GPP TSG-RAN WG1 Meeting #123

- National Institute of Standards and Technology (NIST), “[R1-2600846 Discussion on ISAC Performance Evaluation](https://www.3gpp.org/ftp/meetings_3gpp_sync/RAN1/Docs/R1-2600846.zip),” 3GPP TSG-RAN WG1 Meeting #124



---

## Contributing

Feedback and contributions are welcome! Please open an issue or contact the maintainer.



## Contact

Steve Blandino  
NIST Communications Technology Laboratory  
[steve.blandino@nist.gov](mailto:steve.blandino@nist.gov)  
[https://www.nist.gov/people/steve-blandino](https://www.nist.gov/people/steve-blandino)


## Citing 5GNRad

If you use **5GNRad** to generate figures or results for a publication or standards contribution, please cite the associated manuscript and the software repository.

**Manuscript**
- S. Blandino et al., "Detecting Airborne Objects with 5G NR Radars," MILCOM 2025 - 2025 IEEE Military Communications Conference (MILCOM), Los Angeles, CA, USA, 2025, pp. 1260-1265, doi: 10.1109/MILCOM64451.2025.11310754.

**Software**
- National Institute of Standards and Technology (NIST) Communications Technology Laboratory, *5G New Radio Radar (5GNRad)* [Computer software]. Available: https://github.com/usnistgov/5GNRad

<details>
<summary>BibTeX</summary>

```bibtex
@INPROCEEDINGS{11310754,
  author={Blandino, Steve and Golmie, Nada and Sahoo, Anirudha and Nguyen, Thao and Ropitault, Tanguy and Griffith, David and Sonny, Amala},
  booktitle={MILCOM 2025 - 2025 IEEE Military Communications Conference (MILCOM)}, 
  title={Detecting Airborne Objects with 5G NR Radars}, 
  year={2025},
  volume={},
  number={},
  pages={1260-1265},
  keywords={Three-dimensional displays;Radar clutter;Simulation;Airborne radar;Urban areas;Radar detection;Estimation;Autonomous aerial vehicles;Integrated sensing and communication;Clutter;ISAC;5G mobile communication;6G mobile communication;3GPP Standards;Target Detection},
  doi={10.1109/MILCOM64451.2025.11310754}}

@software{nist_5gnrad,
  author  = {{National Institute of Standards and Technology (NIST) Communications Technology Laboratory}},
  title   = {5G New Radio Radar (5GNRad)},
  url     = {https://github.com/usnistgov/5GNRad},
  note    = {Computer software}
}
# 5GNRad: 5G New Radio Radar

> MATLAB end to end radar processing chain for 5G NR.

**5GNRad** is a MATLAB simulation framework modeling the radar sensing capability of a monostatic 5G NR Base Station (BS), reusing the Positioning Reference Signal (PRS) defined in the 3GPP TS 38.211 standard to estimate the delay and angle of airborne targets. The tool supports propagation scenarios compliant with TR 38.901 and includes recent Release 19 RCS and clutter models.

---

## Features

- End-to-end radar processing: PRS generation, OFDM transmission, beamforming, clutter suppression, range-Doppler estimation.
- 3GPP-compliant channel models: UMi, UMa, and RMa with LoS/NLoS and multipath.
- Realistic sensing scenarios with configurable target trajectories.

---

## Folder Structure

```
./
├─ main.m                        # Entry point of the simulation
├─ examples/                    # Example scenarios
│   └─ <scenarioName>/
│      └─ Input/                # Scenario-specific configuration files
│          ├─ bsConfig.txt
│          ├─ prsConfig.txt
│          ├─ sensConfig.txt
│          ├─ simulationConfig.txt
│          └─ targetConfig.txt
├─ channel/          # Pre-generated channel files (JSON)
│   ├─ backgroundChannel_6GHz_3GPP_38.901_UMa_NLOS.json
│   └─ targetChannel_6GHz_UMaAV.json
├─ src/                        # Source code for PRS, sensing, and simulation
└─ README.md                   # This file
```

---

## Installation

Simply clone or download the repository. No extra installation is required.

---

## Requirements

- MATLAB R2025a (tested)
- [MATLAB 5G Toolbox](https://www.mathworks.com/products/5g.html)

---

## How to Run

1. Open `main.m` in MATLAB.
2. Set the scenario name:
   ```matlab
   scenarioNameStr = "UMi-Av25";
   ```
3. Run the script.

You can modify input parameters via the configuration files located in `examples/<scenario>/Input/`.

---

## Documentation

Documentation is available in the [docs/5G_NR_Radar_doc.pdf](docs/5G_NR_Radar_doc.pdf) file.

---

## Reference

- S. Blandino et al., “[Detecting Airborne Objects with 5G NR Radars](https://arxiv.org/pdf/2505.24763),” IEEE MILCOM 2025.

- National Institute of Standards and Technology (NIST), “[R1-2505684 Discussion on ISAC Performance Evaluation](https://www.3gpp.org/ftp/tsg_ran/WG1_RL1/TSGR1_122/Docs/R1-2505684.zip),” 3GPP TSG-RAN WG1 Meeting #122

---

## Contributing

Feedback and contributions are welcome! Please open an issue or contact the maintainer.

---

## Contact

Steve Blandino  
NIST Communications Technology Laboratory  
[steve.blandino@nist.gov](mailto:steve.blandino@nist.gov)  
[https://www.nist.gov/people/steve-blandino](https://www.nist.gov/people/steve-blandino)
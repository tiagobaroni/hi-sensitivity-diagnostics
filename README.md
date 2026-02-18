# HI Sensitivity Diagnostics

[![License: GPL-3.0](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)
[![Status](https://img.shields.io/badge/status-beta-orange.svg)]()
[![Ko-Fi](https://img.shields.io/badge/Support-Buy%20Me%20A%20Coffee-ff5f5f?style=flat&logo=ko-fi)](https://ko-fi.com/tiagobaroni)

**A physically-motivated diagnostic tool for integration-time planning and sensitivity validation in 21cm HI radio observations.**

This project implements an audit-oriented methodology to distinguish between **thermal-noise-limited performance** and **instrumental/systematic floor limitations** in spectral HI observations.

The goal is not absolute calibration, but operational sensitivity diagnostics and regime classification.

Instead of blindly extrapolating integration time via the Radiometer Equation, this tool evaluates whether additional integration is physically justified.

---

## Scientific Rationale

This tool is strictly based on spectral (line-free) RMS diagnostics. Total power metrics are not used for sensitivity scaling decisions.

The classical Radiometer Equation predicts:

σ ∝ 1 / √τ

However, real-world radio observations are often limited by systematic effects such as:

- Gain drift  
- Baseline ripple (standing waves)  
- RFI contamination  
- Spectral smoothing correlations  
- Instrumental instability  

This tool formalizes the decomposition:

σ_obs² = A/τ + σ_floor²

Where:

- A/τ → Thermal component  
- σ_floor → Instrumental/systematic floor  

The system estimates both terms and classifies the dominant noise regime.

---

## Core Capabilities

### 1. Spectral RMS vs Total Power Separation

The system explicitly distinguishes between:

- `spectral_linefree_rms` → Used for integration decisions  
- `total_power_metric` → Diagnostic only  

```text
# WARNING:
# total_power is diagnostic only.
# Integration-time decisions must rely on spectral_linefree_rms.
```
This prevents false scaling decisions caused by drift.

### 2. Thermal vs Instrumental Regime Classification
- Computes Rσ = σ_obs / σ_th
- Estimates slope of log(σ) vs log(τ)
- Fits model (assuming approximate stationarity of the instrumental floor):

σ² = A/τ + σ_floor²

Outputs:
- Estimated instrumental floor;
- Integration efficiency;
- Regime classification;

### 3. Drift Detection
Drift detection via:
- Mean power per subintegration;
- Low-order polynomial temporal fit;
- Drift flagging;

If drift is detected:
- Thermal scaling validity is flagged;
- Baseline or gain-stability correction is recommended before further integration decisions.

### 4. Allan Variance
Implements:
``` python
 compute_allan_variance(power_series)
```

Used to determine:
- Optimal subintegration time 𝜏𝑠
- Stability regime transitions
- Whether increasing 𝜏𝑠 is worse than increasing N

### 5. Cold Sky Reference Validation
Allows comparison against:
- Known cold-sky reference positions
- External survey spectra (e.g., LAB survey)
- Beam-convolved reference RMS

Computes **reference deviation (%)** and configurable threshold (default: 20%).

---

## Methodological Assumptions

- Spectral RMS is computed after baseline removal and masking of line-free regions.
- Instrumental floor is assumed quasi-stationary over the integration window.
- Correlated spectral smoothing effects are implicitly absorbed in Δν_eff.

---

## Project Structure
```text
hi-sensitivity-diagnostics/
│
├── hi-sensitivity-diagnostics.py   # Main CLI diagnostic tool
├── options.py                      # Runtime configuration container
├── gh_processing_utils.py          # Noise modeling & statistical utilities
│
├── docs/
│   └── Sensitivity_Plan_example.pdf
│
├── pyproject.toml
├── requirements.txt
└── README.md
```

---

## Installation
Requirements:
- Python ≥ 3.10
- numpy
- scipy
- astropy
- matplotlib
- fpdf

Clone and install:
```
git clone https://github.com/tiagobaroni/hi-sensitivity-diagnostics.git
cd hi-sensitivity-diagnostics
pip install .
```

Or install directly:
```
pip install git+https://github.com/tiagobaroni/hi-sensitivity-diagnostics.git
```

---

## Usage

Basic example:
``` python
python hi-sensitivity-diagnostics.py 20 0.8 20 --n_max 100 --steps 10
```

With Allan variance:
``` python
python hi-sensitivity-diagnostics.py 20 0.8 20 --allan
```

With cold sky reference:
``` python
python hi-sensitivity-diagnostics.py 20 0.8 20 --cold_sky reference.fits
```

## Output
The tool generates:
- Regime classification report
- Estimated instrumental floor
- Thermal scaling projection
- Allan variance plot (optional)
- Cold sky deviation report (optional)
- PDF diagnostic summary

---

## License
This project is distributed under the GNU General Public License v3.0 (GPL-3.0). See the LICENSE file for details.


<div align="center">
  <p>If this tool helps your research or learning, consider supporting the development!</p>
  <a href="https://ko-fi.com/tiagobaroni" target="_blank">
    <img
      src="https://storage.ko-fi.com/cdn/kofi2.png"
      height="36"
      style="border:0px;height:36px;"
      alt="Buy Me a Coffee at ko-fi.com"
    />
  </a>
</div>
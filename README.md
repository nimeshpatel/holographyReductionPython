# GLT Near-Field Holography Reduction Pipeline
(This README file was created by Claude, based on Section 5 of our SPIE holography paper).

A Python pipeline for reducing near-field holography data from the
[Greenland Telescope (GLT)](https://www.cfa.harvard.edu/facilities-technology/radio-telescopes/greenland-telescope)
at Pituffik Space Base, Greenland.

This package is a Python reimplementation of the
[SMA HOLIS](https://iweb.cfa.harvard.edu/sma/memos/116.pdf) Fortran/C holography
reduction package, extended with boresight calibration, 2-D phase unwrapping,
interactive GUIs, and automated pipeline scripting.

Reference: Patel et al. (2026), *Near-field holography system for the Greenland
Telescope*, SPIE Astronomical Telescopes and Instrumentation, paper 14153-24.

---

## Table of Contents

- [Overview](#overview)
- [System Description](#system-description)
- [Pipeline Steps](#pipeline-steps)
- [Installation](#installation)
- [Usage](#usage)
- [Parameter Files](#parameter-files)
- [Output Files](#output-files)
- [Map Sizes](#map-sizes)
- [Mathematical Background](#mathematical-background)
- [Known Issues and Corrections](#known-issues-and-corrections)
- [Authors](#authors)
- [References](#references)

---

## Overview

The GLT holography system measures the complex beam pattern of the 12-m antenna
by scanning a 94.5 GHz beacon transmitter located on South Mountain ridge,
approximately 2.34 km from the telescope at an elevation of ~2.9°. The 86 GHz
science receiver and a reference receiver mounted behind the subreflector
simultaneously measure amplitude and phase via a vector voltmeter (VVM) and ADC.
The measured far-field beam pattern is Fourier transformed to recover the aperture
illumination and surface error map of the primary reflector.

**Key specifications:**

| Parameter | Value |
|-----------|-------|
| Frequency | 94.5 GHz (λ = 3.172 mm) |
| Dish diameter | 12 m |
| Beacon distance | 2.34 km |
| Beacon elevation | ~2.9° |
| Map grid | 32×32, 64×64, 96×96, 128×128, or 160×160 pixels |
| Pixel spacing | 41 arcsec |
| Aperture resolution (128×128) | 14 cm |
| Typical surface RMS | ~35–40 µm |

---

## System Description

```
Transmitter (94.5 GHz, 10 dBm, ~6° FWHM beam)
        |
        | 2.34 km
        ↓
 GLT 12-m Dish
        |
   ┌────┴────┐
   │         │
86 GHz    Reference receiver
receiver   (behind subreflector)
   │         │
   └────┬────┘
        │ 21.4 MHz IF (after down-conversion)
        ↓
   Vector Voltmeter (HP 8508A)
        │
        ↓ Amplitude + Phase
   ADC + GPS/IRIG-B timestamping
        │
        ↓
   [Time, Az, El, Amplitude, Phase]
```

The antenna scans in on-the-fly mode: continuous azimuth sweeps at constant
elevation, stepping in elevation between rows. Encoder readings (Az, El) and VVM
readings (amplitude, phase) are time-stamped with IRIG-B from GPS and collated
into a time-series data file.

---

## Pipeline Steps

The pipeline is implemented as a bash script (`holoMap.sh`) that chains the
following Python programs. A flowchart is provided in Fig. 4 of the reference paper.

### Step 1 — `detect_raster_start.py`
Automatically detects the start of the raster scan by identifying a sustained
positive azimuth slope combined with flat elevation. Discards pre-scan transient
data and writes a trimmed output file with timestamps.

**Input:** Raw ADC data file (`holoADC-*.txt`)
**Output:** `trimmed_with_time.txt`

### Step 2 — `boresight_cal.py` *(optional)*
Applies boresight phase-drift calibration. The antenna briefly returns to the
beacon position (~3 s) between each azimuth row. A cubic spline is fitted to
the sequence of boresight phase measurements and subtracted from each raster row
to correct for slow instrumental phase drifts.

**Input:** `trimmed_with_time.txt`
**Output:** `calibrated.txt`, `*_boresight_data`

### Step 3 — `regrid_holo.py`
Regrids the irregularly sampled on-the-fly data onto a regular N×N grid
(32×32, 64×64, or 128×128 pixels, each 41 arcsec). Within each grid cell,
amplitude and phase are combined using **complex vector averaging**:

```
z̄ = (1/N) Σ aₙ exp(iφₙ)
A = |z̄|,   Φ = arg(z̄)
```

This correctly handles phase wrapping during averaging.

**Input:** `trimmed_with_time.txt` (or `calibrated.txt`), `regrid_NxN.prm`
**Output:** `rgin.dat` (j k amplitude phase)

### Step 4 (Optional- run this only if missing data are found in overview plot) — `fix_missing_cell.py`
Interpolates any empty grid cells left by boresight gaps or missing data.

### Step 5 — `preprocess.py`
Applies three phase corrections to the beam-map data before the Fourier
transform:

1. **Near-field (Fresnel) correction** — corrects for spherical wavefront
   curvature from the finite-range beacon:
   `φ_NF = π r² / (λ R)`

2. **Reference-plane correction** — corrects for the path-length difference
   between the signal and reference receiver planes relative to the antenna
   pivot point:
   `Δφ_ref = (2π/λ)(d₁ − d₂)(1 − cos θ)`

3. **Boresight phase drift correction** (mode 1) — spline-interpolated phase
   drift from interleaved boresight measurements is subtracted row by row.

**Input:** `rgin.dat`, `preprocess_N.prm`
**Output:** `ampout.dat`, `phaseout.dat`

### Step 5b — `plot_beam_maps.py`
Plots the pre-FFT beam amplitude and phase maps for inspection.

**Output:** `*_beam_maps.pdf`

### Step 6 — `holis_aber2.py` *(first pass)*
Fourier transforms the corrected far-field maps to the aperture plane and fits
and removes large-scale optical aberrations from the **wrapped** aperture phase:

```
E_aperture(x,y) = F[E_far-field(u,v)]
```

After the FFT, an eight-term aberration basis is fitted and subtracted:
```
Φ_aber = a₀ + a₁x + a₂y + a₃φ_defocus + a₄(2xy) + a₅(x²−y²) + a₆(x³−3xy²) + a₇(3x²y−y³)
```
comprising DC offset, x/y tilt, defocus (Ruze formula), 45° and 0° astigmatism,
and x/y coma. The fitted coefficients are saved to `holis_fit.json`.

**Input:** `ampout.dat`, `phaseout.dat`, `withphase_aber.prm`
**Output:** `Ep.dat` (wrapped residual aperture phase), `holis_fit.json`

### Step 7 — `unwrap2d.py`
Performs 2-D phase unwrapping on the aperture phase map to resolve 2π
ambiguities and produce a continuous, single-valued phase surface. This is
essential when large-scale phase gradients (from pointing offsets, defocus, etc.)
cause the aperture phase to wrap over 2π.

**Input:** `Ep.dat`
**Output:** `tk.dat` (unwrapped aperture phase)

### Step 8 — `holis_aber2.py --unwrap` *(second pass)*
With the continuous unwrapped phase available, performs a refined Fourier
transform and eight-term aberration fit on the unambiguous phase surface.
The residual phase after aberration removal is converted to surface displacement:

```
ε = (λ/4π) Φ_residual
```

with a geometric projection from half-path-length error to normal surface
displacement: `ε_⊥ = ε / cos α`, where `cos α = (1 + r²/4f²)^{-1/2}`.

**Input:** `tk.dat`, `withphase_aber.prm`
**Output:** `Epr.dat` (residual surface error, radians), `Ea_um.dat` (aperture amplitude)

### Step 9 — Save results
Copies all output files to `results/` with a dataset-specific filename prefix
derived from the observation timestamp.

### Steps 10a/10b — `glt_dish_map.py`
Renders the aperture illumination amplitude map and surface error map as PDF
figures with panel boundary overlays. Converts phase to surface displacement
in microns and applies the cos α geometric projection.

**Output:** `*_illumination.pdf`, `*.pdf` (surface error map), `*_rms.json`

### Step 11 — `create_summary_pdf.py`
Combines the illumination and surface error maps side-by-side with dataset
metadata (filename, surface RMS, fitted aberration coefficients table) into a
single summary PDF.

**Output:** `*_summary.pdf`

---

## Installation

### Requirements

```
Python >= 3.8
numpy
scipy
matplotlib
pandas
```

### Install dependencies

```bash
pip install numpy scipy matplotlib pandas
```

### Clone the repository

```bash
git clone git@github.com:nimeshpatel/holographyReductionPython.git
cd holographyReductionPython
```

---

## Usage

### Full pipeline (128×128 map, no boresight calibration)

```bash
bash holoMap.sh holoADC-AzEl-20260519_123929.txt holoADC-20260519 0
```

Arguments:
1. Input raw data file
2. Output prefix (used for results filenames)
3. Boresight mode (0 = none, 1 = spline correction, 2/3 = other modes)

### Individual steps

```bash
# Detect raster start
python detect_raster_start.py holoADC-AzEl-20260519_123929.txt -o trimmed.txt

# Regrid to 128×128
python regrid_holo.py trimmed.txt regrid_128x128.prm

# Preprocess
python preprocess.py   # reads preprocess.prm

# First FFT pass
python holis_aber2.py

# 2D phase unwrapping
python unwrap2d.py

# Second FFT pass with unwrapped phase
python holis_aber2.py --unwrap

# Plot surface error map
python glt_dish_map.py results/PREFIX_Epr.dat --vmin -150 --vmax 150
```

### Average multiple maps

```bash
python average_maps.py results/map1_Epr.dat results/map2_Epr.dat ...
bash averageMaps.sh
```

### Interactive GUI

```bash
python holoMapGUIv2.py
```

---

## Parameter Files

### `regrid_128x128.prm`
Single line: `naz nel azstrt azstep0 elstrt elstep az_thresh el_thresh`

Example:
```
128 128 230.909 0.011400 2.143400 0.011400 0.006000 0.006000
```

### `preprocess_128.prm`
Fixed-format parameter file specifying:
- Input/output filenames
- Grid size (ninp, nout)
- Observing frequency (GHz)
- Sampling interval (arcsec)
- Reference plane distances d₁, d₂ (m)
- VVM/DC/ADC gain and offset values
- Input units (engineering or counts)

### `withphase_aber.prm`
Fixed-format parameter file for `holis_aber2.py` specifying:
- Input/output filenames
- Grid dimension, Nyquist sampling rate
- Beacon distance (m)
- Frequency (GHz), sampling interval (arcsec)
- Primary diameter, secondary diameter
- Primary focal length, Cassegrain magnification
- Near-field correction flag
- Defocus correction (mm)
- Masking parameters (outer/inner diameter, quadrupod half-width)
- Aberration fitting flags (8 terms)

Symbolic links `preprocess.prm`, `regrid.prm`, and `withphase_aber.prm` point
to the active parameter files for the current map size.

---

## Output Files

| File | Description |
|------|-------------|
| `rgin.dat` | Regridded beam map (j k amplitude phase) |
| `ampout.dat` | Preprocessed far-field amplitude |
| `phaseout.dat` | Preprocessed far-field phase |
| `Ep.dat` | Wrapped aperture phase after first FFT pass |
| `holis_fit.json` | Fitted aberration coefficients and RMS values |
| `tk.dat` | Unwrapped aperture phase |
| `Epr.dat` | Residual surface error map (radians) |
| `Ea_um.dat` | Aperture illumination amplitude (unmasked) |
| `mask32.dat` | Aperture mask (1=valid, 0=masked) |
| `results/PREFIX_Epr.dat` | Archived surface error map |
| `results/PREFIX_summary.pdf` | Summary page with maps and fit table |
| `results/PREFIX.log` | Full processing log |

---

## Map Sizes

| Grid | Map range (±deg) | Row duration (s) | Map duration | Resolution |
|------|-----------------|-----------------|--------------|------------|
| 32×32 | ±0.182 | 4.9 | ~10 min | 55.9 cm |
| 64×64 | ±0.364 | 9.7 | ~30 min | 27.9 cm |
| 96×96 | ±0.547 | 14.6 | ~45 min | 18.6 cm |
| 128×128 | ±0.729 | 19.4 | ~60 min | 14.0 cm |
| 160×160 | ±0.911 | 24.3 | ~90 min | 11.2 cm |

Pixel spacing is fixed at 41 arcsec; scan rate is 270 arcsec/s.
Aperture-plane resolution = D / (N × s / θ_FWHM), where θ_FWHM = 61.1 arcsec
at 94.5 GHz and D = 12 m.

---

## Mathematical Background

### Fourier transform relationship
```
E_aperture(x,y) = F[ E_far-field(u,v) ]
```

### Near-field (Fresnel) correction
```
φ_NF = π r² / (λ R)
```
where r is the aperture radial coordinate, λ = 3.172 mm, R = 2.34 km.

### Reference-plane correction
```
Δφ_ref = (2π/λ)(d₁ − d₂)(1 − cos θ)
```
where d₁ = 2.180 m (main dish reference plane),
d₂ = 7.727 m (reference horn plane), θ = off-axis angle.

### Defocus correction (Ruze formula)
```
φ_defocus = (4π Δf / λ) [ (r/2fₚ)² / (1+(r/2fₚ)²) + (r/2f_m)² / (1+(r/2f_m)²) ]
```
where fₚ = 4.8 m (primary focal length), f_m = fₚ × magnification.

### Surface error from phase
```
ε = (λ/4π) Φ_residual
ε_⊥ = ε / cos α,   cos α = (1 + r²/4f²)^{-1/2}
```

---

## Known Issues and Corrections

### Bug fix: Complex vector averaging in regrid_holo.py
The original C code (`regrid3.c`) computed amplitude from the **last** sample
in each grid cell rather than the complex average. This has been corrected:

```python
# CORRECTED: use complex average
amp = sqrt(x_avg² + y_avg²)
phase = atan2(y_avg, x_avg)
```

### Bug fix: Scale factor sign in preprocess.py
The original parameter parsing used `abs()` when computing the phase scale
factor, suppressing the sign of the VVM phase gain. Corrected to preserve sign:

```python
scale_ph = 1.0 / (vvm_g_ph * dc_g_ph * adc_g_ph)  # sign preserved
```

### Why two passes of holis_aber2.py?
Large-scale aberrations (pointing tilts, defocus) create phase gradients across
the aperture that wrap over 2π. Fitting smooth polynomial functions to wrapped
phase data fails. The two-pass approach:
1. First pass: FFT → get wrapped aperture phase → first-order aberration removal
2. `unwrap2d.py`: resolve 2π ambiguities → continuous phase surface
3. Second pass: refit aberrations on unambiguous phase → accurate residuals

---

## Authors

Nimesh A. Patel (CfA), Satoki Matsushita (ASIAA), Derek Kubo (ASIAA),
Tirupati K. Sridharan (NRAO), and the GLT team.

Center for Astrophysics | Harvard & Smithsonian, Cambridge MA, USA
Academia Sinica Institute of Astronomy and Astrophysics, Taipei, Taiwan

Contact: npatel@cfa.harvard.edu

---

## References

1. Chen, M.-T. et al., "The Greenland Telescope," PASP 135, 095001 (2023).
2. Baars, J. W. M. et al., "Near-Field Radio Holography of Large Reflector
   Antennas," IEEE Antennas and Propagation Magazine 49, 24–41 (2007).
3. Zhang, X., "HOLIS," SMA Technical Memorandum No. 116 (1996).
   https://iweb.cfa.harvard.edu/sma/memos/116.pdf
4. Sridharan, T. K., "Holography Analysis Software–HOLISv2,"
   SMA Technical Memorandum No. 169 (2020).
   https://iweb.cfa.harvard.edu/sma/memos/169.pdf
5. Patel, N. A. et al., "Near-field holography system for the Greenland
   Telescope," SPIE Astronomical Telescopes and Instrumentation (2026).

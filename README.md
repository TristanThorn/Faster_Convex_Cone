# CC Unmixing

This repository contains MATLAB scripts for photoacoustic (PA) and ultrasound (US) imaging workflows. It includes routines for frequency–domain reconstruction, digital phantom creation and spectral analysis—and now incorporates a **convex‐cone–constrained unmixing** module to correct depth- and background-dependent SO₂ bias observed with simple linear unmixing.

---

## MATLAB Requirements

- **MATLAB R2018a** or later  
- **Image Processing Toolbox** – for `imshow`, `volumeViewer`, `mean2`, `std2`  
- **Statistics and Machine Learning Toolbox** (optional) – for k-means clustering in phantom pruning  

```matlab
addpath(genpath('path/to/cc_unmixing'));
```

---

## Directory Overview

* `Toolbox/`

  * Spectral‐unmixing helpers
  * Convex‐cone routines (`cone_unmix.m`, `build_pa_cones.m`)
  * Geometry & phantom utilities
* `Substance_spectra/`

  * `.mat` files for Hb, HbO₂, CuSO₄, NiSO₄, etc.
* `Acquisition/` – raw data I/O
* `Reconstruction/` – PA/US Fourier-domain recon (`rekon_oa_freqdom.m`, `rekon_us_freqdom1.m`)
* `Phantoms/` – digital phantom generators
* `data/` – NIfTI phantoms and processed `.mat` files
* `scripts/` – runnable examples (`run_linear_unmixing_demo.m`, `run_convex_cone_unmixing_demo.m`)
* `legacy/` – archived experimental scripts
* `Analysis/` – spectrum extraction & comparison
* `Compensation/` – linear unmixing demo (`example_linear_unmixing.m`) and convex-cone unmixing

---

## Common Workflows

### 1. Data Reconstruction

```matlab
% Photoacoustic (OA) reconstruction
reconstruction/rekon_oa_freqdom('rawPA.raw','outPA.mat');

% Ultrasound reconstruction
reconstruction/rekon_us_freqdom1('rawUS.raw','outUS.mat');

% Overlay PA + US
reconstruction/us_pa_overlay('outPA.mat','outUS.mat','overlayRGB.avi');
```

> Uses the Fourier‐domain k-space interpolation method described in Jaeger *et al.* (Inverse Problems, 2007).

### 2. Laser‐Energy Compensation

```matlab
% Correct per‐wavelength, per‐pulse energy fluctuations
reconstruction/reconstruct_laser_sweep('rawSweepFolder','calibratedPA.mat');
```

### 3. Spectrum Extraction

```matlab
[\u03bb, y_obs, \u03c3_noise] = analysis/extract_spectrum('calibratedPA.mat', ROI_mask);
```

* Computes mean ± SD across the ROI at each wavelength.

### 4. Baseline Linear Unmixing

```matlab
sO2_lin = compensation/example_linear_unmixing(\u03bb, y_obs, [\u03b5_HbO2, \u03b5_Hb]);
```

**Observed failure modes**:

* **Uniform intralipid + CuSO₄/NiSO₄ phantom** → SO₂ MAE ≓ 17%
* **Tube diameter 0.3 mm** too small → negligible spectral coloring → biased fit
* **Depth dependence**: MAE ↑ from ∼0.1 at 4 mm to ∼0.4 at 20 mm
* **Human in vivo**: irregular vessels + non-uniform background → unpredictable errors

---

## 5. Convex‐Cone–Constrained Unmixing

### 5.1 Why Convex Cone?

A cone of realistic fluence‐distorted spectra enforces physically plausible mixtures under heterogeneous illumination, mitigating depth- and background-dependent bias.

### 5.2 Phantom Generation

**Requirements:**

1. **Geometry**

   * Import open-source vascular trees (e.g. Visible Human angiogram).
   * Embed 1–2 mm vessel networks in intralipid background.
   * Use voxel size ≤ 100 µm (to respect 5–15 MHz array resolution).

```matlab
phantoms/generate_phantom_network( ...
  'vesselMask.nii',...
  'Substance_spectra/spectrum_CuSO4.mat',...
  'Substance_spectra/spectrum_NiSO4.mat',...
  'phantoms/thick_vessel.raw');
```

2. **Optical layers**

   * Vessels → blood spectra at specified sO₂ levels.
   * Background → intralipid + uniform CuSO₄/NiSO₄.

```matlab
vol = niftiread('phantoms/thick_vessel.raw');
volumeViewer(vol);
% Confirm vessel diameters ~1–2 mm, background ≥5 mm thick.
```

### 5.3 Fluence Sampling

```matlab
% MATLAB wrapper for MCX or diffusion approximation:
sample_fluence('phantoms/thick_vessel.raw','Training_data/F_RS.mat','diffusion',100);
```

* Outputs `F_RS` (n_rays × n_wavelengths) and `sigma_train`.

### 5.4 Building the Cone

```matlab
GRID_SO2 = 0:0.005:1.0;
pa_cones = build_pa_cones(F_RS, sigma_train, GRID_SO2);
```

* Whiten & unit-normalize each ray.
* Optional k-means pruning (in `Toolbox/kmeans_prune.m`).

### 5.5 Cone Unmixing

```matlab
% Whiten observed spectrum
y_w = y_obs ./ sigma_train;

% Solve α ≥ 0, sum(α)=1 to minimize ||y_w – C*α||₂
[α, idx] = cone_unmix(pa_cones, y_w);

% Recover sO₂
sO2_CC = GRID_SO2(idx);
```

**Performance**:

| Method                     | Phantom MAE | Depth dependence (4→20 mm) |
| -------------------------- | ----------: | -------------------------: |
| Linear unmixing            |        17 % |                ↑ 0.1 → 0.4 |
| Generic convex cone (Wu)   |         5 % |                 moderate ↑ |
| **Customized convex cone** |     **3 %** |             ≤0.05 constant |

---

## 6. Analysis & Plotting

* `analysis/compare_depth_errors.m` – bar plots of MAE vs depth for linear vs cone methods.
* `analysis/plot_curve_fits.m` – overlays true vs fitted spectra at varying depths and sO₂.

---

## 7. Example Scripts

Run the demo workflows from the `scripts` folder:

```matlab
% Linear unmixing example
scripts/run_linear_unmixing_demo

% Convex-cone constrained unmixing
scripts/run_convex_cone_unmixing_demo
```

### Python implementation

The repository now contains an example Python version of the convex-cone
unmixing algorithm. The script `python/convex_cone_so2.py` demonstrates how to
compute an sO₂ estimate using NumPy. Environment paths for the external forward
models are loaded via the SIMPA `PathManager` class.

---

## 8. Future Work

* **In vivo extension**: register PA+US anatomy → patient-specific cones.
* **GPU acceleration**: MCXLab or GPU-diffusion for faster fluence sampling.
* **Automated phantom generation** from US/MRI segmentations.

---

## Citation

1. Jaeger, M. *et al.* “Fourier reconstruction in optoacoustic imaging using complex k-space interpolation.” Inverse Problems 23 (2007).
2. Wu, C. *et al.* “Blood oxygenation quantification via convex cone approach.” IEEE TMI (2025).

```matlab
% In your publications, please cite both references above.
```


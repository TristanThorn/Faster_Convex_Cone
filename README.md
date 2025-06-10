# CC Unmixing

This repository contains MATLAB scripts for photoacoustic (PA) and ultrasound (US) imaging workflows. It includes routines for frequency–domain reconstruction, digital phantom creation and spectral analysis.

## MATLAB Requirements

The code was developed for MATLAB R2018a or later. Several scripts rely on the following toolboxes:

- **Image Processing Toolbox** – required for functions such as `imshow`, `volumeViewer`, `mean2`, and `std2`.
- **Statistics and Machine Learning Toolbox** (optional) – used in a few analysis scripts.

Add the repository and the `Toolbox` directory to your MATLAB path before running any scripts:

```matlab
addpath(genpath('path/to/cc_unmixing'));
```

## Directory Overview

- `Toolbox/` – helper functions for spectral unmixing and geometry calculations.
- `Substance_spectra/` – reference absorption spectra for blood, CuSO4, NiSO4, and other materials (`.mat` files).
- Key MATLAB scripts are organized into subfolders:
  - `acquisition/` – raw data loading utilities.
  - `reconstruction/` – PA/US image reconstruction workflows.
  - `phantoms/` – digital phantom generation scripts.
  - `analysis/` – spectral analysis helpers.
  - `compensation/` – angular compensation and unmixing examples.

## Common Workflows

### 1. Data Reconstruction

Raw PA or US frames acquired with a linear array transducer can be reconstructed with the frequency–domain algorithms contained in `reconstruction/rekon_oa_freqdom.m` and `reconstruction/rekon_us_freqdom1.m`. Example usage is shown in:

- `reconstruction/us_pa_overlay.m` – overlays US and PA data from `.raw` files to create RGB frames.
- `reconstruction/reconstruct_laser_spectrum.m` and `reconstruction/reconstruct_laser_sweep.m` – reconstruct wavelength sweeps and save each reconstructed frame to disk.

Update the file paths inside these scripts to point to your raw acquisitions and then run them in MATLAB.

### 2. Phantom Generation

Digital phantoms are created using the scripts in `phantoms/`. They construct layered volumes with embedded cylinders that mimic blood vessels. Example:

```matlab
phantoms/generate_phantom_basic;      % generic layered tissue phantom
phantoms/generate_phantom_human;      % human-like tissue layers
phantoms/generate_phantom_simple;     % simplified phantom geometry
```

The resulting volume is visualized using `volumeViewer` and can be saved to a `.mat` or NIfTI file.

### 3. Spectral Analysis

Spectral unmixing and plotting utilities are provided in:

- `analysis/extract_spectrum.m` – calculates signal-to-noise ratio spectra from reconstructed images.
- `analysis/plot_reference_spectra.m` – plots absorption spectra of Hb/HbO2 or phantom materials.
- `analysis/compare_spectra.m` – compares measured spectra with reference spectra.
- `compensation/example_linear_unmixing.m` – demonstrates simple SO2 estimation via linear unmixing using the functions in `Toolbox`.

These scripts load spectra from the `Substance_spectra` directory and generate figures of the processed spectra.

## Future Work

- Conduct additional phantom experiments using CuSO4 and NiSO4 to mimic HbO2 and Hb, respectively. Mix these chromophores with intralipid to better emulate tissue scattering and absorption and to obtain statistically significant results.
- Reduce the simulation time for the diffusion approximation routines.
- Validate the customized convex cone method with in-vivo experiments.
- Automatically generate PA phantoms based on registered PA and US data.

## Citation

The reconstruction functions are based on the Fourier domain method described by Jaeger *et&nbsp;al.* in
"Fourier reconstruction in optoacoustic imaging using complex k-space interpolation," *Inverse Problems* 23 (2007). Lines within the code reference this algorithm:

```matlab
% This algorithm is based on the paper:
% "Fourier reconstruction in optoacoustic imaging using complex k-space
% interpolation" by Jaeger et al., Inverse Problems 23 (2007)
```

If you use this repository in academic work, please cite the above publication and acknowledge this repository.

---

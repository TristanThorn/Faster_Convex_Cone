#!/usr/bin/env python
"""
Layered phantom → SIMPA → MCX → NIfTI
Two independent detector/illumination configs in one call:
    • Half-ring curved array  (disk dark-field illumination)
    • Linear array            (identical disk illumination)
All volumes & plots saved to Simpa_Out/
"""

# ------------------------------------------------------------------ #
# 0. Imports
# ------------------------------------------------------------------ #
from pathlib import Path
from typing import Callable, List, Tuple
import numpy as np, h5py, matplotlib.pyplot as plt
from simpa import Tags
import simpa as sp
from simpa.utils import Tags
from simpa.utils.libraries.spectrum_library import Spectrum
from simpa.utils.libraries.tissue_library import MolecularComposition
from simpa.core.device_digital_twins.detection_geometries.curved_array \
     import CurvedArrayDetectionGeometry
from simpa.core.device_digital_twins.detection_geometries.linear_array \
     import LinearArrayDetectionGeometry
from simpa.core.device_digital_twins.illumination_geometries.disk_illumination \
     import DiskIlluminationGeometry
from nibabel import nifti1, loadsave


# ------------------------- Visualisation utils -------------------- #
def plot_phantom_cross_section(vol, *, axis='z', idx=None,
                               title="Phantom cross-section",
                               save_path: Path = None): # type: ignore
    ax = {'x': 0, 'y': 1, 'z': 2}[axis]
    if idx is None:
        idx = vol.shape[ax] // 2
    sl = vol[idx,:,:] if ax==0 else vol[:,idx,:] if ax==1 else vol[:,:,idx]
    plt.figure(figsize=(4.5,4))
    plt.imshow(sl.T, cmap='tab20', origin='lower')
    plt.axis('off'); plt.title(title)
    if save_path: plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()

def extract_spectra_from_voxels(vol4d: np.ndarray, masks: List[np.ndarray]):
    return np.vstack([vol4d[m].mean(axis=0) for m in masks])

def plot_vessel_spectra(spec: np.ndarray, λ: np.ndarray,
                        *, labels=None,
                        title="Vessel fluence spectra",
                        save_path: Path = None): # type: ignore
    plt.figure(figsize=(6.2,4))
    for i,s in enumerate(spec):
        plt.plot(λ, s, '-o', label=(labels[i] if labels else f'Vessel {i+1}'))
    plt.xlabel('Wavelength (nm)'); plt.ylabel('Amplitude (a.u.)')
    plt.title(title); plt.legend(); plt.tight_layout()
    if save_path: plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()

def compute_f_rs(out_dir: Path,
                 λ_nm: np.ndarray,
                 tube_masks: list[np.ndarray],
                 bkg_mask: np.ndarray,
                 λ0: int = 800):
    """Return F_RS(λ) curves for each tube/inclusion in out_dir."""
    # locate reference index
    idx0 = np.where(λ_nm == λ0)[0][0]

    # pre-allocate
    spec_tubes = np.zeros((len(tube_masks), len(λ_nm)))
    spec_bkg   = np.zeros(len(λ_nm))

    # --- load & average fluence ---------------------------------
    for i, wl in enumerate(λ_nm):
        vol = loadsave.load(out_dir / f"fluence_{wl}nm.nii.gz").get_fdata() # type: ignore

        for k, m in enumerate(tube_masks):
            spec_tubes[k, i] = vol[m].mean()

        spec_bkg[i] = vol[bkg_mask].mean()

    # --- build F_RS curves --------------------------------------
    f_rs = (spec_tubes / spec_tubes[:, idx0, None]) / (spec_bkg / spec_bkg[idx0])

    return f_rs


def plot_f_rs(λ_nm, f_rs, sNi_list, *, title="Spectral-colouring curves",
              save_path: Path | None = None):
    plt.figure(figsize=(6,4))
    for curve, sNi in zip(f_rs, sNi_list):
        plt.plot(λ_nm, curve, '-o', label=f"$s_{{Ni}}={sNi*100:.0f}\\%$")
    plt.axvline(800, ls='--', lw=0.8, color='grey')
    plt.ylabel(r"$F_{RS}(\lambda,\Delta s_{Ni})$")
    plt.xlabel("Wavelength (nm)")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
    plt.show()

# ------------------------------------------------------------------ #
#  v7.3-MATLAB (HDF5) extinction-spectrum loader                     #
# ------------------------------------------------------------------ #

def load_extinction_spectrum(mat_path: str,
                             dset_name: str,
                             wavelengths: np.ndarray | None = None):
    """
    Load an extinction spectrum stored as a single (N,1) dataset in a v7.3
    MATLAB .mat file.

    Returns
    -------
    λ_file : np.ndarray   # wavelengths [nm]
    eps     : np.ndarray   # extinction [cm⁻¹·M⁻¹]
    """
    with h5py.File(mat_path, 'r') as f:
        eps = np.array(f[dset_name]).ravel()

    # If no wavelength vector is stored, assume the standard 700-900 nm grid
    if wavelengths is None:
        λ_file = np.arange(700, 901, 10)
    else:
        λ_file = np.asarray(wavelengths)

    return λ_file, eps


# ------------------------------------------------------------------ #
#  CuSO₄ / NiSO₄ mixture (Table-II phantom)                          #
# ------------------------------------------------------------------ #


def make_cuso4_ni_mix(*,
                      C_sulfate: float,
                      sNi: float,
                      λ_nm: np.ndarray | None = None) -> MolecularComposition:
    """
    Return a SIMPA MolecularComposition for a homogeneous CuSO₄/NiSO₄ mixture.
    """
    # --- wavelength grid ----------------------------------------------------
    if λ_nm is None:
        λ_nm = np.arange(700, 901, 10)           # 21 wavelengths

    # --- load extinction spectra -------------------------------------------
    λ_cu, eps_cu_raw = load_extinction_spectrum(
        "Substance_Spectra/spectrumCuSO4_extin.mat", "spectrum_extin_CuSO4"
    )
    λ_ni, eps_ni_raw = load_extinction_spectrum(
        "Substance_Spectra/spectrumNiSO4_extin.mat", "spectrum_extin_NiSO4"
    )

    eps_cu = np.interp(λ_nm, λ_cu, eps_cu_raw)
    eps_ni = np.interp(λ_nm, λ_ni, eps_ni_raw)

    # --- Beer–Lambert absorption -------------------------------------------
    eps_mix = (1 - sNi) * eps_cu + sNi * eps_ni          # cm⁻¹  M⁻¹
    μa_cm   = C_sulfate * eps_mix * np.log(10)           # convert to 1/cm

    # --- Table-II scattering distribution ----------------------------------
    rng        = np.random.default_rng()
    μs800_cm   = rng.uniform(7.0, 25.0)                  # cm⁻¹ at 800 nm
    b_exp      = rng.uniform(0.95, 1.15)
    μs_cm      = μs800_cm * (λ_nm / 800.0)**(-b_exp)

    # --- build Spectrum objects --------------------------------------------
    spec_mua = Spectrum("CuNi μa",  λ_nm, μa_cm)
    spec_mus = Spectrum("CuNi μs'", λ_nm, μs_cm)
    spec_g   = Spectrum("g",        λ_nm, np.full_like(λ_nm, 0.9))

    # --- return MolecularComposition for segmentation mapping --------------
    return sp.TISSUE_LIBRARY.generic_tissue(mua=spec_mua,
                                            mus=spec_mus,
                                            g=spec_g)


# ------------------------------------------------------------------ #
# 1. Phantom generator
# ------------------------------------------------------------------ #
def general_phantom_generator(pixel_mm: float,
                              base_dims_mm: Tuple[int,int,int],
                              layers: list, vessels: list,
                              *, hetero_mode=None, hetero_fraction=.5,
                              dtype=np.uint8):
    fac = 1./pixel_mm
    nx,ny,nz = np.round(np.array(base_dims_mm)*fac).astype(int)
    vol = np.ones((nx,ny,nz), dtype=dtype)
    # layers
    z0=0
    for lay in layers:
        t,label = lay['thickness'], lay['value']
        if t is None: vol[:,:,z0:] = label; break
        nzl = int(round(t*fac)); vol[:,:,z0:z0+nzl] = label; z0+=nzl
    # vessels
    yy,zz = np.meshgrid(np.arange(ny), np.arange(nz), indexing='ij')
    for v in vessels:
        cy,cz = int(round(v['center_y']*fac)), int(round(v['center_z']*fac))
        ry,rz = int(round(v['radius_y']*fac)), int(round(v.get('radius_z',v['radius_y'])*fac))
        mask2d = ((yy-cy)**2)/ry**2 + ((zz-cz)**2)/rz**2 <= 1
        mask3d = np.broadcast_to(mask2d,(nx,ny,nz))
        if hetero_mode is None:
            vol[mask3d] = v['value']
        else:
            if hetero_mode == 'alternating':
                step=max(1,int(round(1/hetero_fraction)))
                alt=np.zeros(nx,bool); alt[::step]=True
            elif hetero_mode=='random':
                alt=np.random.default_rng().random(nx)<hetero_fraction
            else:
                raise ValueError("hetero_mode invalid")
            vol[np.logical_and(mask3d,np.broadcast_to(alt[:,None,None],vol.shape))]=v['value']
    return vol


# ------------------------------------------------------------------ #
# 2. SIMPA settings
# ------------------------------------------------------------------ #
def build_simpa_settings(seg, pixel_mm, *, wavelengths, out_dir: Path, class_mapping: dict | None = None):
    T = sp.TISSUE_LIBRARY
    if class_mapping is None:
        mapping = {
            1: T.epidermis(),
            2: T.dermis(),
            3: T.generic_tissue(),
            4: T.muscle(),
            5: T.blood(0.97),
            6: T.blood(0.70),
        }
    else:
        mapping = class_mapping
    sx, sy, sz = (np.array(seg.shape) * pixel_mm).tolist()

    s = sp.Settings()
    s[Tags.SIMULATION_PATH] = str(out_dir)
    s[Tags.VOLUME_NAME]     = "phantom"
    s[Tags.SPACING_MM]      = pixel_mm
    s[Tags.DIM_VOLUME_X_MM] = sx
    s[Tags.DIM_VOLUME_Y_MM] = sy
    s[Tags.DIM_VOLUME_Z_MM] = sz
    s[Tags.WAVELENGTHS]     = list(wavelengths)
    s[Tags.RANDOM_SEED]     = 1

    s.set_volume_creation_settings({
        Tags.INPUT_SEGMENTATION_VOLUME : seg,
        Tags.SEGMENTATION_CLASS_MAPPING: mapping
    })
    s.set_optical_settings({
        Tags.OPTICAL_MODEL_NUMBER_PHOTONS    : int(1e6),
        Tags.LASER_PULSE_ENERGY_IN_MILLIJOULE: 10,
        Tags.OPTICAL_MODEL_BINARY_PATH       : "/Users/cassandrayang/Documents/GitHub/Faster_Convex_Cone/mcxcl"
    })
    pipeline = [sp.SegmentationBasedAdapter(s), sp.MCXAdapter(s)]
    return s, pipeline


# ------------------------------------------------------------------ #
# 3 Build *one* illumination – reuse for both detectors              #
# ------------------------------------------------------------------ #
def make_disk_illum(vol_dim_mm):
    return DiskIlluminationGeometry(
        beam_radius_mm   = 12.5,
        device_position_mm=np.array([vol_dim_mm[0]/2,
                                      vol_dim_mm[1]/2, -0.1])
    )

# ------------------------------------------------------------------ #
# Build detectors (ONE geometry per device → use set_detection_)     #
# ------------------------------------------------------------------ #
def make_half_ring_probe(vol_dim_mm):
    """Return a PhotoacousticDevice with a 128-elt half-ring array."""
    dev = sp.PhotoacousticDevice(device_position_mm=np.zeros(3))

    r_mm, n = 30, 128
    pitch   = np.pi * r_mm / (n - 1)
    half_ring = CurvedArrayDetectionGeometry(
        pitch_mm=pitch,
        radius_mm=r_mm,
        number_detector_elements=n,
        device_position_mm=np.array([vol_dim_mm[0]/2,
                                     vol_dim_mm[1]/2, 0])
    )
    dev.set_detection_geometry(half_ring)           # ← key line
    dev.add_illumination_geometry(make_disk_illum(vol_dim_mm))
    return dev


def make_linear_probe(vol_dim_mm):
    """Return a PhotoacousticDevice with a 128-elt linear array."""
    dev = sp.PhotoacousticDevice(device_position_mm=np.zeros(3))

    lin_len_mm, n = 38, 128
    pitch = lin_len_mm / (n - 1)
    linear = LinearArrayDetectionGeometry(
        pitch_mm=pitch,
        number_detector_elements=n,
        device_position_mm=np.array([vol_dim_mm[0]/2,
                                     vol_dim_mm[1]/2 + 12,   # shift in y
                                     12])                    # shift in z
    )
    dev.set_detection_geometry(linear)             # ← key line
    dev.add_illumination_geometry(make_disk_illum(vol_dim_mm))
    return dev

# ------------------------------------------------------------------ #
# 4. NIfTI writer
# ------------------------------------------------------------------ #
def find_dataset(h5obj, dname):
    if isinstance(h5obj, h5py.Dataset):
        if h5obj.name.split('/')[-1] == dname: # type: ignore
            return h5obj
        raise TypeError
    if dname in h5obj:
        return h5obj[dname]
    for item in h5obj.values():
        if isinstance(item, h5py.Group):
            try:    return find_dataset(item, dname)
            except KeyError: pass
    raise KeyError(dname)

def get_flu(h5file, wl_nm):
    flu_tag = Tags.DATA_FIELD_FLUENCE
    dset = find_dataset(h5file, flu_tag)
    if isinstance(dset, h5py.Group):
        dset = dset[str(wl_nm)]
    return dset[()] # type: ignore

def get_mua(h5file, wl_nm):
    mua_tag = Tags.DATA_FIELD_ABSORPTION_PER_CM
    dset = find_dataset(h5file, mua_tag)
    if isinstance(dset, h5py.Group):
        return dset[str(wl_nm)][()] # type: ignore
    return dset[()]

def write_nii(arr, voxel_mm, fname: Path, *, dtype=None):
    fname.parent.mkdir(parents=True, exist_ok=True)
    if dtype is not None: 
        arr = arr.astype(dtype)
    hdr = nifti1.Nifti1Header()
    hdr.set_data_dtype(arr.dtype)
    hdr.set_zooms((voxel_mm,)*3)
    loadsave.save(nifti1.Nifti1Image(arr, np.eye(4), hdr), str(fname))



# ------------------------------------------------------------------ #
# 5. vessel mask utility
# ------------------------------------------------------------------ #
def vessel_mask(shape, pixel_mm, v):
    nx,ny,nz=shape; fac=1./pixel_mm
    cy,cz=int(round(v['center_y']*fac)), int(round(v['center_z']*fac))
    ry,rz=int(round(v['radius_y']*fac)), int(round(v.get('radius_z',v['radius_y'])*fac))
    yy,zz=np.meshgrid(np.arange(ny), np.arange(nz), indexing='ij')
    m=((yy-cy)**2)/ry**2+((zz-cz)**2)/rz**2<=1
    return np.broadcast_to(m,(nx,ny,nz))

# ------------------------------------------------------------------ #
# 6. Run one full simulation + export + figures                      # 
# ------------------------------------------------------------------ #
# ------------------------------------------------------------------ #
# 6. Run one full simulation + export + figures                      #
# ------------------------------------------------------------------ #
def run_case(case_name: str,
             build_device: Callable[[np.ndarray], sp.PhotoacousticDevice],
             seg: np.ndarray, pixel: float,
             λ: np.ndarray, vessels: list,
             *,
             class_mapping: dict | None = None,
             sNi_list: list[float] | None = None):

    out = Path("Simpa_Out") / case_name
    out.mkdir(parents=True, exist_ok=True)

    # ---------- settings & pipeline ----------------------------------- #
    settings, pipeline = build_simpa_settings(seg, pixel,
                                              wavelengths=λ,
                                              out_dir=out,
                                              class_mapping=class_mapping)
    device = build_device(np.array(seg.shape) * pixel)

    # ----------- quick visual cross-sections -------------------------- #
    plot_phantom_cross_section(seg, axis='z', idx=50,
                               title="Layers near surface",
                               save_path=out / "layers_surface.png")

    mid_x = seg.shape[0] // 2
    plot_phantom_cross_section(seg, axis='x', idx=mid_x,
                               title="Blood-vessel cross-sections",
                               save_path=out / "vessels_cross.png")

    # ---------- vessel masks ------------------------------------------ #
    v_masks = [vessel_mask(seg.shape, pixel, v) for v in vessels]
    v_spec = np.zeros((len(v_masks), len(λ)))

    # ---------- run SIMPA + MCX --------------------------------------- #
    sp.simulate(pipeline, settings, device)
    h5f = out / "phantom.hdf5"
    assert h5f.exists()

    # ---------- export NIfTI & per-λ PNGs ----------------------------- #
    mid_z = int(round(37 / pixel))
    write_nii(seg, pixel, out / "segmentation.nii.gz", dtype=np.uint8)

    with h5py.File(h5f, "r") as f:
        for i, wl in enumerate(λ):
            flu = get_flu(f, str(wl))
            mua = get_mua(f, str(wl))
            write_nii(flu, pixel, out / f"fluence_{wl}nm.nii.gz")
            write_nii(mua, pixel, out / f"mu_a_{wl}nm.nii.gz")

            for k, m in enumerate(v_masks):
                v_spec[k, i] = flu[m].mean() #type: ignore

            # slice PNG
            plt.figure(figsize=(7.4, 2.4))
            plt.subplot(1, 3, 1); plt.imshow(seg[:, :, mid_z], cmap='tab20');     plt.axis('off'); plt.title("Seg")
            plt.subplot(1, 3, 2); plt.imshow(flu[:, :, mid_z], cmap='inferno');   plt.axis('off'); plt.title(f"Flu {wl}") # type: ignore
            plt.subplot(1, 3, 3); plt.imshow(mua[:, :, mid_z], cmap='viridis');   plt.axis('off'); plt.title(f"μa {wl}") # type: ignore
            plt.tight_layout(); plt.savefig(out / f"slice_{wl}nm.png", dpi=300); plt.close()

    plot_vessel_spectra(v_spec, λ,
                        labels=[f"Vessel {i+1}" for i in range(len(v_masks))],
                        save_path=out / "vessel_spectra.png")

    # ---------- optional F_RS curves ---------------------------------- #
    if sNi_list is not None:
        bkg_mask = seg == 1                       # background = label 1
        f_rs = compute_f_rs(out, λ, v_masks, bkg_mask)
        plot_f_rs(λ, f_rs, sNi_list,
                  save_path=out / "F_RS_curves.png")

    print(f"[{case_name}] done → {out.resolve()}")



# ------------------------------------------------------------------ #
# 7. Main – run both cases                                              #
# ------------------------------------------------------------------ #
# ------------------------------------------------------------------ #
# 7. Main – build phantom & run both detector configs                #
# ------------------------------------------------------------------ #
def main():
    pixel = 0.1
    dims  = (61, 61, 31)            # mm (x, y, z)
    λ     = np.arange(700, 901, 10)

    # ----- 5 high-concentration target tubes ------------------------ #
    tubes = [dict(center_y=15 + i*5,
                  center_z=10,
                  radius_y=2,
                  value=5) for i in range(5)]

    # ----- 7 colouring inclusions each with its own sNi -------------- #
    rng       = np.random.default_rng()
    incls     = []
    sNi_list  = []                       # store Ni fraction per inclusion
    for i in range(7):
        frac = rng.uniform(0.0, 1.0)
        incls.append(dict(center_y=10 + i*7,
                          center_z=18,
                          radius_y=3,
                          value=6,
                          sNi=frac))
        sNi_list.append(frac)

    seg = general_phantom_generator(
        pixel_mm=pixel,
        base_dims_mm=dims,
        layers=[dict(thickness=None, value=1)],    # homogeneous slab
        vessels=tubes + incls
    )

    # ----- SIMPA tissue mapping ------------------------------------- #
    T = sp.TISSUE_LIBRARY
    class_mapping = {
        1: T.generic_tissue(),                               # agarose slab
        5: make_cuso4_ni_mix(C_sulfate=0.144, sNi=0.0),      # targets
        6: make_cuso4_ni_mix(C_sulfate=0.014, sNi=np.mean(sNi_list))  # inclusions # type: ignore
    }

    # ----- Run Half-Ring and Linear-Array cases --------------------- #
    vessel_defs = tubes + incls
    all_sNi     = [0.0]*5 + sNi_list

    run_case("Half_Ring",    make_half_ring_probe,
             seg, pixel, λ, vessel_defs,
             class_mapping=class_mapping,
             sNi_list=all_sNi)

    run_case("Linear_Array", make_linear_probe,
             seg, pixel, λ, vessel_defs,
             class_mapping=class_mapping,
             sNi_list=all_sNi)


if __name__ == "__main__":
    main()
#!/usr/bin/env python
"""
Layered phantom → SIMPA half-ring probe → MCX fluence → NIfTI export
Multispectral version: 700 - 900 nm in 10 nm steps, per-λ .nii files + slice preview.
Tested with:  simpa 1.x, numpy, nibabel, h5py   (MCX installed & on PATH)
"""

# ------------------------------------------------------------------ #
# 0. Imports
# ------------------------------------------------------------------ #
from pathlib import Path
from typing import Union
import numpy as np, h5py, matplotlib.pyplot as plt
import simpa as sp
from simpa.utils import Tags
from simpa.core.device_digital_twins.detection_geometries.curved_array \
     import CurvedArrayDetectionGeometry
from simpa.core.device_digital_twins.illumination_geometries.disk_illumination \
     import DiskIlluminationGeometry
from nibabel import nifti1, loadsave


# ------------------------------------------------------------------ #
# 1. Phantom generator (vectorised, broadcast-safe)
# ------------------------------------------------------------------ #
def general_phantom_generator(pixel_size_mm, base_dims, layers, vessels,
                              *, hetero_mode=None, hetero_fraction=0.5,
                              dtype=np.uint8):
    fac = 0.1 / pixel_size_mm
    nx, ny, nz = np.round(np.asarray(base_dims) * fac).astype(int)
    phantom = np.ones((nx, ny, nz), dtype=dtype)

    # -- horizontal layers --
    z0 = 0
    for lay in layers:
        t = lay["thickness"]
        z1 = nz if t in (None, [], np.nan) else min(z0 + int(round(t * fac)), nz)
        phantom[:, :, z0:z1] = lay["value"]
        z0 = z1
        if z0 >= nz:
            break

    # -- vessels --
    yy, zz = np.meshgrid(np.arange(ny), np.arange(nz), indexing="ij")
    for v in vessels:
        cy, cz = int(round(v["center_y"] * fac)), int(round(v["center_z"] * fac))
        ry = int(round(v["radius_y"] * fac))
        rz = int(round(v.get("radius_z", v["radius_y"]) * fac))
        label = v["value"]

        mask2d = ((yy - cy) ** 2) / ry**2 + ((zz - cz) ** 2) / rz**2 <= 1
        mask3d = np.broadcast_to(mask2d, (nx, ny, nz))

        if hetero_mode is None:
            phantom[mask3d] = label
        else:
            if hetero_mode == "alternating":
                step = max(1, int(round(1 / hetero_fraction)))
                alt = np.zeros(nx, bool); alt[::step] = True
            elif hetero_mode == "random":
                alt = np.random.rand(nx) < hetero_fraction
            else:
                raise ValueError("hetero_mode must be None | alternating | random")

            phantom[np.logical_and(mask3d,  alt[:, None, None])] = label + 1
            phantom[np.logical_and(mask3d, ~alt[:, None, None])] = label
    return phantom


# ------------------------------------------------------------------ #
# 2. SIMPA settings (Segmentation-based workflow, MCX only)
# ------------------------------------------------------------------ #
def build_simpa_settings(seg, pixel_size_mm, *, wavelengths,
                         out_dir="simpa_out"):
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    T = sp.TISSUE_LIBRARY
    mapping = {
        1: T.epidermis(),
        2: T.dermis(),
        3: T.generic_tissue(),
        4: T.muscle(),
        5: T.blood(0.97),
        6: T.blood(0.70),                    # heterogeneity label
    }

    sx, sy, sz = (np.array(seg.shape) * pixel_size_mm).tolist()

    s = sp.Settings()
    s[Tags.SIMULATION_PATH] = out_dir
    s[Tags.VOLUME_NAME]     = "phantom"
    s[Tags.SPACING_MM]      = pixel_size_mm
    s[Tags.DIM_VOLUME_X_MM] = sx
    s[Tags.DIM_VOLUME_Y_MM] = sy
    s[Tags.DIM_VOLUME_Z_MM] = sz
    s[Tags.WAVELENGTHS]     = list(wavelengths)
    s[Tags.RANDOM_SEED]     = 1

    s.set_volume_creation_settings({
        Tags.INPUT_SEGMENTATION_VOLUME : seg,
        Tags.SEGMENTATION_CLASS_MAPPING: mapping,
    })

    s.set_optical_settings({
        Tags.OPTICAL_MODEL_NUMBER_PHOTONS    : int(1e6),
        Tags.LASER_PULSE_ENERGY_IN_MILLIJOULE: 10,
        Tags.OPTICAL_MODEL_BINARY_PATH       : "/Users/cassandrayang/Documents/GitHub/Faster_Convex_Cone/mcxcl"
    })

    pipeline = [sp.SegmentationBasedAdapter(s), sp.MCXAdapter(s)]
    return s, pipeline


# ------------------------------------------------------------------ #
# 3. Half-ring probe with dark-field disk source
# ------------------------------------------------------------------ #
class HalfRingDarkFieldProbe(sp.PhotoacousticDevice):
    """128-element half-ring (≈180°) with a 25 mm disk source."""
    def __init__(self, vol_dim_mm):
        super().__init__(device_position_mm=np.zeros(3))
        radius_mm, n = 30, 128
        pitch_mm = np.pi * radius_mm / (n - 1)

        self.set_detection_geometry(
            CurvedArrayDetectionGeometry(
                pitch_mm                = pitch_mm,
                radius_mm               = radius_mm,
                number_detector_elements= n,
                device_position_mm      = np.array([vol_dim_mm[0]/2,
                                                    vol_dim_mm[1]/2, 0])
            )
        )
        self.add_illumination_geometry(
            DiskIlluminationGeometry(
                beam_radius_mm   = 12.5,
                device_position_mm=np.array([vol_dim_mm[0]/2,
                                              vol_dim_mm[1]/2, -0.1])
            )
        )


# ------------------------------------------------------------------ #
# 4. HDF5 helpers
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
    flu_tag = getattr(Tags, "DATA_FIELD_OPTICAL_FLUENCE",
                      getattr(Tags, "DATA_FIELD_FLUENCE", "fluence"))
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

def write_nii(arr, vox, fname, dtype: Union[np.dtype, type] = np.float16):
    affine = np.diag([vox, vox, vox, 1])
    img    = nifti1.Nifti1Image(arr.astype(dtype), affine)
    loadsave.save(img, fname)
    print("→", fname)


# ------------------------------------------------------------------ #
# 5. Main
# ------------------------------------------------------------------ #
def main():
    pixel = 0.1  # mm
    wavelengths = np.arange(700, 901, 10)      # 700–900 nm

    # -- segmentation volume ------------------------------------------------
    seg = general_phantom_generator(
        pixel_size_mm = pixel,
        base_dims     = (610, 200, 650),
        layers        = [
            {"thickness": 40,  "value": 1},
            {"thickness": 20,  "value": 2},
            {"thickness": 40,  "value": 3},
            {"thickness": None,"value": 4}
        ],
        vessels=[
            {"center_y": 65, "center_z": 37, "radius_y": 3, "value": 5}
        ],
        hetero_mode="alternating", hetero_fraction=0.5
    )

    # -- SIMPA settings & run ------------------------------------------------
    settings, pipeline = build_simpa_settings(seg, pixel, wavelengths=wavelengths)
    device = HalfRingDarkFieldProbe(np.array(seg.shape) * pixel)
    sp.simulate(pipeline, settings, device)

    # -- open MCX output -----------------------------------------------------
    h5 = Path(settings[Tags.SIMULATION_PATH], "phantom.hdf5")
    outdir = Path(settings[Tags.SIMULATION_PATH])
    write_nii(seg, pixel, outdir / "segmentation.nii.gz", dtype=np.uint8)

    mid_z = seg.shape[2] // 2  # slice index for quick preview

    with h5py.File(h5) as f:
        for wl in wavelengths:
            flu = get_flu(f, wl)
            mua = get_mua(f, wl)

            write_nii(flu, pixel, outdir / f"fluence_{wl}nm.nii.gz")
            write_nii(mua, pixel, outdir / f"mu_a_{wl}nm.nii.gz")

            # -------- quick visual verification --------
            plt.figure(figsize=(8, 2.5))
            plt.subplot(1, 3, 1)
            plt.imshow(seg[:, :, mid_z], cmap="tab20")
            plt.title("Seg") ; plt.axis("off")

            plt.subplot(1, 3, 2)
            plt.imshow(flu[:, :, mid_z], cmap="inferno") # type: ignore
            plt.title(f"Fluence {wl}nm") ; plt.axis("off")

            plt.subplot(1, 3, 3)
            plt.imshow(mua[:, :, mid_z], cmap="viridis") # type: ignore
            plt.title(f"μa {wl}nm") ; plt.axis("off")

            plt.tight_layout()
            plt.show()

    print("All NIfTI volumes saved to", outdir.resolve())


if __name__ == "__main__":
    main()



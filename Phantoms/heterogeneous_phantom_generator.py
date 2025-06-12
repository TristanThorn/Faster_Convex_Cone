#!/usr/bin/env python
"""
Layered phantom  ➜ SIMPA half-ring probe  ➜ MCX fluence ➜ .nii export
Tested with:  simpa 1.x,  numpy,  nibabel,  h5py   (MCX installed & in PATH)
"""

# --------------------------------------------------------------------- #
# 0. Imports
# --------------------------------------------------------------------- #
from pathlib import Path
import numpy as np, h5py
import matplotlib.pyplot as plt
import simpa as sp
from simpa.utils import Tags
# detection & illumination classes – import directly from the module files
from simpa.core.device_digital_twins.detection_geometries.curved_array \
    import CurvedArrayDetectionGeometry
from simpa.core.device_digital_twins.illumination_geometries.disk_illumination \
    import DiskIlluminationGeometry
from nibabel import nifti1, loadsave


# --------------------------------------------------------------------- #
# 1. Phantom generator (vectorised, broadcast-safe)
# --------------------------------------------------------------------- #
def general_phantom_generator(pixel_size_mm, base_dims, layers, vessels,
                              *, hetero_mode=None, hetero_fraction=0.5,
                              dtype=np.uint8):
    fac = 0.1 / pixel_size_mm
    nx, ny, nz = np.round(np.asarray(base_dims) * fac).astype(int)
    phantom = np.ones((nx, ny, nz), dtype=dtype)

    # ---- horizontal layers ----
    z0 = 0
    for lay in layers:
        t = lay["thickness"]
        z1 = nz if t in (None, [], np.nan) else min(z0 + int(round(t * fac)), nz)
        phantom[:, :, z0:z1] = lay["value"]
        z0 = z1
        if z0 >= nz:
            break

    # ---- vessels ----
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


# --------------------------------------------------------------------- #
# 2. SIMPA settings (Segmentation-based workflow, MCX only)
# --------------------------------------------------------------------- #
def build_simpa_settings(seg, pixel_size_mm, *, wavelengths=(800,),
                         out_dir="simpa_out"):
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    T = sp.TISSUE_LIBRARY
    mapping = {
        1: T.epidermis(),
        2: T.dermis(),
        3: T.generic_tissue(),
        4: T.muscle(),
        5: T.blood(0.97),
        6: T.blood(0.70),
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
    s[Tags.RANDOM_SEED]     = 1            # reproducible run

    s.set_volume_creation_settings({
        Tags.INPUT_SEGMENTATION_VOLUME : seg,
        Tags.SEGMENTATION_CLASS_MAPPING: mapping,
    })

    s.set_optical_settings({
        Tags.OPTICAL_MODEL_NUMBER_PHOTONS    : int(1e6),
        Tags.LASER_PULSE_ENERGY_IN_MILLIJOULE: 10,
        Tags.OPTICAL_MODEL_BINARY_PATH: "mcxlab/mcxcl"
    })

    pipeline = [sp.SegmentationBasedAdapter(s), sp.MCXAdapter(s)]
    return s, pipeline


# --------------------------------------------------------------------- #
# 3. Half-ring probe with dark-field disk source
# --------------------------------------------------------------------- #
class HalfRingDarkFieldProbe(sp.PhotoacousticDevice):
    """
    128-element half-ring (≈180°) with a 25 mm disk source.
    """
    def __init__(self, vol_dim_mm):
        super().__init__(device_position_mm=np.zeros(3))

        # -- curved array ---------------------------------------------------
        radius_mm = 30
        n         = 128
        # pitch so that (n-1)·pitch_angle ≈ π  (half-ring)
        pitch_mm  = np.pi * radius_mm / (n - 1)

        self.set_detection_geometry(
            CurvedArrayDetectionGeometry(
                pitch_mm                = pitch_mm,
                radius_mm               = radius_mm,
                number_detector_elements= n,
                device_position_mm      = np.array([vol_dim_mm[0]/2,
                                                    vol_dim_mm[1]/2, 0])
            )
        )

        # -- disk illumination ---------------------------------------------
        self.add_illumination_geometry(
            DiskIlluminationGeometry(
                beam_radius_mm   = 12.5,                          # 25 mm Ø
                device_position_mm=np.array([vol_dim_mm[0]/2,
                                              vol_dim_mm[1]/2, -0.1])
            )
        )


# --------------------------------------------------------------------- #
# 4. Helpers
# --------------------------------------------------------------------- #
# NiBabel helper using nifti1 + loadsave
def find_dataset(h5obj, dname):
    if isinstance(h5obj, h5py.Dataset):
        if h5obj.name.split('/')[-1] == dname: # type: ignore
            return h5obj
        raise TypeError
    if dname in h5obj:
        return h5obj[dname]
    for item in h5obj.values():
        if isinstance(item, h5py.Group):
            try:
                return find_dataset(item, dname)
            except KeyError:
                pass
    raise KeyError(dname)

def get_fluence(h5file, wl_nm=800):
    # resolve old / new tag spelling
    flu_tag = getattr(Tags, "DATA_FIELD_OPTICAL_FLUENCE",
                      getattr(Tags, "DATA_FIELD_FLUENCE", "fluence"))
    dset = find_dataset(h5file, flu_tag)        # may still be a Group
    if isinstance(dset, h5py.Group):
        dset = dset[str(wl_nm)]                 # wavelength sub-group
    return dset[()]                             # type: ignore # numpy array

def get_mua(h5obj, wavelength_nm=800):
    """
    Return μa [1/cm] for the requested wavelength, regardless of hierarchy.
    """
    mua_tag = Tags.DATA_FIELD_ABSORPTION_PER_CM  # 'mua'

    dset_or_group = find_dataset(h5obj, mua_tag)

    # In recent wheels 'mua' is a *group* with one dataset per wavelength
    if isinstance(dset_or_group, h5py.Group):
        return dset_or_group[str(wavelength_nm)][()]    # type: ignore # → ndarray

    # Old wheels stored it as a single dataset
    return dset_or_group[()]                            # ndarray


def get_segmentation(volume_h5: h5py.File):
    """
    Retrieve the integer segmentation that was written by
    SegmentationBasedAdapter (same trick as above).
    """
    seg_tag = getattr(Tags, "DATA_FIELD_SEGMENTATION", "seg")
    return volume_h5[seg_tag][()] # type: ignore

def write_nii(arr, vox, fname, dtype=np.float32):
    affine = np.diag([vox, vox, vox, 1])
    img    = nifti1.Nifti1Image(arr.astype(dtype), affine)
    loadsave.save(img, fname)     # same as nib.save()
    print("→", fname)



# --------------------------------------------------------------------- #
# 5. Main
# --------------------------------------------------------------------- #
def main():
    pixel = 0.1          # mm

    seg = general_phantom_generator(
        pixel_size_mm=pixel,
        base_dims=(610, 200, 650),
        layers=[{"thickness": 40,  "value": 1},
                {"thickness": 20,  "value": 2},
                {"thickness": 40,  "value": 3},
                {"thickness": None,"value": 4}],
        vessels=[{"center_y": 65, "center_z": 37,
                  "radius_y": 3,  "value": 5}],
        hetero_mode="alternating", hetero_fraction=0.5)

    settings, pipeline = build_simpa_settings(seg, pixel_size_mm=pixel)
    device = HalfRingDarkFieldProbe(np.array(seg.shape) * pixel)

    sp.simulate(pipeline, settings, device)

    # ---------- export ----------
    h5 = Path(settings[Tags.SIMULATION_PATH], "phantom.hdf5")
    with h5py.File(h5) as f:
        flu = get_fluence(f)
        mua = get_mua(f) # type: ignore
    outdir = Path(settings[Tags.SIMULATION_PATH])
    write_nii(seg, pixel, outdir / "segmentation.nii")
    write_nii(flu, pixel, outdir / "fluence_800nm.nii")
    write_nii(mua, pixel, outdir / "mu_a.nii")
    print("Saved NIfTI volumes to", outdir.resolve())


if __name__ == "__main__":
    main()


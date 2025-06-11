"""Digital phantom generators using SIMPA's VolumeCreationModule.

This script reproduces the MATLAB phantom generation scripts in the
``Phantoms`` directory with SIMPA. Only the volume creation step is
performed and the resulting segmentation volume is returned.
"""

from __future__ import annotations

from typing import Iterable, Optional
import numpy as np
import simpa as sp
from simpa import Tags


def _constant_tissue(segmentation_value: int) -> sp.MolecularComposition:
    """Create a dummy tissue with the desired segmentation label."""
    tissue = sp.TISSUE_LIBRARY.constant(1e-10, 1e-10, 0.9)
    tissue.segmentation_type = segmentation_value
    return tissue


def general_phantom_generator(
    pixel_size: float,
    base_dims: Iterable[int],
    layers: list[dict],
    vessels: list[dict],
    volume_name: str = "digital_phantom",
    path_manager: Optional[sp.PathManager] = None,
) -> np.ndarray:
    """Create a layered phantom with optional vessels.

    ``base_dims`` are specified in multiples of 0.1 mm, matching the
    original MATLAB conventions. ``layers`` is a list of dictionaries
    with keys ``thickness`` (in the same units) and ``value``. A thickness
    of ``None`` fills the remaining volume. ``vessels`` defines cylindrical
    or elliptical vessels with centre coordinates and radii in the same
    units.
    """
    if path_manager is None:
        path_manager = sp.PathManager()

    dims_mm = np.asarray(list(base_dims), dtype=float) * 0.1

    settings = sp.Settings({
        Tags.SIMULATION_PATH: path_manager.get_hdf5_file_save_path(),
        Tags.VOLUME_NAME: volume_name,
        Tags.SPACING_MM: pixel_size,
        Tags.DIM_VOLUME_X_MM: float(dims_mm[0]),
        Tags.DIM_VOLUME_Y_MM: float(dims_mm[1]),
        Tags.DIM_VOLUME_Z_MM: float(dims_mm[2]),
        Tags.RANDOM_SEED: 1,
        Tags.WAVELENGTHS: [700],
        Tags.GPU: False,
    })

    tissue_dict = sp.Settings()
    tissue_dict[Tags.BACKGROUND] = sp.define_background_structure_settings(
        molecular_composition=_constant_tissue(0)
    )

    current_z = 0.0
    for idx, layer in enumerate(layers):
        thickness = layer["thickness"]
        if thickness is None:
            thickness_mm = dims_mm[2] - current_z
        else:
            thickness_mm = float(thickness) * 0.1
        tissue_dict[f"layer_{idx}"] = sp.define_horizontal_layer_structure_settings(
            molecular_composition=_constant_tissue(layer["value"]),
            z_start_mm=current_z,
            thickness_mm=thickness_mm,
            priority=idx + 1,
            consider_partial_volume=True,
        )
        current_z += thickness_mm

    for i, vessel in enumerate(vessels):
        cy_mm = vessel["center_y"] * 0.1
        cz_mm = vessel["center_z"] * 0.1
        ry_mm = vessel["radius_y"] * 0.1
        rz_mm = vessel.get("radius_z", vessel["radius_y"]) * 0.1
        start = [0, cy_mm, cz_mm]
        end = [dims_mm[0], cy_mm, cz_mm]
        priority = 10 + i
        if abs(ry_mm - rz_mm) < 1e-9:
            tissue_dict[f"vessel_{i}"] = sp.define_circular_tubular_structure_settings(
                tube_start_mm=start,
                tube_end_mm=end,
                molecular_composition=_constant_tissue(vessel["value"]),
                radius_mm=ry_mm,
                priority=priority,
                consider_partial_volume=True,
            )
        else:
            major = max(ry_mm, rz_mm)
            ecc = np.sqrt(1.0 - (min(ry_mm, rz_mm) / major) ** 2)
            tissue_dict[f"vessel_{i}"] = sp.define_elliptical_tubular_structure_settings(
                tube_start_mm=start,
                tube_end_mm=end,
                molecular_composition=_constant_tissue(vessel["value"]),
                radius_mm=major,
                eccentricity=ecc,
                priority=priority,
                consider_partial_volume=True,
            )

    settings.set_volume_creation_settings({Tags.STRUCTURES: tissue_dict})

    volumes = sp.ModelBasedAdapter(settings).create_simulation_volume()
    return volumes[Tags.DATA_FIELD_SEGMENTATION]


def example_phantom_basic() -> np.ndarray:
    layers = [
        {"thickness": 10, "value": 1},
        {"thickness": 20, "value": 2},
        {"thickness": 30, "value": 3},
        {"thickness": None, "value": 4},
    ]
    vessels = [
        {"center_y": 70, "center_z": 31, "radius_y": 11, "value": 5},
        {"center_y": 190, "center_z": 53, "radius_y": 13, "value": 5},
        {"center_y": 130, "center_z": 47, "radius_y": 17, "value": 5},
        {"center_y": 230, "center_z": 68, "radius_y": 18, "value": 5},
    ]
    dims = (180, 300, 200)
    return general_phantom_generator(0.05, dims, layers, vessels, "phantom_basic")


def example_phantom_human() -> np.ndarray:
    layers = [
        {"thickness": 5, "value": 1},
        {"thickness": 15, "value": 2},
        {"thickness": 38, "value": 3},
        {"thickness": 36, "value": 4},
        {"thickness": None, "value": 3},
    ]
    vessels = [
        {"center_y": 81, "center_z": 46, "radius_y": 9, "radius_z": 4, "value": 5},
        {"center_y": 30, "center_z": 51, "radius_y": 9, "radius_z": 3, "value": 6},
        {"center_y": 150, "center_z": 65, "radius_y": 9, "radius_z": 9, "value": 7},
        {"center_y": 74, "center_z": 74, "radius_y": 10, "radius_z": 11, "value": 8},
        {"center_y": 135, "center_z": 91, "radius_y": 20, "radius_z": 17, "value": 8},
        {"center_y": 105, "center_z": 98, "radius_y": 10, "radius_z": 8, "value": 9},
        {"center_y": 115, "center_z": 147, "radius_y": 15, "radius_z": 10, "value": 10},
    ]
    dims = (120, 200, 200)
    return general_phantom_generator(0.2, dims, layers, vessels, "phantom_human")


if __name__ == "__main__":
    seg = example_phantom_basic()
    print("Generated phantom volume with shape", seg.shape)

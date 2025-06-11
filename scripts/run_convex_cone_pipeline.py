import numpy as np
from scipy.optimize import minimize
import simpa as sp
from simpa import Tags


def _unmix_vector(pa_vec, hb, hbo2):
    pa_vec = np.asarray(pa_vec, dtype=float)
    C = np.vstack([hb, hbo2]).T
    def obj(a):
        return np.linalg.norm(C @ a - pa_vec) ** 2
    cons = (
        {'type': 'eq', 'fun': lambda a: np.sum(a) - 1},
        {'type': 'ineq', 'fun': lambda a: a}
    )
    res = minimize(obj, x0=np.array([0.5, 0.5]), constraints=cons, method='SLSQP')
    return float(res.x[1])


def convex_cone_unmixing(pa_volume, wavelengths):
    spectra = sp.get_simpa_internal_absorption_spectra_by_names([
        Tags.SIMPA_NAMED_ABSORPTION_SPECTRUM_DEOXYHEMOGLOBIN,
        Tags.SIMPA_NAMED_ABSORPTION_SPECTRUM_OXYHEMOGLOBIN
    ])
    hb = np.array([spectra[0].get_value_for_wavelength(w) for w in wavelengths])
    hbo2 = np.array([spectra[1].get_value_for_wavelength(w) for w in wavelengths])
    flat = pa_volume.reshape(pa_volume.shape[0], -1)
    so2 = np.zeros(flat.shape[1])
    for idx in range(flat.shape[1]):
        so2[idx] = _unmix_vector(flat[:, idx], hb, hbo2)
    return so2.reshape(pa_volume.shape[1:])


def run_pipeline(spacing=0.5, path_manager=None):
    if path_manager is None:
        path_manager = sp.PathManager()

    wavelengths = [750, 800, 850]
    settings = sp.Settings({
        Tags.RANDOM_SEED: 4711,
        Tags.VOLUME_NAME: "cc_pipeline",
        Tags.SIMULATION_PATH: path_manager.get_hdf5_file_save_path(),
        Tags.SPACING_MM: spacing,
        Tags.DIM_VOLUME_X_MM: 30,
        Tags.DIM_VOLUME_Y_MM: 30,
        Tags.DIM_VOLUME_Z_MM: 30,
        Tags.WAVELENGTHS: wavelengths,
        Tags.DO_FILE_COMPRESSION: True,
        Tags.GPU: True
    })

    background = sp.Settings()
    background[Tags.MOLECULE_COMPOSITION] = sp.TISSUE_LIBRARY.constant(1e-4, 1e-4, 0.9)
    background[Tags.STRUCTURE_TYPE] = Tags.BACKGROUND
    tissue = sp.Settings()
    tissue[Tags.BACKGROUND] = background
    tissue["muscle"] = sp.define_horizontal_layer_structure_settings(
        z_start_mm=0, thickness_mm=30,
        molecular_composition=sp.TISSUE_LIBRARY.muscle(),
        priority=1, consider_partial_volume=True, adhere_to_deformation=True)
    tissue["vessel"] = sp.define_circular_tubular_structure_settings(
        tube_start_mm=[15, 15, 10], tube_end_mm=[15, 15, 20],
        molecular_composition=sp.TISSUE_LIBRARY.blood(),
        radius_mm=2, priority=3, consider_partial_volume=True)

    settings.set_volume_creation_settings({
        Tags.STRUCTURES: tissue,
        Tags.SIMULATE_DEFORMED_LAYERS: True
    })
    settings.set_optical_settings({
        Tags.OPTICAL_MODEL_NUMBER_PHOTONS: 1e6,
        Tags.OPTICAL_MODEL_BINARY_PATH: path_manager.get_mcx_binary_path()
    })
    settings.set_acoustic_settings({
        Tags.ACOUSTIC_MODEL_BINARY_PATH: path_manager.get_matlab_binary_path()
    })
    settings.set_reconstruction_settings(sp.create_reconstruction_settings())

    pipeline = [
        sp.ModelBasedAdapter(settings),
        sp.MCXAdapter(settings),
        sp.KWaveAdapter(settings),
        sp.DelayAndSumAdapter(settings)
    ]

    device = sp.PhotoacousticDevice(device_position_mm=np.array([15, 15, 0]))
    device.set_detection_geometry(sp.LinearArrayDetectionGeometry())
    device.add_illumination_geometry(sp.SlitIlluminationGeometry(slit_vector_mm=[5, 0, 0]))

    sp.simulate(pipeline, settings, device)

    file_path = settings[Tags.SIMPA_OUTPUT_FILE_PATH]
    pa_volume = np.asarray([
        sp.load_data_field(file_path, Tags.DATA_FIELD_RECONSTRUCTED_DATA, wl)
        for wl in wavelengths
    ])

    so2_map = convex_cone_unmixing(pa_volume, wavelengths)
    sp.save_data_field(so2_map, file_path, "convex_cone_so2")
    print(f"Convex-cone sO2 map saved to {file_path}")


if __name__ == "__main__":
    run_pipeline()

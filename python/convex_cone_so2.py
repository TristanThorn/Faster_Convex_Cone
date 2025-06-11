import numpy as np


def angle_between(a, b):
    a = np.asarray(a)
    b = np.asarray(b)
    dot = np.dot(a, b)
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    if denom == 0:
        return np.pi / 2
    cos = np.clip(dot / denom, -1.0, 1.0)
    return np.arccos(cos)


def convex_cone_so2(pa, spectra_hb, spectra_hbo2, f_base):
    s_grid = np.arange(0, 1.0001, 0.005)
    angles = []
    nearest = []
    for s in s_grid:
        mua = (1 - s) * np.log(10) * spectra_hb + s * np.log(10) * spectra_hbo2
        deposit = f_base * mua
        angles.append(angle_between(pa, deposit))
        nearest.append(deposit)
    angles = np.array(angles)
    nearest = np.stack(nearest, axis=1)
    idx = np.argmin(angles)
    return s_grid[idx], angles, nearest, np.zeros_like(angles)

if __name__ == "__main__":
    # example usage with random data
    pa = np.ones(21)
    spectra_hb = np.random.rand(21)
    spectra_hbo2 = np.random.rand(21)
    f_base = np.ones(21)
    so2, angles, nearest, flags = convex_cone_so2(pa, spectra_hb, spectra_hbo2, f_base)
    print("Estimated sO2:", so2)

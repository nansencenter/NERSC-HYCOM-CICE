"""Spreading of river discharge from point sources over nearby coastal ocean cells.

Shared by the GloFAS river forcing (river-topaz/scripts/GloFAS) and the AHYPE v4b river forcing
(river_ahype_v4b.py). A point source is a river mouth, an outlet or a GEUS ice discharge gate, already
mapped to a model ocean cell. Its weights over the ocean cells are:

  1. exp(-d / alongshore) for ocean cells within `alongshore` of the source cell (distance d), zero beyond;
  2. times exp(-2 * landdist / rradius) for cells within `rradius` of the coast (landdist), zero beyond;
  3. zero for cells whose nearest land is another landmass than the source cell's;
  4. zero for cells separated from the source by more than 15 land cells (cleaning);
  5. Gaussian smoothing over `sigma` grid cells, zero on land;
  6. normalised to a sum of 1, so each source's volume is conserved. If nothing is left, all the
     water stays in the source cell.

This is the algorithm of propagate_river_plume() in the original GloFAS scripts (Acchiardi 2025), but only
the cells within `alongshore` of each source are computed, and sources on the same cell share their weights.

Strategies (SPREADING):
  pointsource  no spreading: all the water in the source cell
  coastal      rradius = 50 km, alongshore = 100 km, sigma = 1 cell
  diffuse      rradius = 180 km, alongshore = 360 km, sigma = 3 cells (the original GloFAS setting)
"""
import time

import numpy as np
import scipy.sparse as sp
from numba import njit
from scipy.ndimage import gaussian_filter, label
from scipy.spatial import cKDTree

SPREADING = {
    'pointsource': None,
    'coastal': dict(rradius=50_000, alongshore=100_000, sigma=1),
    'diffuse': dict(rradius=180_000, alongshore=360_000, sigma=3),
}

LAND_THRESHOLD = 15   # cells further from the source than this many land cells get no weight


def ecef(lat, lon):
    """Cartesian (x, y, z) coordinates in km on the WGS 84 ellipsoid, as transform_coordinates() of the
    GloFAS scripts (single precision, so distances are the same)."""
    a, e2 = 6378.137, 6.69437999014e-3
    lat = np.radians(np.asarray(lat, dtype='f4').ravel())
    lon = np.radians(np.asarray(lon, dtype='f4').ravel())
    r_n = a / np.sqrt(1 - e2 * np.sin(lat) ** 2)
    return np.column_stack((r_n * np.cos(lat) * np.cos(lon), r_n * np.cos(lat) * np.sin(lon),
                            r_n * (1 - e2) * np.sin(lat)))


@njit
def cleaning_discharge(ocean_weight, weight_mask, mini, minj, idm, jdm, land):
    """Zero the weight of cells separated from the source cell (mini, minj) by more than LAND_THRESHOLD
    land cells, counted along a path first in i, then in j (from the GloFAS scripts)."""
    shape = (jdm, idm)
    for idx in weight_mask[0]:
        coords = np.empty(2, dtype=np.int64)
        index = idx
        for i in range(1, -1, -1):
            coords[i] = index % shape[i]
            index = index // shape[i]
        ii, jj = coords
        if ii != mini:
            list_ii = np.arange(ii, mini, 1 * np.sign(mini - ii))
        else:
            list_ii = np.array([ii])
        if jj != minj:
            list_jj = np.arange(jj, minj, 1 * np.sign(minj - jj))
        else:
            list_jj = np.array([jj])
        len1 = list_ii.shape[0] + list_jj.shape[0]
        indices = np.zeros((len1, 2))
        cntr = 0
        for iii in list_ii:
            indices[cntr, :] = iii, list_jj[0]
            cntr += 1
        for jjj in list_jj:
            indices[cntr, :] = list_ii[-1], jjj
            cntr += 1
        count_land = 0
        k = 0
        while count_land <= LAND_THRESHOLD and k < len1:
            if land[int(indices[k, 0]), int(indices[k, 1])]:
                count_land += 1
            k += 1
        if count_land > LAND_THRESHOLD:
            ocean_weight[idx] = 0
    return ocean_weight


def spreading_weights(lon, lat, land, cells, strategy='diffuse', cleaning=True, rradius=None,
                      alongshore=None, sigma=None):
    """Sparse (jdm*idm, n_sources) float32 matrix of spreading weights; column k sums to 1.

    lon, lat, land: (jdm, idm) arrays of the model grid (land: True on land).
    cells: flat index (j*idm + i) of the ocean cell of each source.
    strategy: a key of SPREADING; rradius [m], alongshore [m] and sigma [cells] override its values.
    cleaning: apply step 4 (cleaning_discharge).
    The discharge per cell (m3/s) is weights @ discharge_per_source.
    """
    land = np.asarray(land, dtype=bool)
    jdm, idm = land.shape
    n = jdm * idm
    cells = np.asarray(cells, dtype=np.int64)
    unique_cells, source_col = np.unique(cells, return_inverse=True)

    params = SPREADING[strategy]
    if params is None and rradius is None:
        w = sp.csc_matrix((np.ones(len(unique_cells), np.float32),
                           (unique_cells, np.arange(len(unique_cells)))), shape=(n, len(unique_cells)))
        return w[:, source_col].tocsc()
    params = dict(params or {})
    params.update({k: v for k, v in dict(rradius=rradius, alongshore=alongshore, sigma=sigma).items()
                   if v is not None})
    rradius, alongshore, sigma = params['rradius'], params['alongshore'], params['sigma']

    t0 = time.time()
    land_flat = land.ravel()
    xyz = ecef(lat, lon)                                        # km
    tree = cKDTree(xyz)
    d_land, nn_land = cKDTree(xyz[land_flat]).query(xyz)
    d_land = d_land * 1000.0                                    # m
    d_land[land_flat] = 0.0
    landmass, _ = label(gaussian_filter(land.astype(float), sigma=0.1) > 0.5, structure=np.ones((3, 3), int))
    nearest_landmass = landmass.ravel()[land_flat][nn_land]     # landmass nearest to each cell

    rows, cols, vals = [], [], []
    for k, c in enumerate(unique_cells):
        if k % 500 == 0:
            print(f'  spreading weights: source cell {k}/{len(unique_cells)}', flush=True)
        cand = np.asarray(tree.query_ball_point(xyz[c], alongshore / 1000.0), dtype=np.int64)
        cand = cand[~land_flat[cand] & (d_land[cand] <= rradius) & (nearest_landmass[cand] == nearest_landmass[c])]
        w = np.zeros(n)
        dist = np.linalg.norm(xyz[cand] - xyz[c], axis=1) * 1000.0
        w[cand] = np.exp(-dist / alongshore) * np.exp(-2.0 * d_land[cand] / rradius)
        j, i = divmod(int(c), idm)
        if cleaning:
            w = cleaning_discharge(w, np.nonzero(w), j, i, idm, jdm, land)
        if sigma:
            w = gaussian_filter(w.reshape(jdm, idm), sigma=sigma).ravel()
        w[land_flat] = 0.0
        if w.sum() == 0.0:          # nothing left: keep the water in the source cell
            w[c] = 1.0
        nz = np.flatnonzero(w)
        rows.append(nz)
        cols.append(np.full(nz.size, k))
        vals.append((w[nz] / w[nz].sum()).astype(np.float32))
    w = sp.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(n, len(unique_cells)), dtype=np.float32)
    print(f'  spreading weights (rradius {rradius / 1000:.0f} km, alongshore {alongshore / 1000:.0f} km, '
          f'sigma {sigma}) for {len(unique_cells)} cells in {time.time() - t0:.0f} s')
    return w[:, source_col].tocsc()

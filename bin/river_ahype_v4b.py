#!/usr/bin/env python3
"""
HYCOM river forcing from AHYPE v4b, E-HYPE and GEUS Greenland solid ice discharge.

Monthly freshwater flux (m/s) on the model grid, the sum of three sets of point sources:
  - AHYPE v4b: 751 river-to-ocean outlets, monthly discharge 2002-2022 (time-varying)
  - E-HYPE: European outlets of rivers_ahype-ehype_clim_rev2.dat (monthly climatology)
  - GEUS: Greenland solid ice discharge per flux gate (time-varying, monthly means from
    geus_preprocess.py), added as liquid freshwater
Each outlet or gate goes to the nearest coastal ocean cell (ocean cell with a land cell as
direct neighbour); outlets outside the model domain are dropped and listed. There is no distance
limit, because a large river can enter an estuary the grid does not resolve (the Yenisei is
330 km from the nearest TP2a0.10 coastal cell). All outlets are spread in the same way:
  pointsource  no spreading
  coastal      rradius=50 km, alongshore=100 km, Gaussian sigma=1 cell
  diffuse      rradius=180 km, alongshore=360 km, sigma=3 cells (the GloFAS radii)
The spreading follows propagate_river_plume() of the GloFAS rivers: exponential decay with
distance from the mouth and from the coast, same landmass only, no crossing of more than 15
land cells, optional Gaussian smoothing. The weights of each outlet are normalised, so the
volume is conserved.

The output has one record per month at 00:00 on the 15th, for whole years, plus one padding
record on 15 December before the first year and one on 15 January after the last year. HYCOM
interpolates linearly between records, so with the padding a run can start on 1 January of the
first year and end on 31 December of the last year. A padding record holds the data of its own
month if all sources have it (e.g. December 2004 for a file starting in 2005), otherwise a copy
of the first or last month (outside 2002-2022). Normally called by river_ahype_v4b.sh or
submit_river_ahype_v4b.sh, which write force/rivers/<E>/riverh.[ab].

Usage:
  python river_ahype_v4b.py --spreading coastal --start-year 2002 --end-year 2022 \\
      --grid regional.grid.a --depth depth_TP2a0.10_04.a --out /path/riverh_ahype_v4b
"""
import argparse
import calendar
import datetime
import os
import re
import sqlite3
import sys
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp
import xarray as xr
from scipy.ndimage import binary_dilation, gaussian_filter
from scipy.ndimage import label as ndlabel
from scipy.spatial import cKDTree

_BINDIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_BINDIR, 'river-topaz', 'scripts', 'GloFAS'))
from utils.grid_operations import transform_coordinates, cleaning_discharge  # noqa: E402
import abfile  # noqa: E402

# ── Default data paths (Betzy) ───────────────────────────────────────────────
_AHYPE_V4B = '/cluster/projects/nn2993k/AHYPE_v4b'
DEFAULTS = dict(
    gpkg = f'{_AHYPE_V4B}/arcfresh_riverdischarge_ahype_v4b_river-to-ocean-outlets.gpkg',
    ts   = f'{_AHYPE_V4B}/arcfresh_riverdischarge_2002-2022_monthly_m3s_river-to-ocean-outlets.txt',
    dat  = os.path.join(_BINDIR, '..', 'input', 'rivers_ahype-ehype_clim_rev2.dat'),
    geus = f'{_AHYPE_V4B}/geus_gates_monthly.nc',
)
AHYPE_TABLE = 'arcfresh_riverdischarge_ahype_v4b_river-to-ocean-outlets'

SPREADING_PARAMS = {
    'pointsource': None,
    'coastal': dict(rradius=50_000,  alongshore=100_000, sigma=1),
    'diffuse': dict(rradius=180_000, alongshore=360_000, sigma=3),
}
R_EARTH = 6371.0
MODEL_DAY0 = datetime.date(1900, 12, 31)


# ── Grid ──────────────────────────────────────────────────────────────────────

class Grid:
    """Model grid from regional.grid.[ab] and the land mask of a depth .a file."""

    def __init__(self, grid_path, depth_path):
        bfile = open(grid_path[:-2] + '.b').read()
        self.idm = int(re.search(r'(\d+)\s+\'idm', bfile).group(1))
        self.jdm = int(re.search(r'(\d+)\s+\'jdm', bfile).group(1))
        self.n = self.idm * self.jdm
        self.npad = -(-self.n // 4096) * 4096
        if os.path.getsize(depth_path) % (self.npad * 4):
            sys.exit(f'Error: {depth_path} is not on the {self.idm}x{self.jdm} grid of {grid_path}')
        self.plon = self.read(grid_path, 0).data
        self.plat = self.read(grid_path, 1).data
        self.area = (self.read(grid_path, 9) * self.read(grid_path, 10)).data   # scpx * scpy
        self.land = self.read(depth_path, 0).mask

    def read(self, path, rec):
        with open(path, 'rb') as f:
            f.seek(rec * self.npad * 4)
            a = np.fromfile(f, '>f4', self.n).reshape(self.jdm, self.idm)
        return np.ma.masked_greater(a, 1e29)

    def ecef(self):
        return transform_coordinates(np.column_stack([self.plat.ravel(), self.plon.ravel()]))


def map_outlets(grid, lons, lats):
    """Nearest coastal ocean cell (flat index), distance to it [km], and whether each outlet lies
    inside the model domain.

    The domain is the polygon of the grid rows j=0 and j=jdm-1 in polar stereographic
    coordinates (the columns i=0 and i=idm-1 are left out: in geographic space they zigzag
    through the Aleutians), as in the AHYPE v4b forcing notebook. No distance limit is used.
    """
    from matplotlib.path import Path
    from pyproj import Proj
    nps = Proj('+proj=stere +lat_0=90 +lon_0=0 +lat_ts=90 +ellps=WGS84')
    edge = np.vstack([np.column_stack(nps(grid.plon[0, :], grid.plat[0, :])),
                      np.column_stack(nps(grid.plon[-1, ::-1], grid.plat[-1, ::-1]))])
    inside = Path(edge).contains_points(np.column_stack(nps(lons, lats)))

    coastal = ~grid.land & binary_dilation(grid.land, structure=np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]]))
    coast_idx = np.flatnonzero(coastal)
    chord, nn = cKDTree(grid.ecef()[coast_idx]).query(transform_coordinates(np.column_stack([lats, lons])))
    return coast_idx[nn], 2 * R_EARTH * np.arcsin(np.clip(chord / (2 * R_EARTH), 0, 1)), inside


# ── Spreading weights ─────────────────────────────────────────────────────────

def spreading_weights(grid, flat_idx, rradius, alongshore, sigma):
    """Sparse (n cells, n outlets) matrix of normalised spreading weights."""
    jdm, idm = grid.jdm, grid.idm
    land = grid.land.ravel()
    xyz = grid.ecef()                                    # km
    tree_grid = cKDTree(xyz)
    land_tree = cKDTree(xyz[land])
    d2land, nn_land = land_tree.query(xyz)
    d_to_land = d2land * 1000.0                          # m
    d_to_land[land] = 0.0
    labels, _ = ndlabel(gaussian_filter(grid.land.astype(float), sigma=0.1) > 0.5,
                        structure=np.ones((3, 3), int))
    closest_label = labels.ravel()[land][nn_land]        # landmass nearest to each cell

    rows, cols, vals = [], [], []
    for k, c in enumerate(flat_idx):
        if k % 200 == 0:
            print(f'    outlet {k}/{len(flat_idx)}')
        cand = np.array(tree_grid.query_ball_point(xyz[c], alongshore / 1000.0), dtype=np.int64)
        cand = cand[~land[cand] & (d_to_land[cand] <= rradius) & (closest_label[cand] == closest_label[c])]
        dist = np.linalg.norm(xyz[cand] - xyz[c], axis=1) * 1000.0
        w = np.zeros(grid.n)
        w[cand] = np.exp(-dist / alongshore) * np.exp(-2.0 * d_to_land[cand] / rradius)
        j, i = divmod(int(c), idm)
        w = cleaning_discharge(w, np.nonzero(w), j, i, idm, jdm, grid.land)
        if sigma:
            w = gaussian_filter(w.reshape(jdm, idm), sigma=sigma).ravel()
        w[land] = 0.0
        if w.sum() == 0.0:          # nothing left: keep the water at the mouth
            w[c] = 1.0
        nz = np.flatnonzero(w)
        rows.append(nz); cols.append(np.full(nz.size, k)); vals.append(w[nz] / w[nz].sum())
    return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                         shape=(grid.n, len(flat_idx)))


def flux_field(grid, flat_idx, q_m3s, weights):
    """Flux [m/s] on the grid from per-outlet discharge [m³/s]."""
    if weights is None:
        vol = np.bincount(flat_idx, weights=q_m3s, minlength=grid.n)
    else:
        vol = weights @ q_m3s
    field = vol.reshape(grid.jdm, grid.idm) / grid.area
    field[grid.land] = 0.0
    return field


# ── Sources ───────────────────────────────────────────────────────────────────

def load_ahype(gpkg, ts_path, years):
    """Outlet positions and a function month -> discharge [m³/s]."""
    con = sqlite3.connect(gpkg)
    meta = pd.read_sql(f'SELECT SUBID, POURX, POURY, oceanbasin FROM "{AHYPE_TABLE}"', con)
    con.close()
    meta = meta[meta['oceanbasin'] != 'North Pacific'].reset_index(drop=True)
    ts = pd.read_csv(ts_path, sep='\t', skiprows=1, index_col=0)
    ts.columns = ts.columns.astype(int)
    ts = ts.reindex(columns=meta['SUBID']).fillna(0.0)
    missing = [f'{y}-{m:02d}' for y in years for m in range(1, 13) if f'{y}-{m:02d}' not in ts.index]
    if missing:
        sys.exit(f'Error: months not in {ts_path}: {missing[0]} ... {missing[-1]}')
    return meta['POURX'].values, meta['POURY'].values, lambda y, m: ts.loc[f'{y}-{m:02d}'].values


def load_ehype(dat_path):
    """E-HYPE outlets of the .dat file and a function month -> discharge [m³/s]."""
    lines = open(dat_path).readlines()
    lon, lat, annual_km3, fracs = [], [], [], []
    for i in range(0, len(lines) - 2, 3):
        if 'E-HYPE' not in lines[i]:
            continue
        vals = list(map(float, lines[i + 1].split()))
        annual_km3.append(vals[0]); fracs.append(vals[1:13])        # km³/yr, monthly fractions
        la, lo = map(float, lines[i + 2].split()[:2])
        lat.append(la); lon.append(lo)
    annual_km3, fracs = np.array(annual_km3), np.array(fracs)
    return (np.array(lon), np.array(lat),
            lambda y, m: annual_km3 * fracs[:, m - 1] * 1e9 / (calendar.monthrange(y, m)[1] * 86400))


def load_geus(path, years):
    """GEUS gates and a function month -> discharge [m³/s]."""
    ds = xr.open_dataset(path)
    version = (f'Dataverse {ds.attrs.get("geus_dataverse_version", "unknown")}, '
               f'product version {ds.attrs.get("geus_product_version", "unknown")}')
    print(f'GEUS ice discharge {version} ({ds.attrs.get("geus_source", "")}), {path}')
    q = ds['discharge'].transpose('time', 'gate').to_pandas()
    q.index = q.index.strftime('%Y-%m')
    missing = [f'{y}-{m:02d}' for y in years for m in range(1, 13) if f'{y}-{m:02d}' not in q.index]
    if missing:
        sys.exit(f'Error: months not in {path}: {missing[0]} ... {missing[-1]}')
    return (ds['lon'].values, ds['lat'].values, lambda y, m: q.loc[f'{y}-{m:02d}'].values,
            version)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--spreading', choices=list(SPREADING_PARAMS), default='pointsource')
    p.add_argument('--start-year', type=int, default=2002)
    p.add_argument('--end-year', type=int, default=2022)
    p.add_argument('--out', required=True, help='output file name, without .a/.b')
    p.add_argument('--grid', required=True, help='regional.grid.a')
    p.add_argument('--depth', required=True, help='depth .a file of the experiment (land mask)')
    p.add_argument('--gpkg', default=DEFAULTS['gpkg'], help='AHYPE v4b outlets (.gpkg)')
    p.add_argument('--ts', default=DEFAULTS['ts'], help='AHYPE v4b monthly discharge (.txt)')
    p.add_argument('--dat', default=DEFAULTS['dat'], help='rivers_ahype-ehype_clim_rev2.dat (E-HYPE)')
    p.add_argument('--geus', default=DEFAULTS['geus'], help='GEUS monthly gate discharge from geus_preprocess.py')
    p.add_argument('--no-ehype', action='store_true', help='leave out E-HYPE')
    p.add_argument('--no-geus', action='store_true', help='leave out GEUS solid ice discharge')
    a = p.parse_args()
    years = list(range(a.start_year, a.end_year + 1))

    grid = Grid(a.grid, a.depth)
    print(f'Grid {grid.idm}x{grid.jdm}, land mask of {a.depth}')

    sources = {'AHYPE v4b': load_ahype(a.gpkg, a.ts, years)}
    if not a.no_ehype:
        sources['E-HYPE'] = load_ehype(a.dat)
    geus_version = None
    if not a.no_geus:
        *sources['GEUS'], geus_version = load_geus(a.geus, years)

    # Map every outlet/gate to its coastal cell; drop those outside the model domain
    names, flat_idx, keep = [], [], {}
    for name, (lon, lat, fn) in sources.items():
        idx, dist, keep[name] = map_outlets(grid, lon, lat)
        names += [name] * keep[name].sum()
        flat_idx.append(idx[keep[name]])
        print(f'{name}: {keep[name].sum()} of {len(lon)} outlets inside the domain; distance to the coastal '
              f'cell median {np.median(dist[keep[name]]):.0f} km, max {dist[keep[name]].max():.0f} km')
        if (~keep[name]).any():
            q = np.mean([fn(y, m) for y in years for m in range(1, 13)], axis=0) * 365.25 * 86400 / 1e9
            print(f'  dropped (outside the domain), {q[~keep[name]].sum():.1f} km3/yr:')
            for k in np.flatnonzero(~keep[name]):
                print(f'    lon {lon[k]:8.2f} lat {lat[k]:6.2f}  {q[k]:7.2f} km3/yr')
    flat_idx = np.concatenate(flat_idx)
    names = np.array(names)

    params = SPREADING_PARAMS[a.spreading]
    weights = None
    if params:
        print(f'Spreading weights ({a.spreading}: {params}) for {len(flat_idx)} outlets ...')
        t0 = time.time()
        weights = spreading_weights(grid, flat_idx, **params)
        print(f'  done in {time.time() - t0:.0f} s, {weights.nnz} non-zero weights')

    fout = abfile.ABFileRiverForcing(a.out, 'w', mask=True, idm=grid.idm, jdm=grid.jdm,
                                     cline1=f'AHYPE v4b + E-HYPE + GEUS ice discharge, spreading={a.spreading}',
                                     cline2='river forcing (m/s)' +
                                     (f', GEUS {geus_version}' if geus_version else ''))

    def discharge(y, m):
        return np.concatenate([fn(y, m)[keep[name]] for name, (_, _, fn) in sources.items()])

    def write(field, y, m):
        fout.write_field(field, grid.land, 'rivers', float((datetime.date(y, m, 15) - MODEL_DAY0).days),
                         float(calendar.monthrange(y, m)[1]))

    def write_padding(y, m, edge):
        """Padding record on the 15th of month (y, m): its own data if every source has that
        month, otherwise a copy of the edge month."""
        try:
            q, src = discharge(y, m), (y, m)
        except KeyError:
            q, src = discharge(*edge), edge
        write(flux_field(grid, flat_idx, q, weights), y, m)
        print(f'  padding record {y}-{m:02d}-15 with the data of {src[0]}-{src[1]:02d}')

    print('\nPadding records, so that runs can cover the whole first and last year:')
    write_padding(years[0] - 1, 12, (years[0], 1))

    print('\nAnnual totals [km3/yr]; the gridded total must equal the sum of the sources '
          'and be the same for every spreading:')
    print('  year  ' + '  '.join(f'{n:>10s}' for n in sources) + '     sources    gridded')
    for y in years:
        tot_src = dict.fromkeys(sources, 0.0)
        tot_grid = 0.0
        for m in range(1, 13):
            sec = calendar.monthrange(y, m)[1] * 86400
            q = discharge(y, m)
            for name in sources:
                tot_src[name] += q[names == name].sum() * sec / 1e9
            field = flux_field(grid, flat_idx, q, weights)
            write(field, y, m)
            tot_grid += (field * grid.area).sum() * sec / 1e9
        s = sum(tot_src.values())
        print(f'  {y}  ' + '  '.join(f'{v:10.1f}' for v in tot_src.values()) + f'  {s:10.1f} {tot_grid:10.1f}')
        if abs(tot_grid - s) > 1e-3 * s:
            sys.exit(f'Error: gridded total {tot_grid:.1f} differs from the sources {s:.1f} km3/yr')

    print()
    write_padding(years[-1] + 1, 1, (years[-1], 12))
    fout.close()
    print(f'\nWrote {a.out}.[ab]')


if __name__ == '__main__':
    main()

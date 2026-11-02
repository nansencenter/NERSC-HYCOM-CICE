#!/usr/bin/env python3
"""
Prepare HYCOM TP2 river forcing from AHYPE v4b + E-HYPE + Roshin Greenland ice melt.

Generates monthly time-varying riverh.[ab] files (12 records per year) by combining:
  - AHYPE v4b point-source monthly discharge (2002-2022, time-varying)
  - E-HYPE climatological discharge from rivers_ahype-ehype_clim_rev2.dat
  - Roshin Greenland ice melt climatology (grivers.[ab] on the TP2 grid)

Three spreading strategies are supported:
  pointsource  No spreading; each outlet injects into its single nearest coastal cell.
  coastal      Moderate spreading: rradius=50 km, alongshore=100 km, Gaussian sigma=1 cell.
  diffuse      Full operational spreading: rradius=180 km, alongshore=360 km, sigma=3 cells.

The E-HYPE and Roshin components are climatological and repeat every year unchanged.
All three components are summed into one combined output file per year.

Usage (from experiment directory or any directory with access to the data):
  python <BINDIR>/river_ahype_v4b.py --spreading pointsource --year 2010 --outdir /path/to/out/

Output: <outdir>/riverh_ahype_v4b_<spreading>_<YYYY>.[ab]

Required data (all have defaults for TP2 on Betzy/Fram):
  --grid   regional.grid.[ab]      TP2 grid file
  --depth  depth_TP2a0.10_04.[ab]  TP2 bathymetry
  --gpkg   AHYPE v4b .gpkg         outlet metadata
  --ts     AHYPE v4b timeseries    monthly discharge 2002-2022 (tab-separated, m3/s)
  --dat    rivers_ahype-ehype_clim_rev2.dat  for E-HYPE outlets
  --roshin grivers.[ab]            Roshin Greenland ice melt on TP2 grid
"""
import argparse
import copy
import datetime
import os
import sqlite3
import sys
import time

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.ndimage import label as ndlabel
from scipy.spatial import cKDTree

# ── GloFAS utilities (transform_coordinates, cleaning_discharge) ──────────────
_BINDIR = os.path.dirname(os.path.abspath(__file__))
_GLOFAS = os.path.join(_BINDIR, 'river-topaz', 'scripts', 'GloFAS')
if _GLOFAS not in sys.path:
    sys.path.insert(0, _GLOFAS)
from utils.grid_operations import transform_coordinates, cleaning_discharge  # noqa: E402
import abfile  # noqa: E402

# ── Default data paths (TP2 on Betzy/Fram) ───────────────────────────────────
_NN2993 = '/cluster/projects/nn2993k'
DEFAULTS = dict(
    grid   = f'{_NN2993}/nlo043/TP2a0.10/topo/regional.grid.a',
    depth  = f'{_NN2993}/nlo043/TP2a0.10/topo/depth_TP2a0.10_04.a',
    gpkg   = (f'{_NN2993}/nlo043/rivers_variants/rivers_variants/data/ahype/'
              'arcfresh_riverdischarge_ahype_v4b_river-to-ocean-outlets.gpkg'),
    ts     = (f'{_NN2993}/nlo043/rivers_variants/rivers_variants/data/ahype/'
              'arcfresh_riverdischarge_2002-2022_monthly_m3s_river-to-ocean-outlets.txt'),
    dat    = f'{_NN2993}/ModelInput/AHYPE_all/rivers_ahype-ehype_clim_rev2.dat',
    roshin = f'{_NN2993}/ModelInput/TRIP/triver_Roshin/grivers.a',
)

# ── Grid dimensions ──────────────────────────────────────────────────────────
IDM, JDM = 400, 380
N        = IDM * JDM
NPAD     = -(-N // 4096) * 4096

# ── Spreading parameters ─────────────────────────────────────────────────────
SPREADING_PARAMS = {
    'pointsource': dict(rradius=None, alongshore=None, sigma=None),
    'coastal':     dict(rradius=50_000,  alongshore=100_000, sigma=1),
    'diffuse':     dict(rradius=180_000, alongshore=360_000, sigma=3),
}

# ── Model time epoch ─────────────────────────────────────────────────────────
MODEL_DAY0 = datetime.date(1900, 12, 31)

AHYPE_TABLE = 'arcfresh_riverdischarge_ahype_v4b_river-to-ocean-outlets'
MAX_DIST_KM = 300.0
R_EARTH     = 6371.0


# ── Grid I/O ──────────────────────────────────────────────────────────────────

def _read_a(path, rec):
    with open(path, 'rb') as f:
        f.seek(rec * NPAD * 4)
        a = np.fromfile(f, '>f4', N).reshape(JDM, IDM)
    return np.ma.masked_greater(a, 1e29)


def load_grid(grid_path, depth_path):
    grid_fields = ['plon', 'plat', 'qlon', 'qlat', 'ulon', 'ulat',
                   'vlon', 'vlat', 'pang', 'scpx', 'scpy']
    g = {name: _read_a(grid_path, i) for i, name in enumerate(grid_fields)}
    depth    = _read_a(depth_path, 0)
    land_2d  = depth.mask
    area_2d  = (g['scpx'] * g['scpy']).data
    plon_2d  = g['plon'].data
    plat_2d  = g['plat'].data
    return plon_2d, plat_2d, area_2d, land_2d


# ── Coastal KD-tree ───────────────────────────────────────────────────────────

def build_coastal_tree(plon_2d, plat_2d, land_2d):
    from scipy.ndimage import binary_dilation
    struct      = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
    land_dilated = binary_dilation(land_2d, structure=struct)
    coastal_2d  = (~land_2d) & land_dilated
    coast_idx   = np.where(coastal_2d.flatten())[0]
    coast_lon   = plon_2d.flatten()[coast_idx]
    coast_lat   = plat_2d.flatten()[coast_idx]
    coast_ecef  = transform_coordinates(np.column_stack([coast_lat, coast_lon]))
    tree        = cKDTree(coast_ecef)
    return tree, coast_idx, coastal_2d


def map_outlets(lons, lats, tree, coast_idx, max_dist_km=MAX_DIST_KM):
    """Return (flat_idx, dist_km, in_domain) for each outlet."""
    ecef        = transform_coordinates(np.column_stack([lats, lons]))
    chord, nn   = tree.query(ecef)
    dist_km     = 2 * R_EARTH * np.arcsin(np.clip(chord / (2 * R_EARTH), 0, 1))
    flat_idx    = coast_idx[nn]
    in_domain   = dist_km <= max_dist_km
    return flat_idx, dist_km, in_domain


# ── Spreading weights ─────────────────────────────────────────────────────────

def compute_weights_coastal_diffuse(flat_idx_arr, plon_2d, plat_2d, land_2d,
                                    rradius, alongshore, sigma):
    """
    Return (JDM*IDM, n_outlets) normalised weight array for coastal or diffuse spreading.
    Uses the GloFAS algorithm: exponential decay in offshore and along-shore distance,
    same-landmass check, land-crossing removal via cleaning_discharge, optional Gaussian blur.
    """
    n          = len(flat_idx_arr)
    plon_flat  = plon_2d.flatten()
    plat_flat  = plat_2d.flatten()
    land_flat  = land_2d.flatten()
    ocean_flat = ~land_flat

    print('  Precomputing ECEF coordinates and land labels ...')
    grid_ecef = transform_coordinates(np.column_stack([plat_flat, plon_flat]))
    X, Y, Z   = grid_ecef[:, 0], grid_ecef[:, 1], grid_ecef[:, 2]

    land_ecef      = transform_coordinates(np.column_stack(
        [plat_flat[land_flat], plon_flat[land_flat]]))
    tree_land      = cKDTree(land_ecef)
    d2land, nn_land = tree_land.query(grid_ecef)
    d_to_land       = d2land * 1000.0
    d_to_land[land_flat] = 0.0

    blurred  = gaussian_filter(land_2d.astype(float), sigma=0.1)
    labeled, _ = ndlabel(blurred > 0.5, structure=np.ones((3, 3), int))
    land_labels    = labeled.flatten()[land_flat]
    closest_label  = land_labels[nn_land].reshape(JDM, IDM)

    weights = np.zeros((N, n), dtype=np.float64)

    for k, flat_k in enumerate(flat_idx_arr):
        if k % 100 == 0:
            print(f'    outlet {k}/{n} ...')
        j_k, i_k = divmod(int(flat_k), IDM)
        lon_k = plon_flat[flat_k]
        lat_k = plat_flat[flat_k]
        e_k   = transform_coordinates([[lat_k, lon_k]])[0]
        dist  = np.sqrt((e_k[0] - X)**2 + (e_k[1] - Y)**2 + (e_k[2] - Z)**2) * 1000.0

        w = np.zeros(N)
        w[ocean_flat] = (np.exp(-dist[ocean_flat] / alongshore) *
                         np.exp(-2.0 * d_to_land[ocean_flat] / rradius))
        w[dist > alongshore]    = 0.0
        w[d_to_land > rradius]  = 0.0

        lbl_k = closest_label[j_k, i_k]
        w[closest_label.flatten() != lbl_k] = 0.0

        wm = np.where(w != 0)
        w  = cleaning_discharge(copy.deepcopy(w), wm, j_k, i_k, IDM, JDM,
                                land_2d.astype(bool))
        weights[:, k] = w

    print('  Smoothing and normalising ...')
    if sigma and sigma > 0:
        weights = gaussian_filter(weights.reshape(JDM, IDM, n),
                                  sigma=(sigma, sigma, 0))
        weights = weights.reshape(N, n)
    else:
        weights = weights.reshape(N, n).copy()

    weights[land_flat, :] = 0.0
    col_sums = weights.sum(axis=0)
    col_sums = np.where(col_sums == 0, 1.0, col_sums)
    weights /= col_sums[np.newaxis, :]
    return weights   # (N, n_outlets)


def build_flux_field(flat_idx_arr, q_m3s_arr, area_2d, land_2d, weights=None):
    """
    Return (JDM, IDM) flux field in m/s given per-outlet discharge in m³/s.
    If weights is None: point-source injection.
    If weights is provided: spread discharge using pre-computed weights (N, n_outlets).
    """
    flux = np.zeros(N, dtype=np.float64)
    if weights is None:
        np.add.at(flux, flat_idx_arr, q_m3s_arr)
    else:
        flux = (weights @ q_m3s_arr)   # (N,)

    flux_2d = flux.reshape(JDM, IDM) / area_2d
    flux_2d[land_2d] = 0.0
    return flux_2d   # m/s


# ── E-HYPE climatological patch from .dat ────────────────────────────────────

def load_ehype_dat(dat_path):
    """Parse E-HYPE outlets from the .dat file. Returns DataFrame with
    columns: subid, annual_km3, fracs (array), lat, lon."""
    records = []
    with open(dat_path) as f:
        lines = f.readlines()
    for i in range(0, len(lines) - 2, 3):
        desc = lines[i].strip()
        if 'E-HYPE' not in desc:
            continue
        subid      = int(desc.split('SUBID')[-1])
        vals       = list(map(float, lines[i + 1].split()))
        annual_km3 = vals[0]                      # km³/yr directly
        fracs      = np.array(vals[1:13])          # monthly fractions, sum ≈ 1
        coord      = list(map(float, lines[i + 2].split()))
        records.append(dict(subid=subid, annual_km3=annual_km3,
                            fracs=fracs, lat=coord[0], lon=coord[1]))
    return pd.DataFrame(records)


def ehype_monthly_m3s(ehype_df, month_idx):
    """Return per-outlet discharge (m³/s) for month_idx (0-based Jan=0).
    Converts: annual_km3 * frac_m * 1e9 / (days_m * 86400)."""
    days = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    q    = (ehype_df['annual_km3'].values *
            np.array([row['fracs'][month_idx] for _, row in ehype_df.iterrows()]) *
            1e9 / (days[month_idx] * 86400))
    return q


# ── Roshin Greenland ice melt ─────────────────────────────────────────────────

def load_roshin(roshin_path, idm=IDM, jdm=JDM):
    """Read grivers.[ab] (12 monthly records) on the TP2 grid. Returns (12, JDM, IDM) m/s."""
    n    = idm * jdm
    npad = -(-n // 4096) * 4096
    months = []
    with open(roshin_path, 'rb') as f:
        for m in range(12):
            f.seek(m * npad * 4)
            a = np.fromfile(f, '>f4', n).reshape(jdm, idm)
            months.append(np.ma.masked_greater(a, 1e29).filled(0.0))
    return np.array(months)   # (12, JDM, IDM), m/s


# ── Output: ABFileRiverForcing ────────────────────────────────────────────────

def model_day(date):
    """Days since 1900-12-31 for a datetime.date."""
    return float((date - MODEL_DAY0).days)


def write_riverh(outpath, fields_ms, year, land_2d, cline1='', cline2=''):
    """Write 12 monthly records to <outpath>.[ab] using ABFileRiverForcing.

    fields_ms : (12, JDM, IDM) array in m/s
    Records are placed at the 15th of each month; rdtime = days in month.
    """
    days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    mask = land_2d.astype(np.float32)
    mask[land_2d]  = 1e30
    mask[~land_2d] = 0.0

    fout = abfile.ABFileRiverForcing(outpath, 'w',
                                     idm=IDM, jdm=JDM,
                                     cline1=cline1, cline2=cline2)
    for m in range(12):
        date_mid = datetime.date(year, m + 1, 15)
        dtime1   = model_day(date_mid)
        rdtime   = float(days_per_month[m])
        field    = fields_ms[m].astype(np.float32)
        field[land_2d] = 0.0
        fout.write_field(field, None, 'rivers', dtime1, rdtime)
    fout.close()
    print(f'  Wrote {outpath}.[ab]')


# ── Main ──────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--spreading', choices=['pointsource', 'coastal', 'diffuse'],
                   default='pointsource',
                   help='Spreading strategy (default: pointsource)')
    p.add_argument('--year', type=int, help='Single year to process')
    p.add_argument('--start-year', type=int, default=2002, help='First year (default: 2002)')
    p.add_argument('--end-year',   type=int, default=2022, help='Last year (default: 2022)')
    p.add_argument('--outdir', default='.', help='Output directory (default: .)')
    p.add_argument('--grid',   default=DEFAULTS['grid'],   help='regional.grid.a path')
    p.add_argument('--depth',  default=DEFAULTS['depth'],  help='depth.a path')
    p.add_argument('--gpkg',   default=DEFAULTS['gpkg'],   help='AHYPE v4b .gpkg path')
    p.add_argument('--ts',     default=DEFAULTS['ts'],     help='AHYPE v4b timeseries .txt')
    p.add_argument('--dat',    default=DEFAULTS['dat'],    help='EHYPE .dat path')
    p.add_argument('--roshin', default=DEFAULTS['roshin'], help='Roshin grivers.a path')
    p.add_argument('--no-ehype',  action='store_true', help='Omit E-HYPE patch')
    p.add_argument('--no-roshin', action='store_true', help='Omit Roshin Greenland ice melt')
    return p.parse_args()


def main():
    args   = parse_args()
    spread = args.spreading
    params = SPREADING_PARAMS[spread]

    years  = [args.year] if args.year else list(range(args.start_year, args.end_year + 1))
    os.makedirs(args.outdir, exist_ok=True)

    # ── Load grid ─────────────────────────────────────────────────────────────
    print('Loading grid ...')
    plon_2d, plat_2d, area_2d, land_2d = load_grid(args.grid, args.depth)

    # ── Build coastal KD-tree ──────────────────────────────────────────────────
    print('Building coastal KD-tree ...')
    tree, coast_idx, coastal_2d = build_coastal_tree(plon_2d, plat_2d, land_2d)

    # ── Load AHYPE v4b outlet metadata ─────────────────────────────────────────
    print('Loading AHYPE v4b outlets ...')
    con  = sqlite3.connect(args.gpkg)
    meta = pd.read_sql(
        f'SELECT fid, SUBID, POURX, POURY, oceanbasin FROM "{AHYPE_TABLE}"', con)
    con.close()

    # Remove North Pacific outlets (outside TP2 domain)
    meta = meta[meta['oceanbasin'] != 'North Pacific'].copy()

    flat_idx_v4b, dist_v4b, in_v4b = map_outlets(
        meta['POURX'].values, meta['POURY'].values, tree, coast_idx)
    meta = meta[in_v4b].copy().reset_index(drop=True)
    flat_idx_v4b = flat_idx_v4b[in_v4b]
    print(f'  AHYPE v4b: {len(meta)} outlets in domain')

    # ── Load AHYPE v4b discharge timeseries ────────────────────────────────────
    print('Loading AHYPE v4b timeseries ...')
    ts = pd.read_csv(args.ts, sep='\t', skiprows=1, index_col=0)
    ts.columns = ts.columns.astype(int)
    ts = ts[ts.columns.intersection(meta['SUBID'])]
    # Reorder columns to match meta
    ts = ts[[sid for sid in meta['SUBID'] if sid in ts.columns]]

    # ── Load E-HYPE outlets from .dat ──────────────────────────────────────────
    ehype_df = ehype_flat_idx = None
    if not args.no_ehype:
        print('Loading E-HYPE outlets from .dat ...')
        ehype_df = load_ehype_dat(args.dat)
        flat_idx_eh, dist_eh, in_eh = map_outlets(
            ehype_df['lon'].values, ehype_df['lat'].values, tree, coast_idx)
        ehype_df = ehype_df[in_eh].copy().reset_index(drop=True)
        ehype_flat_idx = flat_idx_eh[in_eh]
        print(f'  E-HYPE: {len(ehype_df)} outlets in domain '
              f'({ehype_df["annual_km3"].sum():.0f} km³/yr)')

    # ── Load Roshin Greenland ice melt ─────────────────────────────────────────
    roshin_months = None
    if not args.no_roshin:
        print('Loading Roshin Greenland ice melt ...')
        roshin_months = load_roshin(args.roshin)  # (12, JDM, IDM), m/s
        roshin_ann = float(
            (roshin_months.mean(0) * area_2d).sum() * 86400 * 365.25 / 1e9)
        print(f'  Roshin: {roshin_ann:.0f} km³/yr')

    # ── Compute spreading weights (once, reused for all years) ────────────────
    weights_v4b = weights_eh = None
    if spread != 'pointsource':
        r    = params['rradius']
        a    = params['alongshore']
        sig  = params['sigma']
        print(f'Computing {spread} spreading weights for AHYPE v4b '
              f'({len(meta)} outlets, r={r//1000} km, along={a//1000} km, σ={sig}) ...')
        t0 = time.time()
        weights_v4b = compute_weights_coastal_diffuse(
            flat_idx_v4b, plon_2d, plat_2d, land_2d, r, a, sig)
        print(f'  Done in {time.time()-t0:.0f} s')

        if ehype_df is not None:
            print(f'Computing {spread} spreading weights for E-HYPE ({len(ehype_df)} outlets) ...')
            t0 = time.time()
            weights_eh = compute_weights_coastal_diffuse(
                ehype_flat_idx, plon_2d, plat_2d, land_2d, r, a, sig)
            print(f'  Done in {time.time()-t0:.0f} s')

    # ── Process each year ─────────────────────────────────────────────────────
    for year in years:
        print(f'\nProcessing {year} ...')
        fields = np.zeros((12, JDM, IDM), dtype=np.float64)

        for m in range(12):
            month_str = f'{year}-{m+1:02d}'
            if month_str not in ts.index:
                print(f'  WARNING: {month_str} not in timeseries, skipping')
                continue

            # AHYPE v4b
            q_v4b = ts.loc[month_str].reindex(meta['SUBID']).fillna(0.0).values  # m³/s
            fields[m] += build_flux_field(flat_idx_v4b, q_v4b, area_2d, land_2d,
                                          weights=weights_v4b)

            # E-HYPE climatological patch
            if ehype_df is not None:
                q_eh = ehype_monthly_m3s(ehype_df, m)
                fields[m] += build_flux_field(ehype_flat_idx, q_eh, area_2d, land_2d,
                                              weights=weights_eh)

            # Roshin Greenland ice melt
            if roshin_months is not None:
                fields[m] += roshin_months[m]

        # Write output
        outname = os.path.join(args.outdir,
                               f'riverh_ahype_v4b_{spread}_{year}')
        ann_km3 = float((fields.mean(0) * area_2d).sum() * 86400 * 365.25 / 1e9)
        cline1  = f'AHYPE v4b + E-HYPE + Roshin, spreading={spread}, year={year}'
        cline2  = f'Annual total: {ann_km3:.0f} km3/yr'
        write_riverh(outname, fields, year, land_2d, cline1=cline1, cline2=cline2)
        print(f'  Annual total: {ann_km3:.0f} km³/yr')

    print('\nDone.')


if __name__ == '__main__':
    main()

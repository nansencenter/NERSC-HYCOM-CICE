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
All three components are summed into one file covering the whole period, one record per
month at 00:00 on the 15th. Normally called by river_ahype_v4b.sh or
submit_river_ahype_v4b.sh, which write the result to force/rivers/<E>/riverh.[ab].

Usage:
  python river_ahype_v4b.py --spreading coastal --start-year 2002 --end-year 2022 \
      --grid regional.grid.a --depth depth_TP2a0.10_04.a --out /path/riverh_ahype_v4b

Data inputs (defaults for TP2 on nn2993k):
  --gpkg   AHYPE v4b .gpkg         outlet metadata
  --ts     AHYPE v4b timeseries    monthly discharge 2002-2022 (tab-separated, m3/s)
  --dat    rivers_ahype-ehype_clim_rev2.dat  for E-HYPE outlets
  --roshin Roshin Greenland ice melt (12 monthly records) on the model grid
"""
import argparse
import calendar
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

# ── Default data paths (Betzy) ───────────────────────────────────────────────
_AHYPE_V4B = '/cluster/projects/nn2993k/AHYPE_v4b'
DEFAULTS = dict(
    gpkg   = f'{_AHYPE_V4B}/arcfresh_riverdischarge_ahype_v4b_river-to-ocean-outlets.gpkg',
    ts     = f'{_AHYPE_V4B}/arcfresh_riverdischarge_2002-2022_monthly_m3s_river-to-ocean-outlets.txt',
    dat    = os.path.join(_BINDIR, '..', 'input', 'rivers_ahype-ehype_clim_rev2.dat'),
    roshin = '/cluster/projects/nn2993k/TRIP/triver_Roshin',
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


def check_grid_size(path_a):
    """Stop unless the i/jdm line of the .b file next to path_a is IDM x JDM."""
    import re
    with open(path_a[:-2] + '.b') as f:
        m = re.search(r'i/jdm\s*=\s*(\d+)\s+(\d+)', f.read())
    if not m or (int(m.group(1)), int(m.group(2))) != (IDM, JDM):
        sys.exit(f'Error: {path_a} is not on the {IDM}x{JDM} grid')


def load_grid(grid_path, depth_path):
    check_grid_size(depth_path)
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


def ehype_monthly_m3s(ehype_df, year, month_idx):
    """Per-outlet discharge (m³/s) for month_idx (0-based): annual_km3 * frac_m * 1e9 / month length."""
    fracs = np.stack(ehype_df['fracs'].values)[:, month_idx]
    seconds = days_in_month(year, month_idx + 1) * 86400
    return ehype_df['annual_km3'].values * fracs * 1e9 / seconds


# ── Roshin Greenland ice melt ─────────────────────────────────────────────────

def load_roshin(roshin_path):
    """Read the Roshin ice melt (12 monthly records, m/s) on the model grid -> (12, JDM, IDM).

    roshin_path is a .a file, or a directory: then the one grivers*.a on the model grid is
    used (grivers_tp2.a for TP2a0.10 in $TRIP_PATH/triver_Roshin), as in
    merge_TRIP_Hype_abfile_4Greenland.py.
    """
    if os.path.isdir(roshin_path):
        import glob
        import re
        found = []
        for bfile in sorted(glob.glob(os.path.join(roshin_path, 'grivers*.b'))):
            m = re.search(r'i/jdm\s*=\s*(\d+)\s+(\d+)', open(bfile).read())
            if m and (int(m.group(1)), int(m.group(2))) == (IDM, JDM):
                found.append(bfile[:-2] + '.a')
        if len(found) != 1:
            sys.exit(f'Error: need exactly one grivers*.[ab] on the {IDM}x{JDM} grid in '
                     f'{roshin_path}, found {len(found)}')
        roshin_path = found[0]
    print(f'  {roshin_path}')
    check_grid_size(roshin_path)
    return np.array([_read_a(roshin_path, m).filled(0.0) for m in range(12)])


# ── Output: ABFileRiverForcing ────────────────────────────────────────────────

def model_day(date):
    """Days since 1900-12-31 for a datetime.date."""
    return float((date - MODEL_DAY0).days)


def days_in_month(year, month):
    return calendar.monthrange(year, month)[1]


# ── Main ──────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--spreading', choices=['pointsource', 'coastal', 'diffuse'],
                   default='pointsource',
                   help='Spreading strategy (default: pointsource)')
    p.add_argument('--start-year', type=int, default=2002, help='First year (default: 2002)')
    p.add_argument('--end-year',   type=int, default=2022, help='Last year (default: 2022)')
    p.add_argument('--out',   required=True, help='Output file name, without .a/.b')
    p.add_argument('--grid',  required=True, help='regional.grid.a path')
    p.add_argument('--depth', required=True, help='depth .a path (experiment topography)')
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

    years  = list(range(args.start_year, args.end_year + 1))

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
    missing = [f'{y}-{m:02d}' for y in years for m in range(1, 13)
               if f'{y}-{m:02d}' not in ts.index]
    if missing:
        sys.exit(f'Error: months not in {args.ts}: {missing[0]} ... {missing[-1]}')

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

    # ── Write one record per month, at 00:00 on the 15th ──────────────────────
    fout = abfile.ABFileRiverForcing(
        args.out, 'w', mask=True, idm=IDM, jdm=JDM,
        cline1=f'AHYPE v4b + E-HYPE + Roshin, spreading={spread}',
        cline2='river forcing (m/s)')
    # Annual totals do not depend on the spreading: compare them between strategies
    print('\nAnnual totals [km3/yr]:')
    for year in years:
        ann_km3 = 0.0
        for m in range(12):
            q_v4b = ts.loc[f'{year}-{m+1:02d}'].reindex(meta['SUBID']).fillna(0.0).values
            field = build_flux_field(flat_idx_v4b, q_v4b, area_2d, land_2d, weights=weights_v4b)
            if ehype_df is not None:
                q_eh = ehype_monthly_m3s(ehype_df, year, m)
                field += build_flux_field(ehype_flat_idx, q_eh, area_2d, land_2d,
                                          weights=weights_eh)
            if roshin_months is not None:
                field += roshin_months[m]
            ndays = days_in_month(year, m + 1)
            fout.write_field(field, land_2d, 'rivers',
                             model_day(datetime.date(year, m + 1, 15)), float(ndays))
            ann_km3 += float((field * area_2d).sum()) * ndays * 86400 / 1e9
        print(f'  {year}: {ann_km3:.1f}')
    fout.close()
    print(f'\nWrote {args.out}.[ab]')


if __name__ == '__main__':
    main()

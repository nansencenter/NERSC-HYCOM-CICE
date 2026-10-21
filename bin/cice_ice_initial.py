#!/usr/bin/env python
"""
Create a CICE cold-start file (ice_initial.nc) from reanalysis sea ice fields.

CICE (ice_ic='read', see read_state_var in ice_init.F90) reads:
   aice : ice concentration
   hi   : ice thickness over the ice-covered part of the cell (m)
   sst, sss : only read when sss_data_type='hycom' (ocn_data_hycom_init)

Concentration and thickness can come from different sources:
   --aice-source glorys : siconc from a GLORYS12 daily file (--glorys-file)
   --aice-source topaz  : siconc from a TOPAZ4b daily file (--topaz-file)
   --aice-source piomas : PIOMAS area (monthly means, interpolated in time)
   --hi-source   topaz  : TOPAZ4b sithick/siconc (sithick is the grid cell mean)
   --hi-source   glorys : GLORYS12 sithick (already thickness over ice-covered part)
   --hi-source   piomas : PIOMAS heff/area (heff is the grid cell mean)

Thickness is computed on the source grid where the source concentration
exceeds amin, interpolated, and filled by nearest neighbour so that every
point with ice in the chosen concentration field gets a thickness.

SST/SSS are taken from the top layer of a HYCOM archive (typically the nesting
file valid at the start date), so that they are consistent with the ocean
initial state.

Horizontal interpolation is done on points (inverse distance weighting of the k
nearest source points, distances as chords between unit vectors on the
sphere), not in index space. This avoids problems with periodic i-indices,
displaced poles (PIOMAS) and folds or poles of the target grid, so it can be
used for TP2 and TP5 alike.
"""
import argparse
import datetime
import logging
import os
from typing import Dict, List, Optional, Tuple

import abfile
import netCDF4
import numpy
import scipy.spatial

# Set up logger
_loglevel=logging.INFO
logger = logging.getLogger(__name__)
logger.setLevel(_loglevel)
formatter = logging.Formatter("%(asctime)s - %(name)10s - %(levelname)7s: %(message)s")
ch = logging.StreamHandler()
ch.setLevel(_loglevel)
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.propagate=False

_earth_radius_km = 6371.

# Search radius (km) for inverse distance weighting, per source (~ 3-4 source grid cells)
_radius_km = {"glorys": 25., "topaz": 50., "piomas": 100.}


def unit_vectors(lon: numpy.ndarray, lat: numpy.ndarray) -> numpy.ndarray:
   """ Unit vectors on the sphere.

   Args:
      lon: Longitudes in degrees.
      lat: Latitudes in degrees, same shape as lon.

   Returns:
      Array of shape lon.shape + (3,) with Cartesian coordinates on the unit sphere.
   """
   lon = numpy.radians(numpy.asarray(lon, dtype=float))
   lat = numpy.radians(numpy.asarray(lat, dtype=float))
   return numpy.stack([numpy.cos(lat)*numpy.cos(lon),
                       numpy.cos(lat)*numpy.sin(lon),
                       numpy.sin(lat)], axis=-1)


def interpolate(src_lon: numpy.ndarray, src_lat: numpy.ndarray, field: numpy.ndarray, valid: numpy.ndarray,
                dst_lon: numpy.ndarray, dst_lat: numpy.ndarray, k: int = 4,
                radius_km: Optional[float] = None) -> numpy.ndarray:
   """ Inverse distance weighting of the k nearest valid source points.

   Distances are chords between unit vectors, so the result does not depend on the
   index structure (periodicity, poles, folds) of either grid.

   Args:
      src_lon, src_lat: Source point coordinates in degrees.
      field: Source values, same shape as src_lon.
      valid: Boolean mask of source points to use.
      dst_lon, dst_lat: Destination point coordinates in degrees.
      k: Number of nearest source points.
      radius_km: Only use source points within this distance (None: no limit).

   Returns:
      Field on destination points, NaN where no valid source point is within radius_km.
   """
   valid = numpy.asarray(valid, dtype=bool)
   src_xyz = unit_vectors(src_lon[valid], src_lat[valid])
   src_val = numpy.asarray(field, dtype=float)[valid]
   dst_xyz = unit_vectors(dst_lon, dst_lat).reshape(-1, 3)
   tree = scipy.spatial.cKDTree(src_xyz)
   k = min(k, src_val.size)
   ub = numpy.inf if radius_km is None else radius_km/_earth_radius_km
   dist, ind = tree.query(dst_xyz, k=k, distance_upper_bound=ub)
   if k == 1:
      dist = dist[:, None]
      ind = ind[:, None]
   found = numpy.isfinite(dist)
   w = numpy.where(found, 1./numpy.maximum(dist, 1e-12), 0.)
   vals = numpy.where(found, src_val[numpy.minimum(ind, src_val.size-1)], 0.)
   wsum = w.sum(axis=1)
   out = numpy.full(wsum.shape, numpy.nan)
   ok = wsum > 0.
   out[ok] = (w[ok]*vals[ok]).sum(axis=1)/wsum[ok]
   return out.reshape(numpy.asarray(dst_lon).shape)


def read_regular(fname: str, varnames: List[str],
                 latmin: Optional[float] = None) -> Tuple[numpy.ndarray, numpy.ndarray, Dict[str, numpy.ndarray]]:
   """ Read 2D fields from a daily file on a regular lon/lat grid (Copernicus Marine subset).

   Args:
      fname: netCDF file with 1D longitude/latitude coordinates.
      varnames: Variables to read (first time record if there is a time dimension).
      latmin: Only read latitudes >= latmin.

   Returns:
      2D longitudes, 2D latitudes and a dict of fields (NaN where missing).
   """
   nc = netCDF4.Dataset(fname, "r")
   lon = nc.variables["longitude"][:]
   lat = nc.variables["latitude"][:]
   jj = slice(None) if latmin is None else (lat >= latmin)
   fields = {}
   for v in varnames:
      fld = nc.variables[v]
      fld = fld[0, jj, :] if fld.ndim == 3 else fld[jj, :]
      fields[v] = numpy.ma.filled(fld.astype(float), numpy.nan)
   nc.close()
   lon2d, lat2d = numpy.meshgrid(lon, lat[jj])
   return lon2d, lat2d, fields


def read_piomas(piomas_dir: str,
                date: datetime.datetime) -> Tuple[numpy.ndarray, numpy.ndarray, Dict[str, numpy.ndarray]]:
   """ Linear interpolation in time of PIOMAS monthly means (area, heff).

   Args:
      piomas_dir: Directory with PIOMAS_YYYY.nc.
      date: Time to interpolate to; must be bracketed by records of the adjacent years' files.

   Returns:
      2D longitudes, 2D latitudes and a dict with "area" and "heff" (NaN over land).
   """
   records = []
   for year in (date.year-1, date.year, date.year+1):
      fname = os.path.join(piomas_dir, "PIOMAS_%04d.nc" % year)
      if not os.path.exists(fname):
         continue
      nc = netCDF4.Dataset(fname, "r")
      times = netCDF4.num2date(nc.variables["time"][:], nc.variables["time"].units,
                               calendar=getattr(nc.variables["time"], "calendar", "standard"),
                               only_use_cftime_datetimes=False, only_use_python_datetimes=True)
      for irec, t in enumerate(times):
         records.append((t, fname, irec))
      nc.close()
   records.sort()
   before = [r for r in records if r[0] <= date]
   after = [r for r in records if r[0] > date]
   if not before or not after:
      raise ValueError("Date %s not bracketed by PIOMAS records in %s" % (date, piomas_dir))
   r0, r1 = before[-1], after[0]
   w1 = (date - r0[0]).total_seconds()/(r1[0] - r0[0]).total_seconds()
   logger.info("PIOMAS: %s (weight %.3f) and %s (weight %.3f)" % (r0[0], 1.-w1, r1[0], w1))
   fields = {}
   for w, (t, fname, irec) in ((1.-w1, r0), (w1, r1)):
      nc = netCDF4.Dataset(fname, "r")
      lon = numpy.asarray(nc.variables["longitude"][:])
      lat = numpy.asarray(nc.variables["latitude"][:])
      for v in ("area", "heff"):
         fld = numpy.ma.filled(nc.variables[v][irec, :, :].astype(float), numpy.nan)
         fields[v] = fields.get(v, 0.) + w*fld
      nc.close()
   return lon, lat, fields


def source_fields(source: str, args: argparse.Namespace,
                  latmin: float) -> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray]:
   """ Read concentration and ice thickness from one source on its own grid.

   Args:
      source: "glorys", "topaz" or "piomas".
      args: Command line arguments (file names, date).
      latmin: Minimum latitude to read for regular-grid sources.

   Returns:
      Longitudes, latitudes, concentration (clipped to [0,1], NaN over land) and thickness
      over the ice-covered part of the cell (NaN where undefined).
   """
   if source == "glorys":
      lon, lat, f = read_regular(args.glorys_file, ["siconc", "sithick"], latmin)
      aice, hi = f["siconc"], f["sithick"]
      # GLORYS siconc is missing over ice-free ocean, not only over land. Use
      # the land mask of a GLORYS physics file (zos) on the same grid to set
      # ice-free ocean to 0
      m_lon, m_lat, m = read_regular(args.glorys_mask_file, ["zos"], latmin)
      jj = numpy.isin(numpy.round(m_lat[:, 0], 4), numpy.round(lat[:, 0], 4))
      if jj.sum() != lat.shape[0] or not numpy.allclose(m_lon[jj], lon):
         raise ValueError("Grid of %s does not match %s" % (args.glorys_mask_file, args.glorys_file))
      ocean = numpy.isfinite(m["zos"][jj])
      logger.info("glorys: %d ocean points with missing siconc set to 0" % (ocean & numpy.isnan(aice)).sum())
      aice = numpy.where(ocean, numpy.nan_to_num(aice, nan=0.), numpy.nan)
   elif source == "topaz":
      lon, lat, f = read_regular(args.topaz_file, ["siconc", "sithick"], latmin)
      aice = f["siconc"]
      hi = f["sithick"]/numpy.where(aice > 0., aice, numpy.nan)
   elif source == "piomas":
      lon, lat, f = read_piomas(args.piomas_dir, args.date)
      aice = f["area"]
      hi = f["heff"]/numpy.where(aice > 0., aice, numpy.nan)
   else:
      raise ValueError("Unknown source %s" % source)
   return lon, lat, numpy.clip(aice, 0., 1.), hi


def fill_nearest(field: numpy.ndarray, valid: numpy.ndarray, lon: numpy.ndarray,
                 lat: numpy.ndarray) -> numpy.ndarray:
   """ Fill points where valid is False with the nearest valid value.

   Args:
      field: 2D field.
      valid: Boolean mask of valid points.
      lon, lat: Point coordinates in degrees.

   Returns:
      Filled copy of field.
   """
   out = numpy.array(field, dtype=float)
   bad = ~valid
   if bad.any():
      tree = scipy.spatial.cKDTree(unit_vectors(lon[valid], lat[valid]))
      d, i = tree.query(unit_vectors(lon[bad], lat[bad]))
      out[bad] = field[valid][i]
   return out


def main(args: argparse.Namespace) -> None:
   """ Build ice_initial.nc on the HYCOM/CICE grid and write it to args.outfile.

   Args:
      args: Parsed command line arguments (see the argument parser below).
   """
   amin, hmax = args.amin, args.hmax

   # Target grid and mask
   gfile = abfile.ABFileGrid(os.path.join(args.griddir, "regional.grid"), "r")
   plon = numpy.asarray(gfile.read_field("plon"))
   plat = numpy.asarray(gfile.read_field("plat"))
   cell_area = numpy.asarray(gfile.read_field("scpx"))*numpy.asarray(gfile.read_field("scpy"))
   gfile.close()
   nc = netCDF4.Dataset(args.kmt, "r")
   kmt = numpy.ma.filled(nc.variables["kmt"][:], 0.)
   nc.close()
   ocean = kmt > 0.5
   latmin = numpy.nanmin(plat) - 1.
   logger.info("Target grid %s, %d ocean points" % (str(plon.shape), ocean.sum()))

   # Concentration
   s_lon, s_lat, s_aice, s_hi = source_fields(args.aice_source, args, latmin)
   aice = interpolate(s_lon, s_lat, s_aice, numpy.isfinite(s_aice), plon, plat,
                      radius_km=_radius_km[args.aice_source])
   outside = ocean & numpy.isnan(aice)
   logger.info("aice (%s): %d ocean points without source data within range, set to 0" % (args.aice_source, outside.sum()))
   aice = numpy.where(ocean & numpy.isfinite(aice), aice, 0.)
   aice[aice < amin] = 0.

   # Thickness: from source points with ice, then nearest-neighbour fill so
   # that every point with aice>0 gets a thickness
   if args.hi_source != args.aice_source:
      s_lon, s_lat, s_aice, s_hi = source_fields(args.hi_source, args, latmin)
   s_ice = numpy.isfinite(s_hi) & (s_aice > amin)
   nclip = (s_hi[s_ice] > hmax).sum()
   s_hi = numpy.minimum(s_hi, hmax)
   logger.info("hi (%s): %d source ice points, %d clipped to %.1f m" % (args.hi_source, s_ice.sum(), nclip, hmax))
   hi = interpolate(s_lon, s_lat, s_hi, s_ice, plon, plat, radius_km=_radius_km[args.hi_source])
   need = (aice > 0.) & numpy.isnan(hi)
   if need.any():
      logger.info("hi: %d ice points outside %s ice cover, filled by nearest thickness" % (need.sum(), args.hi_source))
      hi[need] = interpolate(s_lon, s_lat, s_hi, s_ice, plon, plat, k=1)[need]
   hi = numpy.where(aice > 0., hi, 0.)

   # SST/SSS from top layer of HYCOM archive
   arc = abfile.ABFileArchv(args.archv, "r")
   sst = numpy.ma.filled(arc.read_field("temp", 1), numpy.nan).astype(float)
   sss = numpy.ma.filled(arc.read_field("salin", 1), numpy.nan).astype(float)
   arc.close()
   for name, fld in (("sst", sst), ("sss", sss)):
      bad = ocean & ~numpy.isfinite(fld)
      if bad.any():
         logger.warning("%s: %d ocean points are NaN/missing in %s, filled by nearest neighbour (lat %s)" %
                        (name, bad.sum(), args.archv, numpy.array2string(plat[bad], precision=2)))
      good = ocean & numpy.isfinite(fld)
      fld[:] = numpy.where(ocean, fill_nearest(fld, good, plon, plat), 0.)

   # Diagnostics
   area_km2 = (aice*cell_area).sum()*1e-6
   vol_km3 = (aice*hi*cell_area).sum()*1e-9
   logger.info("Ice area   : %.3f million km2" % (area_km2*1e-6))
   logger.info("Ice volume : %.3f thousand km3" % (vol_km3*1e-3))
   logger.info("Mean thickness of ice: %.2f m" % (vol_km3/max(area_km2, 1e-12)*1e3))
   # CICE only places ice where sst <= Tf+0.2 (read_state_var); relevant for sss_data_type='hycom'
   tf = -0.054*sss
   dropped = (aice > amin) & (sst > tf + 0.2)
   logger.info("Ice points with sst > Tf+0.2 (dropped by CICE if sss_data_type='hycom'): %d of %d (%.4f million km2)" %
               (dropped.sum(), (aice > amin).sum(), (aice*cell_area)[dropped].sum()*1e-12))

   # Write
   nc = netCDF4.Dataset(args.outfile, "w", format="NETCDF4")
   nc.createDimension("y", plon.shape[0])
   nc.createDimension("x", plon.shape[1])
   specs = [
      ("aice", aice, "1", "sea_ice_area_fraction", "Sea ice concentration (%s)" % args.aice_source),
      ("hi", hi, "m", "sea_ice_thickness", "Sea ice thickness over ice-covered area (%s)" % args.hi_source),
      ("sst", sst, "degC", "sea_surface_temperature", "Top layer temperature from %s" % os.path.basename(args.archv)),
      ("sss", sss, "1e-3", "sea_surface_salinity", "Top layer salinity from %s" % os.path.basename(args.archv)),
      ("lon", plon, "degrees_east", "longitude", "longitude"),
      ("lat", plat, "degrees_north", "latitude", "latitude"),
      ("kmt", kmt, "1", "", "CICE land mask (%s)" % os.path.basename(args.kmt)),
   ]
   for name, fld, units, stdname, longname in specs:
      var = nc.createVariable(name, "f8", ("y", "x"))
      var.units = units
      if stdname:
         var.standard_name = stdname
      var.long_name = longname
      var[:] = fld
   nc.title = "CICE initial state"
   nc.field_date = args.date.strftime("%Y-%m-%dT%H:%M:%S")
   nc.aice_source = args.aice_source
   nc.hi_source = args.hi_source
   nc.history = "%s: cice_ice_initial.py date=%s aice_source=%s hi_source=%s glorys_file=%s topaz_file=%s piomas_dir=%s archv=%s amin=%g hmax=%g" % (
      datetime.datetime.now().isoformat(timespec="seconds"), nc.field_date, args.aice_source, args.hi_source,
      args.glorys_file, args.topaz_file, args.piomas_dir, args.archv, amin, hmax)
   nc.close()
   logger.info("Wrote %s" % args.outfile)


if __name__ == "__main__":
   class DateTimeParseAction(argparse.Action):
      """ Parse YYYY-mm-ddTHH:MM:SS into a datetime """
      def __call__(self, parser: argparse.ArgumentParser, args: argparse.Namespace, values: str,
                   option_string: Optional[str] = None) -> None:
         setattr(args, self.dest, datetime.datetime.strptime(values, "%Y-%m-%dT%H:%M:%S"))

   parser = argparse.ArgumentParser(description="Create CICE ice_initial.nc from reanalysis sea ice fields")
   parser.add_argument("date", action=DateTimeParseAction, help="Initial time, format YYYY-mm-ddTHH:MM:SS")
   parser.add_argument("archv", help="HYCOM archive (without .a/.b) valid at date, used for sst/sss, e.g. nest/022/archv.1993_244_00")
   parser.add_argument("--kmt", required=True, help="CICE mask file, e.g. topo/kmt_TP2a0.10_04.nc")
   parser.add_argument("--griddir", default="topo", help="Directory with regional.grid.[ab] (default topo)")
   parser.add_argument("--aice-source", choices=["glorys", "topaz", "piomas"], default="glorys")
   parser.add_argument("--hi-source", choices=["topaz", "glorys", "piomas"], default="topaz")
   parser.add_argument("--glorys-file", help="GLORYS12 daily file with siconc, sithick")
   parser.add_argument("--glorys-mask-file", help="GLORYS12 physics file (with zos) on the same grid, for the land mask. " +
                       "Default: /nird/datapeak/NS9481K/MERCATOR_DATA/PHY/YYYY/MERCATOR-PHY-24-YYYY-MM-DD-12.nc")
   parser.add_argument("--topaz-file", help="TOPAZ4b daily file with siconc, sithick")
   parser.add_argument("--piomas-dir", default="/nird/datapeak/NS9481K/PIOMAS", help="Directory with PIOMAS_YYYY.nc")
   parser.add_argument("--amin", type=float, default=0.15, help="Minimum concentration (default 0.15)")
   parser.add_argument("--hmax", type=float, default=8.0, help="Maximum ice thickness in m (default 8)")
   parser.add_argument("-o", "--outfile", default="ice_initial.nc")
   args = parser.parse_args()
   for src in (args.aice_source, args.hi_source):
      if src == "glorys" and not args.glorys_file:
         parser.error("--glorys-file is required for source glorys")
      if src == "topaz" and not args.topaz_file:
         parser.error("--topaz-file is required for source topaz")
   if args.glorys_mask_file is None:
      args.glorys_mask_file = args.date.strftime(
         "/nird/datapeak/NS9481K/MERCATOR_DATA/PHY/%Y/MERCATOR-PHY-24-%Y-%m-%d-12.nc")
   main(args)

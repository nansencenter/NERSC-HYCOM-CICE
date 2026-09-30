#!/usr/bin/env python
"""
Create a HYCOM restart file from a HYCOM archive file, e.g. a GLORYS nesting
file (archv.YYYY_DDD_00.[ab]), to initialise HYCOM from that state.

Follows HYCOM's archv2restart (hycom_ALL/.../archive/src_2.2.35):
   u, v, dp, temp, saln : u-vel., v-vel., thknss, temp, salin (both time levels).
                          Archive velocities are baroclinic, as are restart velocities.
   ubavg, vbavg         : u_btrop, v_btrop (all 3 time levels)
   pbavg                : (srfhgt - montg1)/thref (all 3 time levels)
   dpmixl               : mix_dpth (both time levels)
   pbot, psikk, thkk    : from a template restart file of the same configuration

montg1 in the archive must be consistent with psikk/thkk of the template restart
(see calc_montg1.py), otherwise pbavg is wrong.

Records are written in the same order as in the template restart. Cells that are
data void in the template (land-only tiles) are written as data void, other land
cells as 0, as HYCOM does. NaN values in ocean cells of the archive are replaced
by the nearest valid ocean value (in index space).

The restart is valid at the archive time. Its name follows the HYCOM convention
restart.YYYY_DDD_HH_SSSS, taken from the archive name archv.YYYY_DDD_HH.
"""
import argparse
import logging
import os
import re
from typing import Dict, List, Optional, Tuple

import abfile
import modeltools.hycom
import numpy
import scipy.ndimage

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

thref = 1.e-3


def fill_nan(fld: numpy.ndarray, ocean: numpy.ndarray) -> Tuple[numpy.ndarray, numpy.ndarray]:
   """ Replace NaN in ocean cells by the nearest valid ocean value (in index space).

   Args:
      fld: 2D field, NaN where missing.
      ocean: 2D boolean ocean mask.

   Returns:
      Filled copy of fld (land set to 0) and boolean mask of the ocean cells that were filled.
   """
   fld = numpy.array(fld, dtype=float)
   good = ocean & numpy.isfinite(fld)
   bad = ocean & ~good
   if bad.any():
      ind = scipy.ndimage.distance_transform_edt(~good, return_distances=False, return_indices=True)
      fld[bad] = fld[ind[0][bad], ind[1][bad]]
   fld[~ocean] = 0.
   return fld, bad


def main(archv: str, template: str, blkdat: str, outdir: str, zero_velocity: bool = False,
         zero_ssh: bool = False, outname: Optional[str] = None) -> str:
   """ Create a HYCOM restart file from a HYCOM archive file.

   Args:
      archv: Archive file name without .a/.b.
      template: Restart file of the same configuration without .a/.b; provides pbot,
         psikk, thkk, the record layout and the data void (land-only tile) mask.
      blkdat: blkdat.input of the experiment (idm, jdm, kdm, iexpt, baclin, thbase).
      outdir: Output directory.
      zero_velocity: Set baroclinic and barotropic velocities to zero.
      zero_ssh: Set pbavg to zero.
      outname: Output name without .a/.b. Default restart.YYYY_DDD_HH_0000 from archive name.

   Returns:
      Output file name without .a/.b.
   """
   bp = modeltools.hycom.BlkdatParser(blkdat)
   idm, jdm, kdm = int(bp["idm"]), int(bp["jdm"]), int(bp["kdm"])
   iexpt = int(bp["iexpt"])
   baclin = float(bp["baclin"])
   thbase = float(bp["thbase"])

   # Template restart: header, record order, pbot/psikk/thkk
   tmpl = abfile.ABFileRestart(template, "r", idm=idm, jdm=jdm)
   if abs(tmpl._thbase - thbase) > 1e-6:
      raise ValueError("thbase in %s (%g) differs from blkdat.input (%g)" % (template, tmpl._thbase, thbase))
   records = [tmpl._fields[i] for i in sorted(tmpl._fields)]
   void = numpy.ma.getmaskarray(tmpl.read_field("pbot", 0))  # land-only tiles
   fixed: Dict[Tuple[str, int], numpy.ndarray] = {}
   for rec in records:
      if rec["field"] in ("pbot", "psikk", "thkk"):
         fixed[(rec["field"], rec["tlevel"])] = numpy.ma.filled(tmpl.read_field(rec["field"], rec["k"], rec["tlevel"]), 0.)

   # Archive
   arc = abfile.ABFileArchv(archv, "r")
   ocean = ~numpy.ma.getmaskarray(arc.read_field("thknss", 1))
   if (ocean & void).any():
      raise ValueError("%d ocean cells of %s are data void in template %s" % ((ocean & void).sum(), archv, template))
   logger.info("Archive %s: %d ocean cells" % (archv, ocean.sum()))

   filled: Dict[str, List[Tuple[int, numpy.ndarray]]] = {}

   def get(name: str, k: int = 0, velocity: bool = False) -> numpy.ndarray:
      """ Read archive field name at layer k (0 for 2D fields), with NaN filled.

      Velocities (u/v grids) keep their archive values with NaN set to 0; other fields
      have NaN ocean cells filled by the nearest ocean value and land set to 0.
      """
      fld = numpy.ma.filled(arc.read_field(name, k), numpy.nan)
      if velocity:
         # u/v live on their own grids, keep archive values, only remove NaN
         return numpy.nan_to_num(fld, nan=0.)
      fld, bad = fill_nan(fld, ocean)
      if bad.any():
         filled.setdefault(name, []).append((k, bad))
      return fld

   dp = numpy.array([get("thknss", k) for k in range(1, kdm+1)])
   pbot = fixed[("pbot", 1)]
   rel = numpy.abs(dp.sum(axis=0) - pbot)[ocean]/pbot[ocean]
   logger.info("sum(dp) vs template pbot: max relative difference %.3g" % rel.max())
   if rel.max() > 1e-4:
      raise ValueError("Layer thicknesses of %s do not add up to pbot of %s (max rel. diff %.3g)" % (archv, template, rel.max()))

   fields = {
      "dp": dp,
      "temp": numpy.array([get("temp", k) for k in range(1, kdm+1)]),
      "saln": numpy.array([get("salin", k) for k in range(1, kdm+1)]),
      "u": numpy.array([get("u-vel.", k, velocity=True) for k in range(1, kdm+1)]),
      "v": numpy.array([get("v-vel.", k, velocity=True) for k in range(1, kdm+1)]),
      "ubavg": get("u_btrop", velocity=True),
      "vbavg": get("v_btrop", velocity=True),
      "pbavg": (get("srfhgt") - get("montg1"))/thref,
      "dpmixl": get("mix_dpth"),
   }
   fields["pbavg"][~ocean] = 0.
   for name, lst in filled.items():
      cells = numpy.any([bad for k, bad in lst], axis=0)
      j, i = numpy.where(cells)
      logger.warning("%s: NaN in %d ocean cells, %d layers, filled by nearest neighbour; (j,i)=%s" %
                     (name, cells.sum(), len(lst), list(zip(j[:5].tolist(), i[:5].tolist()))))
   if zero_velocity:
      logger.info("Setting u, v, ubavg, vbavg to zero")
      for name in ("u", "v", "ubavg", "vbavg"):
         fields[name][:] = 0.
   if zero_ssh:
      logger.info("Setting pbavg to zero")
      fields["pbavg"][:] = 0.
   logger.info("pbavg range %.1f to %.1f Pa" % (fields["pbavg"][ocean].min(), fields["pbavg"][ocean].max()))

   # Time of archive
   dtime = None
   for d in arc._fields.values():
      dtime = float(d["day"])
      break
   arc.close()
   nstep = int(round(dtime*86400./baclin))
   if outname is None:
      m = re.match(r"archv\.([0-9]{4})_([0-9]{3})_([0-9]{2})", os.path.basename(archv))
      if not m:
         raise ValueError("Cannot derive restart name from %s, use --outname" % archv)
      outname = "restart.%s_%s_%s_0000" % m.groups()
   outfile = os.path.join(outdir, outname)
   for ext in (".a", ".b"):
      if os.path.exists(outfile + ext):
         raise ValueError("%s exists, will not overwrite" % (outfile + ext))
   logger.info("Restart time: model day %.5f, nstep %d (baclin %g s)" % (dtime, nstep, baclin))

   out = abfile.ABFileRestart(outfile, "w", mask=True, idm=idm, jdm=jdm)
   out.write_header(iexpt, tmpl._iversn, tmpl._yrflag, tmpl._sigver, nstep, dtime, thbase)
   for rec in records:
      name, k, tl = rec["field"], rec["k"], rec["tlevel"]
      if (name, tl) in fixed:
         fld = fixed[(name, tl)]
      elif name in fields:
         fld = fields[name][k-1] if fields[name].ndim == 3 else fields[name]
      else:
         raise ValueError("Field %s in template %s is not supported" % (name, template))
      fld = numpy.where(void, 0., fld)
      out.write_field(fld, void, name, k, tl)
   out.close()
   tmpl.close()
   logger.info("Wrote %s.[ab]" % outfile)
   return outfile


if __name__ == "__main__":
   parser = argparse.ArgumentParser(description="Create a HYCOM restart file from a HYCOM archive file")
   parser.add_argument("archv", help="Archive file without .a/.b, e.g. ../nest/022/archv.1993_244_00")
   parser.add_argument("template", help="Restart file of the same configuration (without .a/.b), for pbot, psikk, thkk and record layout")
   parser.add_argument("--blkdat", default="blkdat.input", help="blkdat.input of the experiment (default ./blkdat.input)")
   parser.add_argument("--outdir", default=".", help="Output directory (default .)")
   parser.add_argument("--outname", help="Output name without .a/.b (default restart.YYYY_DDD_HH_0000 from archive name)")
   parser.add_argument("--zero-velocity", action="store_true", help="Set baroclinic and barotropic velocities to zero")
   parser.add_argument("--zero-ssh", action="store_true", help="Set pbavg (barotropic pressure, i.e. SSH) to zero")
   args = parser.parse_args()
   main(args.archv, args.template, args.blkdat, args.outdir, zero_velocity=args.zero_velocity,
        zero_ssh=args.zero_ssh, outname=args.outname)

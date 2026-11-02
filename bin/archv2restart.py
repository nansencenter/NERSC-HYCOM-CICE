#!/usr/bin/env python
"""
Create a HYCOM restart file from a HYCOM archive file, e.g. a GLORYS nesting
file (archv.YYYY_DDD_00.[ab]), to initialise HYCOM from that state. No template
restart file is needed: everything comes from the archive and blkdat.input.

Fields (time levels as in a restart written by HYCOM):
   u, v, dp, temp, saln : u-vel., v-vel., thknss, temp, salin (both time levels).
                          Archive velocities are baroclinic, as are restart velocities.
   ubavg, vbavg         : u_btrop, v_btrop (all 3 time levels)
   pbot, thkk, psikk    : computed from the archive as at a HYCOM cold start (inicon.F90),
                          i.e. as if the archive state were the initial state of the run:
                             pbot  = sum of the layer thicknesses
                             montg(1) = 0,  montg(k+1) = montg(k) - p(k+1)*(thstar(k+1)-thstar(k))*thref**2
                             thkk  = thstar(kdm)  (virtual potential density in the bottom layer)
                             psikk = montg(kdm)   (Montgomery potential in the bottom layer)
   pbavg                : barotropic pressure consistent with srfhgt of the archive (all 3
                          time levels), using the decomposition of calc_montg1.py:
                             montg1 = montg1pb*pbavg + montg1c,  srfhgt = montg1 + thref*pbavg
                             => pbavg = (srfhgt - montg1c)/(montg1pb + thref)
                          With psikk/thkk from above, montg1c = 0.
   dpmixl               : mix_dpth (both time levels)

psikk and thkk are kept unchanged by HYCOM for the whole run (and by cycle_wrap.py).
montg1 in the nesting files used for the boundary forcing must therefore be computed
with this restart file as reference (calc_montg1.py); a warning is given if montg1 in
the archive is not consistent with it.

Only kapref=0 (no thermobaric correction, thstar = sigma(T,S) - thbase) is supported.
Land cells are written as 0. NaN values in ocean cells of the archive are replaced by
the nearest valid ocean value (in index space).

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

from calc_montg1 import montg1_pb, montg1_no_pb   # bin/calc_montg1.py

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
g = 9.806


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


def reference_state(thstar: numpy.ndarray, p: numpy.ndarray,
                    ocean: numpy.ndarray) -> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:
   """ pbot, thkk and psikk as computed by HYCOM at a cold start (inicon.F90).

   The initial state is taken to be at rest with a flat sea surface (montg = 0 in the
   top layer); the Montgomery potential is integrated downward through the layers.

   Args:
      thstar: Virtual potential density minus thbase, shape (kdm, jdm, idm).
      p: Layer interface pressures, shape (kdm+1, jdm, idm), p[0] = 0.
      ocean: 2D ocean mask.

   Returns:
      pbot, thkk, psikk, each of shape (jdm, idm), 0 on land.
   """
   kdm = thstar.shape[0]
   montg = numpy.zeros(thstar.shape[1:])
   for k in range(kdm-1):
      montg = montg - p[k+1]*(thstar[k+1] - thstar[k])*thref**2
   pbot = numpy.where(ocean, p[kdm], 0.)
   thkk = numpy.where(ocean, thstar[kdm-1], 0.)
   psikk = numpy.where(ocean, montg, 0.)
   return pbot, thkk, psikk


def barotropic_pressure(thstar: numpy.ndarray, p: numpy.ndarray, srfhgt: numpy.ndarray,
                        psikk: numpy.ndarray, thkk: numpy.ndarray, pbot: numpy.ndarray,
                        ocean: numpy.ndarray) -> Tuple[numpy.ndarray, numpy.ndarray]:
   """ Barotropic pressure pbavg consistent with srfhgt, psikk and thkk (as in calc_montg1.py).

   Args:
      thstar: Virtual potential density minus thbase, shape (kdm, jdm, idm).
      p: Layer interface pressures, shape (kdm+1, jdm, idm).
      srfhgt: Sea surface height times g (m2/s2), shape (jdm, idm).
      psikk: Montgomery potential in the bottom layer.
      thkk: Virtual potential density in the bottom layer.
      pbot: Bottom pressure.
      ocean: 2D ocean mask.

   Returns:
      pbavg (Pa) and the corresponding montg1 (m2/s2), 0 on land.
   """
   montg1pb = montg1_pb(thstar, p, numpy.where(ocean, pbot, 1.))
   montg1c = montg1_no_pb(psikk, thkk, thstar, p)
   pbavg = numpy.where(ocean, (srfhgt - montg1c)/(montg1pb + thref), 0.)
   montg1 = numpy.where(ocean, montg1pb*pbavg + montg1c, 0.)
   return pbavg, montg1


def restart_records(kdm: int, kapnum: int) -> List[Tuple[str, int, int]]:
   """ Records (field, layer, time level) of a HYCOM 2.3 physics restart, in file order.

   Args:
      kdm: Number of layers.
      kapnum: Number of thermobaric reference states (1 for kapref>=0).

   Returns:
      List of (field name, layer k (0 for 2D), time level).
   """
   recs: List[Tuple[str, int, int]] = []
   for name in ("u", "v", "dp", "temp", "saln"):
      recs += [(name, k, tl) for tl in (1, 2) for k in range(1, kdm+1)]
   for name in ("ubavg", "vbavg", "pbavg"):
      recs += [(name, 0, tl) for tl in (1, 2, 3)]
   recs += [("pbot", 0, 1)]
   recs += [("psikk", 0, tl) for tl in range(1, kapnum+1)]
   recs += [("thkk", 0, tl) for tl in range(1, kapnum+1)]
   recs += [("dpmixl", 0, tl) for tl in (1, 2)]
   return recs


def main(archv: str, blkdat: str, outdir: str, sigver: int, zero_velocity: bool = False,
         zero_ssh: bool = False, outname: Optional[str] = None) -> str:
   """ Create a HYCOM restart file from a HYCOM archive file.

   Args:
      archv: Archive file name without .a/.b.
      blkdat: blkdat.input of the experiment (idm, jdm, kdm, iexpt, iversn, yrflag,
         baclin, thbase, thflag, kapref).
      outdir: Output directory.
      sigver: Equation of state version (SIGVER in EXPT.src), written to the header.
      zero_velocity: Set baroclinic and barotropic velocities to zero.
      zero_ssh: Set pbavg to zero (flat sea surface).
      outname: Output name without .a/.b. Default restart.YYYY_DDD_HH_0000 from archive name.

   Returns:
      Output file name without .a/.b.
   """
   bp = modeltools.hycom.BlkdatParser(blkdat)
   idm, jdm, kdm = int(bp["idm"]), int(bp["jdm"]), int(bp["kdm"])
   iexpt, iversn, yrflag = int(bp["iexpt"]), int(bp["iversn"]), int(bp["yrflag"])
   baclin = float(bp["baclin"])
   thbase = float(bp["thbase"])
   thflag = int(bp["thflag"])
   kapref = int(bp["kapref"])
   if kapref != 0:
      raise NotImplementedError("kapref=%d in %s: only kapref=0 is supported" % (kapref, blkdat))
   kapnum = 1

   # Archive
   arc = abfile.ABFileArchv(archv, "r")
   ocean = ~numpy.ma.getmaskarray(arc.read_field("thknss", 1))
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

   fields = {
      "dp": numpy.array([get("thknss", k) for k in range(1, kdm+1)]),
      "temp": numpy.array([get("temp", k) for k in range(1, kdm+1)]),
      "saln": numpy.array([get("salin", k) for k in range(1, kdm+1)]),
      "u": numpy.array([get("u-vel.", k, velocity=True) for k in range(1, kdm+1)]),
      "v": numpy.array([get("v-vel.", k, velocity=True) for k in range(1, kdm+1)]),
      "ubavg": get("u_btrop", velocity=True),
      "vbavg": get("v_btrop", velocity=True),
      "dpmixl": get("mix_dpth"),
   }
   srfhgt = get("srfhgt")
   montg1_arc = get("montg1")
   dtime = float(next(iter(arc._fields.values()))["day"])
   arc.close()
   for name, lst in filled.items():
      cells = numpy.any([bad for k, bad in lst], axis=0)
      j, i = numpy.where(cells)
      logger.warning("%s: NaN in %d ocean cells, %d layers, filled by nearest neighbour; (j,i)=%s" %
                     (name, cells.sum(), len(lst), list(zip(j[:5].tolist(), i[:5].tolist()))))

   # Reference state (as at a HYCOM cold start) and barotropic pressure
   sig = modeltools.hycom.Sigma(thflag)
   thstar = numpy.asarray(sig.sig(fields["temp"], fields["saln"]), dtype=float) - thbase
   dp = fields["dp"]
   p = numpy.concatenate([numpy.zeros((1, jdm, idm)), numpy.cumsum(dp, axis=0)])
   fields["pbot"], fields["thkk"], fields["psikk"] = reference_state(thstar, p, ocean)
   fields["pbavg"], montg1 = barotropic_pressure(thstar, p, srfhgt, fields["psikk"], fields["thkk"],
                                                 fields["pbot"], ocean)
   logger.info("pbot range %.4g to %.4g Pa; psikk range %.4g to %.4g m2/s2; thkk range %.4g to %.4g" %
               (fields["pbot"][ocean].min(), fields["pbot"][ocean].max(), fields["psikk"][ocean].min(),
                fields["psikk"][ocean].max(), fields["thkk"][ocean].min(), fields["thkk"][ocean].max()))

   # Compare with montg1 in the archive, excluding cells where NaN were filled
   nanfilled = numpy.zeros(ocean.shape, dtype=bool)
   for lst in filled.values():
      for k, bad in lst:
         nanfilled |= bad
   cmp = ocean & ~nanfilled
   dmontg = numpy.abs(montg1 - montg1_arc)[cmp].max()
   logger.info("montg1 in archive vs montg1 implied by this restart: max difference %.3g m2/s2 (%.3g m SSH)" %
               (dmontg, dmontg/g))
   if dmontg/g > 0.001:
      logger.warning("montg1 in %s is not consistent with psikk/thkk of the new restart. The nesting files used "
                     "for the boundary forcing must be corrected with calc_montg1.py, using the new restart as "
                     "reference" % archv)

   if zero_velocity:
      logger.info("Setting u, v, ubavg, vbavg to zero")
      for name in ("u", "v", "ubavg", "vbavg"):
         fields[name][:] = 0.
   if zero_ssh:
      logger.info("Setting pbavg to zero")
      fields["pbavg"][:] = 0.
   logger.info("pbavg range %.1f to %.1f Pa" % (fields["pbavg"][ocean].min(), fields["pbavg"][ocean].max()))

   # Write
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

   out = abfile.ABFileRestart(outfile, "w", mask=False, idm=idm, jdm=jdm)
   out.write_header(iexpt, iversn, yrflag, sigver, nstep, dtime, thbase)
   for name, k, tl in restart_records(kdm, kapnum):
      fld = fields[name]
      fld = fld[k-1] if fld.ndim == 3 else fld
      out.write_field(numpy.where(ocean, fld, 0.) if name not in ("u", "v", "ubavg", "vbavg") else fld,
                      None, name, k, tl)
   out.close()
   logger.info("Wrote %s.[ab]" % outfile)
   return outfile


if __name__ == "__main__":
   parser = argparse.ArgumentParser(description="Create a HYCOM restart file from a HYCOM archive file")
   parser.add_argument("archv", help="Archive file without .a/.b, e.g. ../nest/022/archv.1993_244_00")
   parser.add_argument("--blkdat", default="blkdat.input", help="blkdat.input of the experiment (default ./blkdat.input)")
   parser.add_argument("--sigver", type=int, default=None,
                       help="Equation of state version for the header (default: SIGVER from the environment, set by EXPT.src)")
   parser.add_argument("--outdir", default=".", help="Output directory (default .)")
   parser.add_argument("--outname", help="Output name without .a/.b (default restart.YYYY_DDD_HH_0000 from archive name)")
   parser.add_argument("--zero-velocity", action="store_true", help="Set baroclinic and barotropic velocities to zero")
   parser.add_argument("--zero-ssh", action="store_true", help="Set pbavg (barotropic pressure, i.e. SSH) to zero")
   args = parser.parse_args()
   sigver = args.sigver
   if sigver is None:
      if "SIGVER" not in os.environ:
         parser.error("SIGVER not set: source EXPT.src or use --sigver")
      sigver = int(os.environ["SIGVER"])
   main(args.archv, args.blkdat, args.outdir, sigver, zero_velocity=args.zero_velocity,
        zero_ssh=args.zero_ssh, outname=args.outname)

#!/usr/bin/env python
"""
Check that the initial files of an INITFLG="--init-ice" start match the experiment.

Called by expt_preprocess.sh in the scratch directory, after the HYCOM restart, blkdat.input,
regional.depth and (with CICE) cice_kmt.nc and ice_initial.nc are in place.

Checked against the experiment, for any HYCOM restart:
   - grid size (idm, jdm): size of the .a file
   - number of layers (kdm): number of dp records
   - topography: pbot = depth*onem (HYCOM sets pbot from the depths)
   - equation of state: sigver and thbase in the restart header
   - CICE land mask: kmt in ice_initial.nc (written by cice_ice_initial.py) against cice_kmt.nc

Restarts built by archv2restart.py also have a copy of the blkdat.input they were built with
(restart.YYYY_DDD_HH_0000.blkdat). If it is there, the layers (sigma) and thflag are checked
as well; otherwise a warning is given.

Exit status 1 if anything does not match, otherwise 0.
"""
import argparse
import logging
import os
import sys
from typing import List

import abfile
import modeltools.hycom
import netCDF4
import numpy

# Set up logger
_loglevel=logging.INFO
logger = logging.getLogger("check_init_files")
logger.setLevel(_loglevel)
formatter = logging.Formatter("%(asctime)s - %(name)10s - %(levelname)7s: %(message)s")
ch = logging.StreamHandler()
ch.setLevel(_loglevel)
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.propagate=False

onem = 9806.         # HYCOM: pressure (Pa) per m of depth
_rtol_depth = 1e-4   # relative tolerance for pbot/onem vs depth
_tol_sigma = 1e-6

# Hybrid coordinate settings: differences give a warning only (HYCOM adjusts the layers)
_hybrid_keys = ["nhybrd", "nsigma", "dp00", "dp00x", "dp00f", "ds00", "ds00x", "ds00f", "isotop"]


def check_restart(restart: str, blkdat: str, depth: str, sigver: int, errors: List[str]) -> None:
   """ Check grid, layers, topography and equation of state of a HYCOM restart. """
   bp = modeltools.hycom.BlkdatParser(blkdat)
   idm, jdm, kdm = int(bp["idm"]), int(bp["jdm"]), int(bp["kdm"])

   # Grid size from the size of the .a file (records padded to multiples of 4096 words)
   nrec = sum(1 for line in open(restart + ".b") if "layer,tlevel,range" in line)
   n = idm*jdm
   reclen = 4*((n + 4095)//4096)*4096
   size = os.path.getsize(restart + ".a")
   if size != nrec*reclen:
      errors.append("%s.a has %d bytes, expected %d for %d records on idm,jdm=%d,%d (blkdat.input): "
                    "the restart is for another grid" % (restart, size, nrec*reclen, nrec, idm, jdm))
      return

   r = abfile.ABFileRestart(restart, "r", idm=idm, jdm=jdm)

   # Layers
   kdm_r = sum(1 for d in r._fields.values() if d["field"] == "dp" and d["tlevel"] == 1)
   if kdm_r != kdm:
      errors.append("restart has %d layers, blkdat.input has kdm=%d" % (kdm_r, kdm))

   # Equation of state
   if r._sigver != sigver:
      errors.append("restart has sigver=%d, experiment has SIGVER=%d (EXPT.src)" % (r._sigver, sigver))
   if abs(r._thbase - float(bp["thbase"])) > 1e-6:
      errors.append("restart has thbase=%g, blkdat.input has thbase=%g" % (r._thbase, float(bp["thbase"])))

   # Topography: pbot = depth*onem in the ocean, 0 (or NaN) on land
   d = abfile.ABFileBathy(depth, "r", idm=idm, jdm=jdm, mask=True).read_field("depth")
   d = numpy.ma.filled(d, 0.)
   ocean = d > 0.
   pbot = numpy.ma.filled(r.read_field("pbot", 0), numpy.nan)
   r.close()
   ocean_r = numpy.isfinite(pbot) & (pbot > 0.)
   nmask = (ocean != ocean_r).sum()
   if nmask > 0:
      errors.append("land mask of the restart differs from %s in %d cells: the restart is for another "
                    "topography" % (depth, nmask))
   else:
      rel = numpy.abs(pbot[ocean]/onem - d[ocean])/d[ocean]
      if rel.max() > _rtol_depth:
         errors.append("pbot of the restart differs from the depths in %s by up to %.3g m: the restart is for "
                       "another topography" % (depth, numpy.abs(pbot[ocean]/onem - d[ocean]).max()))
   logger.info("Restart %s: grid %dx%d, %d layers, sigver %d, thbase %g, topography checked against %s" %
               (restart, idm, jdm, kdm_r, r._sigver, r._thbase, depth))


def check_blkdat(restart: str, blkdat: str, errors: List[str], warnings: List[str]) -> None:
   """ Compare the layers with the blkdat.input copy written by archv2restart.py, if present. """
   ref = restart + ".blkdat"
   if not os.path.exists(ref):
      warnings.append("%s not found: cannot check that the layers (sigma) and thflag of the restart match "
                      "blkdat.input. Restarts built by archv2restart.py come with this file" % ref)
      return
   bp = modeltools.hycom.BlkdatParser(blkdat)
   bpr = modeltools.hycom.BlkdatParser(ref)
   if int(bp["thflag"]) != int(bpr["thflag"]):
      errors.append("thflag=%d in blkdat.input, but the restart was built with thflag=%d (%s)" %
                    (int(bp["thflag"]), int(bpr["thflag"]), ref))
   sig, sigr = numpy.array(bp["sigma"]), numpy.array(bpr["sigma"])
   if sig.shape != sigr.shape:
      errors.append("%d layers (sigma) in blkdat.input, but the restart was built with %d (%s)" %
                    (sig.size, sigr.size, ref))
   else:
      k = numpy.where(numpy.abs(sig - sigr) > _tol_sigma)[0]
      if k.size:
         errors.append("layer target densities (sigma) in blkdat.input differ from those the restart was built "
                       "with (%s), first in layer %d: %g vs %g" % (ref, k[0]+1, sig[k[0]], sigr[k[0]]))
   for key in _hybrid_keys:
      try:
         a, b = bp[key], bpr[key]
      except KeyError:
         continue
      if a != b:
         warnings.append("%s=%s in blkdat.input, but %s in %s (hybrid coordinate settings; HYCOM adjusts the "
                         "layers)" % (key, a, b, ref))
   logger.info("Layers (sigma) and thflag checked against %s" % ref)


def check_ice(ice_initial: str, kmt: str, errors: List[str], warnings: List[str]) -> None:
   """ Compare the CICE land mask in ice_initial.nc with the one of the experiment. """
   with netCDF4.Dataset(kmt, "r") as nc:
      k = numpy.ma.filled(nc.variables["kmt"][:], 0.) > 0.5
   with netCDF4.Dataset(ice_initial, "r") as nc:
      if "kmt" not in nc.variables:
         warnings.append("%s has no kmt field: cannot check its land mask. Files built by cice_ice_initial.py "
                         "have it" % ice_initial)
         return
      ki = numpy.ma.filled(nc.variables["kmt"][:], 0.) > 0.5
   if ki.shape != k.shape:
      errors.append("%s is on a %s grid, %s on %s" % (ice_initial, ki.shape, kmt, k.shape))
   elif (ki != k).any():
      errors.append("land mask (kmt) of %s differs from %s in %d cells: the file is for another topography" %
                    (ice_initial, kmt, (ki != k).sum()))
   else:
      logger.info("CICE land mask of %s checked against %s" % (ice_initial, kmt))


def main(restart: str, blkdat: str, depth: str, sigver: int, ice_initial: str = None, kmt: str = None) -> int:
   errors: List[str] = []
   warnings: List[str] = []
   check_restart(restart, blkdat, depth, sigver, errors)
   check_blkdat(restart, blkdat, errors, warnings)
   if ice_initial is not None and kmt is not None:
      check_ice(ice_initial, kmt, errors, warnings)
   for w in warnings:
      logger.warning(w)
   for e in errors:
      logger.error(e)
   if errors:
      logger.error("The initial files do not match the experiment (see above)")
      return 1
   logger.info("Initial files match the experiment")
   return 0


if __name__ == "__main__":
   parser = argparse.ArgumentParser(description="Check that the initial files of an --init-ice start match the experiment")
   parser.add_argument("restart", help="HYCOM restart without .a/.b, e.g. restart.1993_244_00_0000")
   parser.add_argument("--blkdat", default="blkdat.input", help="blkdat.input of the experiment (default ./blkdat.input)")
   parser.add_argument("--depth", default="regional.depth", help="Depth file without .a/.b (default regional.depth)")
   parser.add_argument("--sigver", type=int, default=None, help="SIGVER of the experiment (default: from the environment)")
   parser.add_argument("--ice-initial", help="CICE initial state, e.g. ice_initial.nc")
   parser.add_argument("--kmt", help="CICE land mask of the experiment, e.g. cice_kmt.nc")
   args = parser.parse_args()
   sigver = args.sigver
   if sigver is None:
      if "SIGVER" not in os.environ:
         parser.error("SIGVER not set: source EXPT.src or use --sigver")
      sigver = int(os.environ["SIGVER"])
   sys.exit(main(args.restart, args.blkdat, args.depth, sigver, args.ice_initial, args.kmt))

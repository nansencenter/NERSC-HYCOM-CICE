#!/usr/bin/env python
"""
Create BGC initialisation files (relax.<tracer>.[ab]) from a BGC nesting archive,
e.g. archv_fabm.1993_244_00 (CMEMS BGC reanalysis), for starting BGC together with
physics from a restart file (INITFLG="--init-ice" with ntracr>0).

HYCOM initialises the tracers when it starts from a restart without tracers
(ntracr<0 in blkdat.input, see initrc in trcupd.F90): first all FABM tracers are set
to their default values (fabm.yaml), then every tracer with a file relax.<tracer>.a is
set from that file, for the start month, remapped from the relax.intf layers to the
model layers. These files therefore have the format of the monthly tracer climatology:
12 months, each on the relax.intf interfaces of that month.

Here, all 12 months hold the same state, the BGC archive valid at the start date,
remapped conservatively from the layers of the (physics) archive to the relax.intf
layers of each month. HYCOM then remaps them back to the model layers, which at the
start are the layers of the physics archive.

The BGC archive must be on the layers of the physics archive (as written by
nemo2archvz_native.py) and valid at the same time. Values > 1e29 are treated as
missing; NaN in ocean cells are replaced by the nearest ocean value (index space).
"""
import argparse
import logging
import os
from typing import List, Tuple

import abfile
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

_missing = 1.e29  # BGC archives use 2**99 for missing values


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


def remap_layers(val: numpy.ndarray, dp_src: numpy.ndarray, dp_dst: numpy.ndarray) -> numpy.ndarray:
   """ Conservative remapping of a layered field between two sets of layer thicknesses.

   Each destination layer gets the thickness-weighted mean of the source layers it
   overlaps (in pressure). Destination layers without overlap (zero thickness, or below
   the source column) get the value of the source layer at their top interface.

   Args:
      val: Source values, shape (kdm, jdm, idm).
      dp_src: Source layer thicknesses (pressure units), same shape.
      dp_dst: Destination layer thicknesses (pressure units), same shape.

   Returns:
      Remapped values, shape (kdm, jdm, idm).
   """
   kdm = val.shape[0]
   p_src = numpy.concatenate([numpy.zeros((1,) + val.shape[1:]), numpy.cumsum(dp_src, axis=0)])
   p_dst = numpy.concatenate([numpy.zeros((1,) + val.shape[1:]), numpy.cumsum(dp_dst, axis=0)])
   out = numpy.zeros(val.shape)
   for kd in range(kdm):
      top, bot = p_dst[kd], p_dst[kd+1]
      acc = numpy.zeros(val.shape[1:])
      cov = numpy.zeros(val.shape[1:])
      for ks in range(kdm):
         ovl = numpy.clip(numpy.minimum(bot, p_src[ks+1]) - numpy.maximum(top, p_src[ks]), 0., None)
         acc += ovl*val[ks]
         cov += ovl
      # source layer containing the top interface of the destination layer
      ks_top = numpy.clip((p_src[1:] <= top[None]).sum(axis=0), 0, kdm-1)
      nearest = numpy.take_along_axis(val, ks_top[None], axis=0)[0]
      out[kd] = numpy.where(cov > 0., acc/numpy.where(cov > 0., cov, 1.), nearest)
   return out


def relax_layers(relax_int: str, month: int, kdm: int, idm: int, jdm: int,
                 pbot: numpy.ndarray) -> Tuple[numpy.ndarray, List[float]]:
   """ Layer thicknesses of the relax.intf layers for one month, as used by HYCOM.

   HYCOM uses the interfaces of relax.intf as the tops of layers 1..kdm and the model
   bottom as the bottom of layer kdm (initrc in trcupd.F90).

   Args:
      relax_int: relax_int file without .a/.b.
      month: Month (1-12).
      kdm, idm, jdm: Grid dimensions.
      pbot: Bottom pressure of the model (sum of layer thicknesses), shape (jdm, idm).

   Returns:
      Layer thicknesses, shape (kdm, jdm, idm), and the layer densities from the .b file.
   """
   rf = abfile.ABFileRelax(relax_int, "r", idm=idm, jdm=jdm)
   pw = numpy.empty((kdm + 1, jdm, idm))
   dens = []
   for k in range(1, kdm+1):
      fld = rf.read_field("int", k, month)
      if fld is None:
         raise ValueError("No field int, layer %d, month %d in %s" % (k, month, relax_int))
      pw[k-1] = numpy.ma.filled(fld, 0.)
      dens.append([d["dens"] for d in rf._fields.values() if d["k"] == k and d["month"] == month][0])
   rf.close()
   pw[kdm] = pbot
   pw = numpy.minimum(numpy.maximum.accumulate(pw, axis=0), pbot[None])
   return numpy.diff(pw, axis=0), dens


def main(archv: str, bgc_archv: str, relax_int: str, outdir: str) -> List[str]:
   """ Write relax.<tracer>.[ab] for every tracer in the BGC archive.

   Args:
      archv: Physics archive without .a/.b (layer thicknesses, ocean mask, time).
      bgc_archv: BGC archive without .a/.b, on the layers of archv, valid at the same time.
      relax_int: relax_int file of the experiment (relax/<E>/relax_int) without .a/.b.
      outdir: Output directory.

   Returns:
      Names of the tracers written.
   """
   arc = abfile.ABFileArchv(archv, "r")
   idm, jdm = arc.idm, arc.jdm
   kdm = max(arc.fieldlevels)
   dtime = float(next(iter(arc._fields.values()))["day"])
   ocean = ~numpy.ma.getmaskarray(arc.read_field("thknss", 1))
   dp = numpy.array([numpy.ma.filled(arc.read_field("thknss", k), 0.) for k in range(1, kdm+1)])
   dp[:, ~ocean] = 0.
   arc.close()
   pbot = dp.sum(axis=0)
   logger.info("Physics archive %s: model day %.5f, %d layers, %d ocean cells" % (archv, dtime, kdm, ocean.sum()))

   barc = abfile.ABFileArchv(bgc_archv, "r")
   bday = float(next(iter(barc._fields.values()))["day"])
   if abs(bday - dtime) > 1e-6:
      raise ValueError("BGC archive %s is valid at model day %.5f, physics archive at %.5f" % (bgc_archv, bday, dtime))
   names = sorted(barc.fieldnames)
   logger.info("Tracers in BGC archive: %s" % ", ".join(names))

   layers = [relax_layers(relax_int, m, kdm, idm, jdm, pbot) for m in range(1, 13)]

   for name in names:
      fname = os.path.join(outdir, "relax.%s" % name)
      for ext in (".a", ".b"):
         if os.path.exists(fname + ext):
            raise ValueError("%s exists, will not overwrite" % (fname + ext))
      val = numpy.empty(dp.shape)
      cells = numpy.zeros(ocean.shape, dtype=bool)
      for k in range(1, kdm+1):
         fld = numpy.ma.filled(barc.read_field(name, k), numpy.nan).astype(float)
         fld[fld > _missing] = numpy.nan
         val[k-1], bad = fill_nan(fld, ocean)
         cells |= bad
      if cells.any():
         j, i = numpy.where(cells)
         logger.warning("%s: NaN in %d ocean cells of BGC archive, filled by nearest neighbour; (j,i)=%s" %
                        (name, cells.sum(), list(zip(j[:5].tolist(), i[:5].tolist()))))

      afile = abfile.AFile(idm, jdm, fname + ".a", "w", mask=True)
      with open(fname + ".b", "w") as bfile:
         bfile.write("%-79s\n" % "BGC initial state for HYCOM (relax format, 12 identical months)")
         bfile.write("%-79s\n" % ("From %s" % os.path.basename(bgc_archv))[:79])
         bfile.write("%-79s\n" % "Layered on relax.intf interfaces of each month")
         bfile.write("%-79s\n" % name)
         bfile.write("i/jdm = %4d %4d\n" % (idm, jdm))
         for m in range(1, 13):
            dp_relax, dens = layers[m-1]
            rval = remap_layers(val, dp, dp_relax)
            for k in range(1, kdm+1):
               fld = numpy.where(ocean, rval[k-1], 0.)
               hmin, hmax = afile.writerecord(fld, ~ocean)
               bfile.write("%s: month,layer,dens,range = %02d  %02d %6.3f %15.7E %15.7E\n" %
                           ("trc", m, k, dens[k-1], hmin, hmax))
      afile.close()
      logger.info("Wrote %s.[ab]" % fname)
   barc.close()
   return names


if __name__ == "__main__":
   parser = argparse.ArgumentParser(description="Create BGC initialisation files relax.<tracer>.[ab] from a BGC nesting archive")
   parser.add_argument("archv", help="Physics archive without .a/.b, e.g. ../nest/022/archv.1993_244_00")
   parser.add_argument("bgc_archv", help="BGC archive without .a/.b, e.g. ../nest/022/archv_fabm.1993_244_00")
   parser.add_argument("--relax-int", default="../relax/${E}/relax_int",
                       help="relax_int of the experiment without .a/.b (default ../relax/$E/relax_int)")
   parser.add_argument("--outdir", default=".", help="Output directory (default .)")
   args = parser.parse_args()
   main(args.archv, args.bgc_archv, os.path.expandvars(args.relax_int), args.outdir)

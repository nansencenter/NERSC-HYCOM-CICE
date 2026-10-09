#!/usr/bin/env python
"""
Wrap HYCOM and CICE restart files from the end of one spin-up cycle to the start
of the next.

The restart pair valid at from_date (end of cycle n) is copied to the data
directory of cycle n+1 as a restart pair valid at to_date (start of cycle n+1).
The model state is unchanged; only the time information is shifted back by
(from_date - to_date), so the files look as if the model had written them at to_date:

   HYCOM restart.YYYY_DDD_HH_0000.[ab] : .a copied unchanged; nstep and dtime in the
                                          .b header shifted (HYCOM checks nstep against
                                          the start time in limits)
   CICE  cice/iced.YYYY-MM-DD-SSSSS.nc : global attributes istep1 and time shifted,
                                          nyr, month, mday, sec set to to_date

Both dates must be at the same time of day, and both old and new times are checked
against the HYCOM (blkdat.input) and CICE (ice_in) calendars.
"""
import argparse
import datetime
import logging
import os
import re
import shutil
from typing import Tuple

import f90nml
import modeltools.hycom
import netCDF4

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


def hycom_dtime(date: datetime.datetime, yrflag: int) -> float:
   """ HYCOM model day of a date, computed as in hycom_limits.py.

   Args:
      date: Date and time.
      yrflag: HYCOM calendar flag from blkdat.input.

   Returns:
      Model day (dtime).
   """
   fday, fhour, fsecond = modeltools.hycom.datetime_to_ordinal(date, yrflag)
   return modeltools.hycom.dayfor(date.year, fday, fhour, yrflag) + fsecond/86400.


def hycom_restart_name(date: datetime.datetime, prefix: str = "restart") -> str:
   """ HYCOM restart file name (without .a/.b) valid at date, as in expt_preprocess.sh.

   Args:
      date: Date and time.
      prefix: Restart file prefix (nmrsti in blkdat.input, default "restart").

   Returns:
      e.g. restart.1998_001_00_0000
   """
   hsec = date.minute*60 + date.second
   return "%s.%04d_%03d_%02d_%04d" % (prefix, date.year, date.timetuple().tm_yday, date.hour, hsec)


def cice_restart_name(date: datetime.datetime, restart_file: str = "iced") -> str:
   """ CICE restart file name valid at date, as in expt_preprocess.sh.

   Args:
      date: Date and time.
      restart_file: restart_file in ice_in (default "iced").

   Returns:
      e.g. iced.1998-01-01-00000.nc
   """
   dsec = date.hour*3600 + date.minute*60 + date.second
   return "%s.%s-%05d.nc" % (restart_file, date.strftime("%Y-%m-%d"), dsec)


def wrap_hycom(src: str, dst: str, from_date: datetime.datetime, to_date: datetime.datetime,
               yrflag: int) -> Tuple[int, float]:
   """ Copy a HYCOM restart pair and shift nstep/dtime in the .b header.

   Args:
      src: Source restart name without .a/.b, valid at from_date.
      dst: Destination restart name without .a/.b, valid at to_date.
      from_date: Valid time of src.
      to_date: Valid time of dst.
      yrflag: HYCOM calendar flag.

   Returns:
      New nstep and dtime.
   """
   with open(src + ".b") as f:
      lines = f.readlines()
   m = re.match(r"(RESTART2?: nstep,dtime(?:,thbase)? =)\s*(\S+)\s+(\S+)(.*)$", lines[1].rstrip("\n"))
   if not m:
      raise ValueError("Unexpected second header line in %s.b: %s" % (src, lines[1]))
   nstep, dtime = int(m.group(2)), float(m.group(3))
   rest = m.group(4).split()

   dtime_from, dtime_to = hycom_dtime(from_date, yrflag), hycom_dtime(to_date, yrflag)
   if abs(dtime - dtime_from) > 1e-6:
      raise ValueError("%s.b has dtime %.5f, expected %.5f for %s" % (src, dtime, dtime_from, from_date))
   steps_per_day = nstep/dtime
   if abs(steps_per_day - round(steps_per_day)) > 1e-6:
      raise ValueError("nstep/dtime = %.6f in %s.b is not an integer number of steps per day" % (steps_per_day, src))
   shift = dtime_from - dtime_to
   nstep_new = nstep - int(round(shift*steps_per_day))
   dtime_new = dtime_to

   line = "%s %12d%19.10f" % (m.group(1), nstep_new, dtime_new)
   if rest:
      line += "%24.13f" % float(rest[0])
   lines[1] = line + "\n"
   shutil.copyfile(src + ".a", dst + ".a")
   with open(dst + ".b", "w") as f:
      f.writelines(lines)
   logger.info("HYCOM: %s -> %s, nstep %d -> %d, dtime %.5f -> %.5f" %
               (os.path.basename(src), dst, nstep, nstep_new, dtime, dtime_new))
   return nstep_new, dtime_new


def wrap_cice(src: str, dst: str, from_date: datetime.datetime, to_date: datetime.datetime,
              dt: float, year_init: int) -> Tuple[int, float]:
   """ Copy a CICE netCDF restart file and shift its time attributes.

   Args:
      src: Source restart file, valid at from_date.
      dst: Destination restart file, valid at to_date.
      from_date: Valid time of src.
      to_date: Valid time of dst.
      dt: CICE time step (s), from ice_in.
      year_init: year_init from ice_in; time is counted in seconds from year_init-01-01.

   Returns:
      New istep1 and time.
   """
   t0 = datetime.datetime(year_init, 1, 1)
   time_from = (from_date - t0).total_seconds()
   time_to = (to_date - t0).total_seconds()
   shift = time_from - time_to
   if abs(shift/dt - round(shift/dt)) > 1e-9:
      raise ValueError("Wrap interval %g s is not a multiple of the CICE time step %g s" % (shift, dt))

   shutil.copyfile(src, dst)
   nc = netCDF4.Dataset(dst, "a")
   istep1, time = int(nc.getncattr("istep1")), float(nc.getncattr("time"))
   if abs(time - time_from) > 1e-3:
      nc.close()
      os.remove(dst)
      raise ValueError("%s has time %.1f s, expected %.1f s for %s (year_init=%d)" % (src, time, time_from, from_date, year_init))
   istep1_new = istep1 - int(round(shift/dt))
   nc.setncattr("istep1", netCDF4.numpy.int32(istep1_new))
   nc.setncattr("time", netCDF4.numpy.float64(time_to))
   if "nyr" in nc.ncattrs():
      nc.setncattr("nyr", netCDF4.numpy.int32(to_date.year - year_init + 1))
      nc.setncattr("month", netCDF4.numpy.int32(to_date.month))
      nc.setncattr("mday", netCDF4.numpy.int32(to_date.day))
      nc.setncattr("sec", netCDF4.numpy.int32(to_date.hour*3600 + to_date.minute*60 + to_date.second))
   nc.close()
   logger.info("CICE : %s -> %s, istep1 %d -> %d, time %.0f -> %.0f" %
               (os.path.basename(src), dst, istep1, istep1_new, time, time_to))
   return istep1_new, time_to


def main(from_dir: str, to_dir: str, from_date: datetime.datetime, to_date: datetime.datetime,
         blkdat: str, ice_in: str) -> None:
   """ Wrap the HYCOM and CICE restarts at from_date in from_dir to to_date in to_dir.

   Args:
      from_dir: Data directory of the finished cycle.
      to_dir: Data directory of the next cycle (created if necessary).
      from_date: End date of the finished cycle.
      to_date: Start date of the next cycle.
      blkdat: blkdat.input of the experiment (yrflag).
      ice_in: CICE namelist of the experiment (dt, year_init, restart_dir, restart_file).
   """
   if (from_date - to_date).total_seconds() % 86400 != 0:
      raise ValueError("from_date and to_date must be at the same time of day")
   yrflag = int(modeltools.hycom.BlkdatParser(blkdat)["yrflag"])
   nml = f90nml.read(ice_in)["setup_nml"]
   restart_dir = nml["restart_dir"].strip("./").strip("/")

   src_h = os.path.join(from_dir, hycom_restart_name(from_date))
   dst_h = os.path.join(to_dir, hycom_restart_name(to_date))
   src_c = os.path.join(from_dir, restart_dir, cice_restart_name(from_date, nml["restart_file"]))
   dst_c = os.path.join(to_dir, restart_dir, cice_restart_name(to_date, nml["restart_file"]))
   for f in (src_h + ".a", src_h + ".b", src_c):
      if not os.path.exists(f):
         raise FileNotFoundError(f)
   for f in (dst_h + ".a", dst_h + ".b", dst_c):
      if os.path.exists(f):
         raise FileExistsError("%s exists, will not overwrite" % f)
   os.makedirs(os.path.dirname(dst_c), exist_ok=True)

   wrap_hycom(src_h, dst_h, from_date, to_date, yrflag)
   wrap_cice(src_c, dst_c, from_date, to_date, float(nml["dt"]), int(nml["year_init"]))


if __name__ == "__main__":
   class DateTimeParseAction(argparse.Action):
      """ Parse YYYY-mm-ddTHH:MM:SS into a datetime """
      def __call__(self, parser: argparse.ArgumentParser, args: argparse.Namespace, values: str,
                   option_string: str = None) -> None:
         setattr(args, self.dest, datetime.datetime.strptime(values, "%Y-%m-%dT%H:%M:%S"))

   parser = argparse.ArgumentParser(description="Wrap HYCOM/CICE restarts from the end of one spin-up cycle to the start of the next")
   parser.add_argument("from_dir", help="Data dir of the finished cycle, e.g. data/cycle_01")
   parser.add_argument("to_dir", help="Data dir of the next cycle, e.g. data/cycle_02")
   parser.add_argument("from_date", action=DateTimeParseAction, help="End of finished cycle, YYYY-mm-ddTHH:MM:SS")
   parser.add_argument("to_date", action=DateTimeParseAction, help="Start of next cycle, YYYY-mm-ddTHH:MM:SS")
   parser.add_argument("--blkdat", default="blkdat.input", help="blkdat.input (default ./blkdat.input)")
   parser.add_argument("--ice-in", default="ice_in", help="CICE namelist (default ./ice_in)")
   args = parser.parse_args()
   main(args.from_dir, args.to_dir, args.from_date, args.to_date, args.blkdat, args.ice_in)

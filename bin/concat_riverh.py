"""Join GloFAS river forcing chunks (riverh_*.[ab]) into one riverh.[ab], with checks.

The GloFAS job array in docs/forcing.md runs hycom_river_hifre.py for consecutive
periods ("chunks", one per year) in parallel. Each chunk except the last runs one day
into the next chunk, because a day's records are only written once the following day is
processed. Where chunks overlap, the records of the later chunk are used, so the result
is the same as one continuous river_synoptic.sh run.

Records are copied byte for byte. With --from/--to only records whose model day lies
in [from, to] (dates, both inclusive) are kept, e.g. to cut out a single year.

The script stops with an error, and removes its output, if
  - the chunks are not on the same grid, or a chunk is shorter than its .b file says,
  - the records are not evenly spaced in time (a missing or incomplete chunk, or a chunk
    from another period),
  - the output does not start/end on the days given with --from/--to,
  - a record contains NaN, negative values or values above MAXFLUX (land is ignored).

Usage:
   python concat_riverh.py [--from YYYY-MM-DD] [--to YYYY-MM-DD] out chunk1.a chunk2.a ...
   (out without .a/.b; chunks in time order)
"""
import argparse
import datetime
import os
import re
import sys

import numpy as np

DTIME_RE = re.compile(r"dtime1,range\s*=\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)")
MODEL_DAY0 = datetime.date(1900, 12, 31)
LAND = 1e29        # HYCOM land/padding values are 2**100
MAXFLUX = 1.e-3    # [m/s], far above any real river forcing


def read_b(bfile: str) -> tuple[list[str], list[str], int, int]:
    """Return (header lines, record lines, idm, jdm) of a river forcing .b file."""
    lines = open(bfile).read().splitlines()
    nhead = next(i for i, line in enumerate(lines) if DTIME_RE.search(line))
    m = re.search(r"i/jdm\s*=\s*(\d+)\s+(\d+)", "\n".join(lines[:nhead]))
    if not m:
        sys.exit("Error: could not find i/jdm in %s" % bfile)
    return lines[:nhead], lines[nhead:], int(m.group(1)), int(m.group(2))


def dtime(line: str) -> tuple[float, float]:
    """Return (model day, record interval in days) of a .b record line."""
    m = DTIME_RE.search(line)
    return float(m.group(1)), float(m.group(2))


def model_day(date: str) -> float:
    """Convert YYYY-MM-DD to a model day (days since 1900-12-31)."""
    return float((datetime.date.fromisoformat(date) - MODEL_DAY0).days)


def date_of(t: float) -> str:
    """Convert a model day to YYYY-MM-DD HH:MM."""
    return (datetime.datetime(1900, 12, 31) + datetime.timedelta(days=t)).strftime("%Y-%m-%d %H:%M")


def check_values(data: bytes, n: int, where: str) -> None:
    """Stop if the first n values of a record have NaN, negative or too large ocean values."""
    v = np.frombuffer(data, dtype=">f4", count=n)
    if np.isnan(v).any():
        sys.exit("Error: NaN in %s" % where)
    ocean = v[v < LAND]
    if ocean.size and (ocean.min() < 0. or ocean.max() > MAXFLUX):
        sys.exit("Error: river flux in %s outside 0 to %g m/s (min %g, max %g)"
                 % (where, MAXFLUX, ocean.min(), ocean.max()))


def concat(out: str, chunks: list[str], first_day: str | None, last_day: str | None) -> None:
    """Join the chunks into out.[ab]; see the module docstring for the checks."""
    tmin = model_day(first_day) if first_day else -1e30
    tmax = model_day(last_day) + 1 if last_day else 1e30

    bs = [read_b(c[:-2] + ".b") for c in chunks]
    header, _, idm, jdm = bs[0]
    if any((b[2], b[3]) != (idm, jdm) for b in bs):
        sys.exit("Error: chunks are not all on the same %dx%d grid" % (idm, jdm))
    n = idm * jdm
    nbytes = ((n + 4095) // 4096) * 4096 * 4

    nout = 0
    first = last = dt = None
    with open(out + ".a", "wb") as fa, open(out + ".b", "w") as fb:
        fb.write("\n".join(header) + "\n")
        for k, (chunk, (_, recs, _, _)) in enumerate(zip(chunks, bs)):
            # records at or after the start of the next chunk come from the next chunk
            tnext = dtime(bs[k + 1][1][0])[0] if k + 1 < len(bs) else 1e30
            with open(chunk, "rb") as fc:
                for i, line in enumerate(recs):
                    t, dt_rec = dtime(line)
                    if t >= tnext or t < tmin or t >= tmax:
                        continue
                    dt = dt_rec if dt is None else dt
                    if abs(dt_rec - dt) > 1e-6:
                        sys.exit("Error: record interval changes from %g to %g days at %s, record %d"
                                 % (dt, dt_rec, chunk, i + 1))
                    if last is not None and abs(t - last - dt) > 1e-6:
                        sys.exit("Error: gap or overlap between %s and %s (%s, record %d): records are "
                                 "missing, a chunk is incomplete or from another period"
                                 % (date_of(last), date_of(t), chunk, i + 1))
                    fc.seek(i * nbytes)
                    data = fc.read(nbytes)
                    if len(data) != nbytes:
                        sys.exit("Error: %s is shorter than its .b file says" % chunk)
                    check_values(data, n, "%s, record %d (%s)" % (chunk, i + 1, date_of(t)))
                    fa.write(data)
                    fb.write(line + "\n")
                    first = t if first is None else first
                    last = t
                    nout += 1
    if nout == 0:
        sys.exit("Error: no records selected")
    if first_day and abs(first - tmin) > 1e-6:
        sys.exit("Error: output starts on %s, not on %s" % (date_of(first), first_day))
    if last_day and last < tmax - 1 - 1e-6:
        sys.exit("Error: output ends on %s, before %s" % (date_of(last), last_day))
    print("Wrote %s.[ab]: %d records every %g days, %s to %s, all checks passed"
          % (out, nout, dt, date_of(first), date_of(last)))


def main() -> None:
    """Parse arguments, join the chunks and remove the output if a check fails."""
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--from", dest="first", help="first day to keep, YYYY-MM-DD")
    p.add_argument("--to", dest="last", help="last day to keep (inclusive), YYYY-MM-DD")
    p.add_argument("out", help="output name, without .a/.b")
    p.add_argument("chunks", nargs="+", help="chunk .a files, in time order")
    a = p.parse_args()
    try:
        concat(a.out, a.chunks, a.first, a.last)
    except SystemExit:
        for f in (a.out + ".a", a.out + ".b"):
            if os.path.exists(f):
                os.remove(f)
        raise


if __name__ == "__main__":
    main()

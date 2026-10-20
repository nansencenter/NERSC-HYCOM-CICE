"""Join GloFAS river forcing chunks (riverh_*.[ab]) into one riverh.[ab].

river_synoptic_slurm.sh runs hycom_river_hifre.py for consecutive periods
("chunks") in parallel. Each chunk except the last runs one day into the next chunk,
because a day's records are only written once the following day is processed.
Where chunks overlap, the records of the later chunk are used, so the result is the
same as one continuous river_synoptic.sh run.

Records are copied byte for byte. With --from/--to only records whose model day lies
in [from, to] (dates, both inclusive) are kept, e.g. to cut out a single year.

Usage:
   python concat_riverh.py [--from YYYY-MM-DD] [--to YYYY-MM-DD] out chunk1.a chunk2.a ...
   (out without .a/.b; chunks in time order)
"""
import argparse
import datetime
import re
import sys

DTIME_RE = re.compile(r"dtime1,range\s*=\s*([-+0-9.eE]+)")
MODEL_DAY0 = datetime.date(1900, 12, 31)


def read_b(bfile: str) -> tuple[list[str], list[str], int, int]:
    """Return (header lines, record lines, idm, jdm) of a river forcing .b file."""
    lines = open(bfile).read().splitlines()
    nhead = next(i for i, line in enumerate(lines) if DTIME_RE.search(line))
    m = re.search(r"i/jdm\s*=\s*(\d+)\s+(\d+)", "\n".join(lines[:nhead]))
    if not m:
        sys.exit("Error: could not find i/jdm in %s" % bfile)
    return lines[:nhead], lines[nhead:], int(m.group(1)), int(m.group(2))


def dtime(line: str) -> float:
    """Return the model day (days since 1900-12-31) of a .b record line."""
    return float(DTIME_RE.search(line).group(1))


def model_day(date: str) -> float:
    """Convert YYYY-MM-DD to a model day (days since 1900-12-31)."""
    return float((datetime.date.fromisoformat(date) - MODEL_DAY0).days)


def main() -> None:
    """Parse arguments and write <out>.[ab]."""
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--from", dest="first", help="first day to keep, YYYY-MM-DD")
    p.add_argument("--to", dest="last", help="last day to keep (inclusive), YYYY-MM-DD")
    p.add_argument("out", help="output name, without .a/.b")
    p.add_argument("chunks", nargs="+", help="chunk .a files, in time order")
    a = p.parse_args()

    tmin = model_day(a.first) if a.first else -1e30
    tmax = model_day(a.last) + 1 if a.last else 1e30

    bs = [read_b(c[:-2] + ".b") for c in a.chunks]
    header, _, idm, jdm = bs[0]
    if any((b[2], b[3]) != (idm, jdm) for b in bs):
        sys.exit("Error: chunks are not all on the same %dx%d grid" % (idm, jdm))
    nbytes = ((idm * jdm + 4095) // 4096) * 4096 * 4

    nout = 0
    first = last = None
    with open(a.out + ".a", "wb") as fa, open(a.out + ".b", "w") as fb:
        fb.write("\n".join(header) + "\n")
        for k, (chunk, (_, recs, _, _)) in enumerate(zip(a.chunks, bs)):
            # records at or after the start of the next chunk come from the next chunk
            tnext = dtime(bs[k + 1][1][0]) if k + 1 < len(bs) else 1e30
            with open(chunk, "rb") as fc:
                for i, line in enumerate(recs):
                    t = dtime(line)
                    if t >= tnext or t < tmin or t >= tmax:
                        continue
                    if last is not None and t <= last:
                        sys.exit("Error: records not increasing in time at %s, record %d" % (chunk, i + 1))
                    fc.seek(i * nbytes)
                    data = fc.read(nbytes)
                    if len(data) != nbytes:
                        sys.exit("Error: %s is shorter than its .b file says" % chunk)
                    fa.write(data)
                    fb.write(line + "\n")
                    first = t if first is None else first
                    last = t
                    nout += 1
    if nout == 0:
        sys.exit("Error: no records selected")
    print("Wrote %s.[ab]: %d records, model days %.4f to %.4f" % (a.out, nout, first, last))


if __name__ == "__main__":
    main()

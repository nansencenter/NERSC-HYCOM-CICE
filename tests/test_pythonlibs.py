"""
Phase 1 tests — Python library readers against real TP0 data.

All tests run against topo files committed in the repo
(TP0a1.00/topo/regional.grid.[ab] and depth_TP0a1.00_01.[ab]).
No generated or synthetic data; no skips for missing files.
"""

import datetime
import os
import sys

import numpy as np
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "pythonlibs", "abfile"))
sys.path.insert(0, os.path.join(REPO_ROOT, "pythonlibs", "modeltools"))

import importlib.util

import abfile.abfile as abf


def _load_submodule(pkg_dir, rel_path):
    """Load a modeltools submodule directly, bypassing the package __init__
    (which pulls in cfunits / other optional deps we don't need here)."""
    path = os.path.join(REPO_ROOT, "pythonlibs", "modeltools", pkg_dir, rel_path)
    name = rel_path.replace(os.sep, ".").removesuffix(".py")
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_blkdat    = _load_submodule("modeltools/hycom", "_blkdat.py")
_timetools = _load_submodule("modeltools/hycom", "_timetools.py")

BlkdatParser        = _blkdat.BlkdatParser
datetime_to_ordinal = _timetools.datetime_to_ordinal
dayfor              = _timetools.dayfor
forday_datetime     = _timetools.forday_datetime


# ── ABFileGrid ────────────────────────────────────────────────────────────────

class TestABFileGrid:
    def test_dimensions(self, tp0_topo):
        grid = abf.ABFileGrid(tp0_topo["grid"], "r")
        assert grid.idm == tp0_topo["idm"]
        assert grid.jdm == tp0_topo["jdm"]

    def test_plon_shape_and_range(self, tp0_topo):
        grid = abf.ABFileGrid(tp0_topo["grid"], "r")
        plon = grid.read_field("plon")
        assert plon.shape == (tp0_topo["jdm"], tp0_topo["idm"])
        assert float(plon.min()) >= -180.0
        assert float(plon.max()) <= 360.0

    def test_plat_shape_and_range(self, tp0_topo):
        grid = abf.ABFileGrid(tp0_topo["grid"], "r")
        plat = grid.read_field("plat")
        assert plat.shape == (tp0_topo["jdm"], tp0_topo["idm"])
        assert float(plat.min()) >= -90.0
        assert float(plat.max()) <=  90.0

    def test_plon_minmax_from_bfile_header(self, tp0_topo):
        """Cross-check array min/max against the .b header by reading the file directly."""
        import re
        grid = abf.ABFileGrid(tp0_topo["grid"], "r")
        plon = grid.read_field("plon")
        with open(tp0_topo["grid"] + ".b") as f:
            for line in f:
                m = re.match(r"plon:\s+min,max\s*=\s*(\S+)\s+(\S+)", line)
                if m:
                    bmin, bmax = float(m.group(1)), float(m.group(2))
                    np.testing.assert_allclose(float(plon.min()), bmin, atol=1e-4)
                    np.testing.assert_allclose(float(plon.max()), bmax, atol=1e-4)
                    return
        raise AssertionError("plon not found in .b header")

    def test_plat_minmax_from_bfile_header(self, tp0_topo):
        import re
        grid = abf.ABFileGrid(tp0_topo["grid"], "r")
        plat = grid.read_field("plat")
        with open(tp0_topo["grid"] + ".b") as f:
            for line in f:
                m = re.match(r"plat:\s+min,max\s*=\s*(\S+)\s+(\S+)", line)
                if m:
                    bmin, bmax = float(m.group(1)), float(m.group(2))
                    np.testing.assert_allclose(float(plat.min()), bmin, atol=1e-4)
                    np.testing.assert_allclose(float(plat.max()), bmax, atol=1e-4)
                    return
        raise AssertionError("plat not found in .b header")


# ── ABFileBathy ───────────────────────────────────────────────────────────────

class TestABFileBathy:
    def test_reads_depth_field(self, tp0_topo):
        bathy = abf.ABFileBathy(
            tp0_topo["bathy"], "r",
            idm=tp0_topo["idm"], jdm=tp0_topo["jdm"],
        )
        depth = bathy.read_field("depth")
        assert depth is not None
        assert depth.shape == (tp0_topo["jdm"], tp0_topo["idm"])

    def test_ocean_depths_positive(self, tp0_topo):
        bathy = abf.ABFileBathy(
            tp0_topo["bathy"], "r",
            idm=tp0_topo["idm"], jdm=tp0_topo["jdm"],
        )
        depth = bathy.read_field("depth")
        ocean = depth.compressed() if hasattr(depth, "compressed") else depth.ravel()
        assert float(ocean.min()) > 0.0

    def test_has_land_points(self, tp0_topo):
        bathy = abf.ABFileBathy(
            tp0_topo["bathy"], "r",
            idm=tp0_topo["idm"], jdm=tp0_topo["jdm"],
        )
        depth = bathy.read_field("depth")
        assert hasattr(depth, "mask") and depth.mask.any(), \
            "Expected masked land points in TP0 bathy"

    def test_depth_matches_bfile_minmax(self, tp0_topo):
        bathy = abf.ABFileBathy(
            tp0_topo["bathy"], "r",
            idm=tp0_topo["idm"], jdm=tp0_topo["jdm"],
        )
        depth = bathy.read_field("depth")
        bmin, bmax = bathy.bminmax("depth")
        ocean = depth.compressed() if hasattr(depth, "compressed") else depth.ravel()
        np.testing.assert_allclose(float(ocean.min()), bmin, atol=1e-2)
        np.testing.assert_allclose(float(ocean.max()), bmax, atol=1e-2)


# ── BlkdatParser ─────────────────────────────────────────────────────────────

class TestBlkdatParser:
    def test_ice_flag(self, tp0_blkdat):
        bp = BlkdatParser(tp0_blkdat)
        # iceflg is not in BlkdatParser._integer_fields so it comes back as a string
        assert int(bp["iceflg"]) == 2, "TP0 expt_01.0 must use ESMF/CICE coupling (iceflg=2)"

    def test_grid_dimensions(self, tp0_blkdat):
        bp = BlkdatParser(tp0_blkdat)
        assert bp["idm"] == 100
        assert bp["jdm"] == 110

    def test_no_bgc(self, tp0_blkdat):
        bp = BlkdatParser(tp0_blkdat)
        assert bp["ntracr"] == 0, "TP0 expt_01.0 should have no BGC tracers"

    def test_baroclinic_timestep(self, tp0_blkdat):
        bp = BlkdatParser(tp0_blkdat)
        assert bp["baclin"] == pytest.approx(1200.0)

    def test_sigma_version(self, tp0_blkdat):
        # thflag=0 means sigma-0 coordinate
        bp = BlkdatParser(tp0_blkdat)
        assert bp["thflag"] == 0


# ── Time tools ────────────────────────────────────────────────────────────────

class TestTimeTools:
    """Pure function tests — no file I/O, no data fixtures."""

    def test_datetime_to_ordinal_known_date(self):
        # 2006-07-09 12:00 UTC is day 190, hour 12 of year 2006 (yrflag=3)
        dt = datetime.datetime(2006, 7, 9, 12)
        iday, ihour, isecond = datetime_to_ordinal(dt, yrflag=3)
        assert iday == 190
        assert ihour == 12
        assert isecond == 0

    def test_dayfor_returns_float(self):
        dtime = dayfor(2006, 190, 12, 3)
        assert isinstance(dtime, float)
        assert dtime > 0.0

    def test_dayfor_datetime_to_ordinal_roundtrip(self):
        dt = datetime.datetime(2010, 3, 15, 6)
        iday, ihour, isecond = datetime_to_ordinal(dt, yrflag=3)
        dtime = dayfor(dt.year, iday, ihour, 3)
        dt2 = forday_datetime(dtime, 3)
        assert dt2.year  == dt.year
        assert dt2.month == dt.month
        assert dt2.day   == dt.day

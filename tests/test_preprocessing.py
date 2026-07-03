"""
Phase 1 tests — preprocessing script smoke tests.

Covers:
  - nemo2archvz_regular.py  (GLORYS → HYCOM archive in z-coordinates)
  - hycom_atmfor.py         (ERA5 → HYCOM atmospheric forcing files)

Both use synthetic NetCDF fixtures (conftest.glorys_fixtures / era5_fixtures)
that match the real GLORYS12/Mercator and ERA5 6-hourly file formats exactly.
All tests are skipped when the required heavy dependencies are not installed
(scipy, pyproj, metpy, netCDF4, cfunits).
"""

import importlib.util
import os
import subprocess
import sys
import types
from typing import Any

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0, os.path.join(REPO_ROOT, "pythonlibs", "abfile"))
sys.path.insert(0, os.path.join(REPO_ROOT, "pythonlibs", "modeltools"))
BIN_DIR   = os.path.join(REPO_ROOT, "bin")
TOPO_DIR  = os.path.join(REPO_ROOT, "TP0a1.00", "topo")
EXPT_DIR  = os.path.join(REPO_ROOT, "TP0a1.00", "expt_01.0")

_SCRIPT_ENV = {
    **os.environ,
    "PYTHONPATH": os.pathsep.join(filter(None, [
        os.path.join(REPO_ROOT, "pythonlibs", "abfile"),
        os.path.join(REPO_ROOT, "pythonlibs", "modeltools"),
        os.environ.get("PYTHONPATH", ""),
    ])),
}

_PREPROC_DEPS = ["netCDF4", "scipy", "pyproj", "cfunits", "metpy"]


def _require_preproc_deps() -> None:
    """Skip the calling test if any heavy preprocessing dependency is absent."""
    missing = [d for d in _PREPROC_DEPS if not _can_import(d)]
    if missing:
        pytest.skip(f"preprocessing deps not installed: {', '.join(missing)}")


def _can_import(name: str) -> bool:
    """Return True if *name* is importable in the current Python environment."""
    try:
        __import__(name)
        return True
    except ImportError:
        return False


def _load_script(name: str) -> types.ModuleType:
    """Import a bin/ script as a module without invoking its __main__ block."""
    path = os.path.join(BIN_DIR, name)
    spec = importlib.util.spec_from_file_location(name.removesuffix(".py"), path)
    mod  = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ── nemo2archvz_regular.py ────────────────────────────────────────────────────

class TestNemo2ArchvzRegular:
    def test_import(self):
        """Script imports without error — catches missing dependencies."""
        _require_preproc_deps()
        _load_script("nemo2archvz_regular.py")

    def test_end_to_end(self, glorys_fixtures, tmp_path):
        """
        nemo2archvz_regular.py runs on synthetic GLORYS data and writes
        an archv.YYYY_DDD_HH.[ab] file pair to the working directory.

        Invocation:
            nemo2archvz_regular.py <mask_bathy> <physics> --coord_file <coord>
        """
        _require_preproc_deps()
        (tmp_path / "data").mkdir()  # script writes output to ./data/archv.*

        result = subprocess.run(
            [
                sys.executable,
                os.path.join(BIN_DIR, "nemo2archvz_regular.py"),
                glorys_fixtures["mask"],
                glorys_fixtures["physics"],
                "--coord_file", glorys_fixtures["coord"],
            ],
            capture_output=True, text=True,
            env=_SCRIPT_ENV, cwd=str(tmp_path),
        )
        assert result.returncode == 0, (
            f"nemo2archvz_regular.py exited with code {result.returncode}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
        ab_files = list((tmp_path / "data").glob("archv.*.a"))
        assert ab_files, (
            f"No archv.*.a produced in {tmp_path}/data\nstderr:\n{result.stderr}"
        )


# ── hycom_atmfor.py ───────────────────────────────────────────────────────────

class TestHycomAtmfor:
    def test_import(self):
        """Script imports without error — catches missing dependencies."""
        _require_preproc_deps()
        _load_script("hycom_atmfor.py")

    def test_end_to_end(self, era5_fixtures, tmp_path):
        """
        hycom_atmfor.py runs on synthetic ERA5 data and writes
        forcing.*.[ab] files to the working directory.

        The script reads regional.grid.[ab] and blkdat.input from the CWD,
        so we symlink the committed TP0 topo and experiment files.

        Invocation:
            hycom_atmfor.py 1993-01-01T06:00:00 1993-01-01T06:00:00 <xml> ERA5_test
        """
        _require_preproc_deps()

        # symlink grid and blkdat into tmp_path so the script finds them
        for name in ("regional.grid.a", "regional.grid.b"):
            (tmp_path / name).symlink_to(os.path.join(TOPO_DIR, name))
        (tmp_path / "blkdat.input").symlink_to(
            os.path.join(EXPT_DIR, "blkdat.input")
        )

        result = subprocess.run(
            [
                sys.executable,
                os.path.join(BIN_DIR, "hycom_atmfor.py"),
                "1993-01-01T06:00:00",
                "1993-01-01T06:00:00",
                era5_fixtures["xml"],
                "ERA5_test",
            ],
            capture_output=True, text=True,
            env=_SCRIPT_ENV, cwd=str(tmp_path),
        )
        assert result.returncode == 0, (
            f"hycom_atmfor.py exited with code {result.returncode}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
        forcing_files = list(tmp_path.glob("forcing.*.a"))
        assert forcing_files, (
            f"No forcing.*.a produced in {tmp_path}\nstderr:\n{result.stderr}"
        )

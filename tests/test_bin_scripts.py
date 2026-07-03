"""
Phase 1 tests — bin/ script import and CLI smoke tests.

Checks that the scripts in bin/ are importable (i.e. their dependencies are
installed) and that the ones with pure-function logic produce correct output
for known inputs against committed TP0 data.
"""

import importlib.util
import os
import subprocess
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN_DIR   = os.path.join(REPO_ROOT, "bin")

sys.path.insert(0, os.path.join(REPO_ROOT, "pythonlibs", "abfile"))
sys.path.insert(0, os.path.join(REPO_ROOT, "pythonlibs", "modeltools"))

# Environment for subprocess calls — scripts find pythonlibs via PYTHONPATH
_SCRIPT_ENV = {
    **os.environ,
    "PYTHONPATH": os.pathsep.join(filter(None, [
        os.path.join(REPO_ROOT, "pythonlibs", "abfile"),
        os.path.join(REPO_ROOT, "pythonlibs", "modeltools"),
        os.environ.get("PYTHONPATH", ""),
    ])),
}


def _import_script(name):
    """Import a bin/ script as a module without executing __main__."""
    path = os.path.join(BIN_DIR, name)
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ── import smoke tests ────────────────────────────────────────────────────────

@pytest.mark.parametrize("script", [
    "hycom_date.py",
    "hycom_limits.py",
    "hycom_timesteps.py",
    "namelist_extract.py",
    # hycom_kapref.py / calc_montg1.py import matplotlib — excluded from routine CI
])
def test_script_imports(script):
    """Each script must be importable — catches missing dependencies."""
    _import_script(script)


# ── hycom_date.py ─────────────────────────────────────────────────────────────

class TestHycomDate:
    def test_datetime_to_ordinal_cli(self):
        result = subprocess.run(
            [sys.executable,
             os.path.join(BIN_DIR, "hycom_date.py"),
             "datetime", "ordinal", "2006-07-09T12:00:00"],
            capture_output=True, text=True, env=_SCRIPT_ENV,
        )
        assert result.returncode == 0, result.stderr
        # Output is (year, day-of-year, hour) — expect day 190, hour 12
        assert "190" in result.stdout
        assert "12"  in result.stdout

    def test_ordinal_to_datetime_cli(self):
        result = subprocess.run(
            [sys.executable,
             os.path.join(BIN_DIR, "hycom_date.py"),
             "ordinal", "datetime", "2006", "190", "12"],
            capture_output=True, text=True, env=_SCRIPT_ENV,
        )
        assert result.returncode == 0, result.stderr
        assert "2006-07-09" in result.stdout


# ── namelist_extract.py ───────────────────────────────────────────────────────

class TestNamelistExtract:
    def test_extract_ice_in_value(self, tp0_ice_in):
        result = subprocess.run(
            [sys.executable,
             os.path.join(BIN_DIR, "namelist_extract.py"),
             tp0_ice_in, "setup_nml", "days_per_year"],
            capture_output=True, text=True, env=_SCRIPT_ENV,
        )
        assert result.returncode == 0, result.stderr
        # ice_in for TP0 uses 365 days/year
        assert "365" in result.stdout


# ── hycom_timesteps.py ────────────────────────────────────────────────────────

class TestHycomTimesteps:
    def test_valid_timesteps_for_tp0(self, tp0_blkdat, tmp_path):
        """Script reads blkdat.input and prints valid baroclinic timestep options."""
        result = subprocess.run(
            [sys.executable,
             os.path.join(BIN_DIR, "hycom_timesteps.py"),
             tp0_blkdat, "600", "1800"],  # search range 600–1800 s covers TP0's 1200 s
            capture_output=True, text=True, env=_SCRIPT_ENV, cwd=tmp_path,
        )
        assert result.returncode == 0, result.stderr
        # Script outputs via logging (stderr); 1200 s is valid for TP0
        assert "1200" in result.stderr

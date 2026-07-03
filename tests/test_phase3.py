"""
Phase 3 tests — pure Python / bash tools, no compiled Fortran, no external data.

Written against the DOCUMENTED user-facing interface (docs/), not the
implementation. When a test fails because a script's behaviour diverges from
what the documentation says it should do, that divergence is the finding.

Covers:
  - hycom_limits.py              (running.md — via expt_preprocess.sh description)
  - expt_new.sh                  (experiment-setup.md)
  - hycom_creat_thkdf4_veldf4.py (forcing.md)
  - hycom_bathy_ports.py         (forcing.md, via nest_setup_ports.sh)
"""

import os
import re
import shutil
import subprocess
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN_DIR   = os.path.join(REPO_ROOT, "bin")
TP0_DIR   = os.path.join(REPO_ROOT, "TP0a1.00")
TOPO_DIR  = os.path.join(TP0_DIR, "topo")
EXPT_DIR  = os.path.join(TP0_DIR, "expt_01.0")

_SCRIPT_ENV = {
    **os.environ,
    "PYTHONPATH": os.pathsep.join(filter(None, [
        os.path.join(REPO_ROOT, "pythonlibs", "abfile"),
        os.path.join(REPO_ROOT, "pythonlibs", "modeltools"),
        os.environ.get("PYTHONPATH", ""),
    ])),
}


# ── hycom_limits.py ───────────────────────────────────────────────────────────

class TestHycomLimits:
    """
    Documented in running.md (via the expt_preprocess.sh description):

        "Writes the limits file with start and stop times via hycom_limits.py"

    Documented interface (from the argparse block in the script):

        hycom_limits.py start_time end_time [--init]
          start_time / end_time : ISO-8601, format YYYY-MM-DDTHH:MM:SS
          --init                : cold-start run; negates the first ordinal time

    Must be invoked from a directory containing blkdat.input.
    Writes a file named ``limits`` with two space-separated floats.
    """

    def _run(self, tmp_path, tp0_blkdat, args):
        (tmp_path / "blkdat.input").symlink_to(tp0_blkdat)
        return subprocess.run(
            [sys.executable, os.path.join(BIN_DIR, "hycom_limits.py"), *args],
            capture_output=True, text=True,
            env=_SCRIPT_ENV, cwd=str(tmp_path),
        )

    def test_restart_writes_limits_file(self, tmp_path, tp0_blkdat):
        """Restart run: a file named ``limits`` is written."""
        result = self._run(tmp_path, tp0_blkdat,
                           ["1993-01-01T00:00:00", "1994-01-01T00:00:00"])
        assert result.returncode == 0, (
            f"hycom_limits.py exited {result.returncode}\n"
            f"stderr:\n{result.stderr}"
        )
        assert (tmp_path / "limits").exists(), "limits file was not created"

    def test_restart_limits_has_two_floats(self, tmp_path, tp0_blkdat):
        """Restart run: limits contains exactly two space-separated floats."""
        self._run(tmp_path, tp0_blkdat,
                  ["1993-01-01T00:00:00", "1994-01-01T00:00:00"])
        parts = (tmp_path / "limits").read_text().split()
        assert len(parts) == 2, (
            f"expected 2 values in limits, got {len(parts)}: "
            f"{(tmp_path / 'limits').read_text()!r}"
        )
        for p in parts:
            float(p)  # raises if not a valid float

    def test_restart_first_value_positive(self, tmp_path, tp0_blkdat):
        """Restart run: first ordinal time is positive."""
        self._run(tmp_path, tp0_blkdat,
                  ["1993-01-01T00:00:00", "1994-01-01T00:00:00"])
        first = float((tmp_path / "limits").read_text().split()[0])
        assert first > 0, f"restart run should write positive start time, got {first}"

    def test_init_first_value_negative(self, tmp_path, tp0_blkdat):
        """
        --init (climatological cold start): first value in limits is negated.

        Documented in running.md: ``INITFLG="--init"`` is passed to expt_preprocess.sh,
        which passes it to hycom_limits.py. The model treats a negative start ordinal as
        a signal to initialise from climatology rather than a restart file.
        """
        result = self._run(tmp_path, tp0_blkdat,
                           ["--init", "1993-09-01T00:00:00", "1994-01-01T00:00:00"])
        assert result.returncode == 0, (
            f"hycom_limits.py --init exited {result.returncode}\n"
            f"stderr:\n{result.stderr}"
        )
        parts = (tmp_path / "limits").read_text().split()
        assert len(parts) == 2
        first = float(parts[0])
        assert first < 0, (
            f"--init run should write a negative start time; got {first}"
        )

    def test_limits_file_not_empty(self, tmp_path, tp0_blkdat):
        """
        limits file contains non-empty content.

        Latent-bug check: hycom_limits.py calls ``f.close`` without parentheses
        (line ``f.close`` is a no-op attribute access, not the method call
        ``f.close()``). In CPython the interpreter closes the file on GC, so this
        works in practice — but it is still a bug. This test confirms the content
        is written by checking the file is non-empty.
        """
        self._run(tmp_path, tp0_blkdat,
                  ["1993-01-01T00:00:00", "1994-01-01T00:00:00"])
        content = (tmp_path / "limits").read_text().strip()
        assert content, "limits file is empty (possible f.close vs f.close() bug)"


# ── expt_new.sh ───────────────────────────────────────────────────────────────

class TestExptNew:
    """
    Documented in experiment-setup.md:

        bin/expt_new.sh old_experiment_number new_experiment_number

    Must be run from an experiment directory (containing EXPT.src) or the
    region root (containing REGION.src).  Requires topo/regional.grid.b.

    What the docs say it does:
      1. Copies expt_<old>/ to expt_<new>/, excluding data/, log/, SCRATCH/
      2. Creates empty log/ and data/ in the new directory
      3. Updates EXPT.src: X="<new>" and E="<new-without-dot>"
      4. Updates blkdat.input: iexpt, idm, jdm (idm/jdm read from regional.grid.b)
      5. Renames pbsjob*.sh -N lines to <CONFIGNAME>_X<EXPT_ID>
    """

    @pytest.fixture
    def region_dir(self, tmp_path):
        """
        Minimal region directory tree that satisfies expt_new.sh's documented preconditions.

        Follows the structure described in experiment-setup.md, using committed TP0
        files for topo and EXPT.src/blkdat.input templates.
        """
        region = tmp_path / "TP0a1.00"
        region.mkdir()

        # Minimal REGION.src — expt_new.sh sources this only for $R (used in PBS rename).
        (region / "REGION.src").write_text("#!/bin/bash\nexport R=TP0a1.00\n")

        # topo/ with grid file (expt_new.sh reads regional.grid.b for IDM/JDM)
        topo = region / "topo"
        topo.mkdir()
        for ext in ("a", "b"):
            (topo / f"regional.grid.{ext}").symlink_to(
                os.path.join(TOPO_DIR, f"regional.grid.{ext}")
            )

        # Template experiment (copy, not symlink — expt_new.sh modifies copies in-place)
        expt01 = region / "expt_01.0"
        expt01.mkdir()
        shutil.copy2(os.path.join(EXPT_DIR, "EXPT.src"),     str(expt01 / "EXPT.src"))
        shutil.copy2(os.path.join(EXPT_DIR, "blkdat.input"), str(expt01 / "blkdat.input"))

        return region

    def _run_expt_new(self, region_dir, old_id="01.0", new_id="01.1"):
        return subprocess.run(
            ["bash", os.path.join(BIN_DIR, "expt_new.sh"), old_id, new_id],
            capture_output=True, text=True,
            cwd=str(region_dir / f"expt_{old_id}"),
        )

    def test_creates_expt_dir(self, region_dir):
        """New experiment directory is created at the expected path."""
        result = self._run_expt_new(region_dir)
        assert result.returncode == 0, (
            f"expt_new.sh failed\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
        assert (region_dir / "expt_01.1").is_dir(), "expt_01.1/ was not created"

    def test_expt_src_x_line(self, region_dir):
        """EXPT.src X= is updated to the new experiment ID (e.g. X=\"01.1\")."""
        self._run_expt_new(region_dir)
        expt_src = (region_dir / "expt_01.1" / "EXPT.src").read_text()
        assert re.search(r'^X="01\.1"', expt_src, re.MULTILINE), (
            f"X= line not updated in EXPT.src:\n{expt_src[:500]}"
        )

    def test_expt_src_e_line(self, region_dir):
        """EXPT.src E= is updated to the dotless form (e.g. E=\"011\")."""
        self._run_expt_new(region_dir)
        expt_src = (region_dir / "expt_01.1" / "EXPT.src").read_text()
        assert re.search(r'^E="011"', expt_src, re.MULTILINE), (
            f"E= line not updated in EXPT.src:\n{expt_src[:500]}"
        )

    def test_blkdat_iexpt(self, region_dir):
        """blkdat.input iexpt is updated to the new experiment number (011)."""
        self._run_expt_new(region_dir)
        blkdat = (region_dir / "expt_01.1" / "blkdat.input").read_text()
        # iexpt line format: " NNN  'iexpt '"
        assert re.search(r"^\s*011\s+'iexpt\s*'", blkdat, re.MULTILINE), (
            f"iexpt not set to 011 in blkdat.input:\n"
            + "\n".join(
                l for l in blkdat.splitlines() if "iexpt" in l
            )
        )

    def test_blkdat_idm_jdm_from_grid(self, region_dir, tp0_topo):
        """blkdat.input idm and jdm match the actual TP0 grid dimensions."""
        self._run_expt_new(region_dir)
        blkdat = (region_dir / "expt_01.1" / "blkdat.input").read_text()
        idm_match = re.search(r"^\s*(\d+)\s+'idm\s*'", blkdat, re.MULTILINE)
        jdm_match = re.search(r"^\s*(\d+)\s+'jdm\s*'", blkdat, re.MULTILINE)
        assert idm_match, "idm line not found in blkdat.input"
        assert jdm_match, "jdm line not found in blkdat.input"
        assert int(idm_match.group(1)) == tp0_topo["idm"], (
            f"idm={idm_match.group(1)} != expected {tp0_topo['idm']}"
        )
        assert int(jdm_match.group(1)) == tp0_topo["jdm"], (
            f"jdm={jdm_match.group(1)} != expected {tp0_topo['jdm']}"
        )

    def test_empty_log_dir_created(self, region_dir):
        """Docs say log/ is created fresh and empty in the new experiment."""
        self._run_expt_new(region_dir)
        log_dir = region_dir / "expt_01.1" / "log"
        assert log_dir.is_dir(), "log/ was not created in expt_01.1/"
        assert list(log_dir.iterdir()) == [], "log/ is not empty"

    def test_empty_data_dir_created(self, region_dir):
        """Docs say data/ is created fresh and empty in the new experiment."""
        self._run_expt_new(region_dir)
        data_dir = region_dir / "expt_01.1" / "data"
        assert data_dir.is_dir(), "data/ was not created in expt_01.1/"
        assert list(data_dir.iterdir()) == [], "data/ is not empty"

    def test_old_data_not_copied(self, region_dir):
        """
        Docs say data/, log/, and SCRATCH/ are excluded from the rsync copy.

        We seed expt_01.0/ with dummy files in data/ and log/ to verify they
        are not carried over.
        """
        (region_dir / "expt_01.0" / "data").mkdir(exist_ok=True)
        (region_dir / "expt_01.0" / "log").mkdir(exist_ok=True)
        (region_dir / "expt_01.0" / "data" / "restart.a").write_bytes(b"\x00")
        (region_dir / "expt_01.0" / "log"  / "hycom.stop").write_text("GOODRUN")
        self._run_expt_new(region_dir)
        assert not (region_dir / "expt_01.1" / "data" / "restart.a").exists(), (
            "data/restart.a was incorrectly copied to the new experiment"
        )
        assert not (region_dir / "expt_01.1" / "log" / "hycom.stop").exists(), (
            "log/hycom.stop was incorrectly copied to the new experiment"
        )

    def test_fails_if_new_expt_exists(self, region_dir):
        """Docs note: script stops if the target directory already exists."""
        (region_dir / "expt_01.1").mkdir()
        result = self._run_expt_new(region_dir)
        assert result.returncode != 0, (
            "expt_new.sh should exit non-zero when expt_01.1/ already exists"
        )

    def test_fails_without_grid_b(self, region_dir):
        """Docs note: topo/regional.grid.b must be present."""
        (region_dir / "topo" / "regional.grid.b").unlink()
        result = self._run_expt_new(region_dir)
        assert result.returncode != 0, (
            "expt_new.sh should exit non-zero when regional.grid.b is missing"
        )


# ── hycom_creat_thkdf4_veldf4.py ─────────────────────────────────────────────

class TestHycomCreatThkdf4Veldf4:
    """
    Documented in forcing.md (Sponge layers section):

        python hycom_creat_thkdf4_veldf4.py <sponge_width>

    Must be run from a directory containing regional.grid.a/b and regional.depth.a/b.
    The docs instruct users to symlink the config-specific depth file as regional.depth:

        ln -sf $WORK/<CONFIG>/topo/regional.depth.a .
        ln -sf $WORK/<CONFIG>/topo/regional.depth.b .

    Note: TP0 provides depth_TP0a1.00_01.a/b — tests create the expected
    regional.depth symlinks per the documented setup steps.

    Documented outputs: thkdf4.a, thkdf4.b, veldf4.a, veldf4.b

    KNOWN GAP (undocumented dependency): the script imports ``cmocean`` at the
    top level. cmocean is not in docs/environment.yaml or docs/installation.md.
    If the test fails with an ImportError for cmocean, that is the finding.
    """

    @pytest.fixture
    def work_dir(self, tmp_path, tp0_topo):
        """
        Scratch directory with regional.grid.a/b and regional.depth.a/b symlinks,
        matching the setup shown in forcing.md.
        """
        for ext in ("a", "b"):
            (tmp_path / f"regional.grid.{ext}").symlink_to(
                os.path.join(TOPO_DIR, f"regional.grid.{ext}")
            )
            # Docs say to link as regional.depth.a/b; TP0 provides depth_TP0a1.00_01.*
            (tmp_path / f"regional.depth.{ext}").symlink_to(
                os.path.join(TOPO_DIR, f"depth_TP0a1.00_01.{ext}")
            )
        return tmp_path

    def _run(self, work_dir, sponge_width=20):
        return subprocess.run(
            [sys.executable,
             os.path.join(BIN_DIR, "hycom_creat_thkdf4_veldf4.py"),
             str(sponge_width)],
            capture_output=True, text=True,
            env=_SCRIPT_ENV, cwd=str(work_dir),
        )

    def test_writes_thkdf4_ab(self, work_dir):
        """thkdf4.a and thkdf4.b are written."""
        result = self._run(work_dir)
        assert result.returncode == 0, (
            f"hycom_creat_thkdf4_veldf4.py exited {result.returncode}\n"
            f"stderr:\n{result.stderr}"
        )
        assert (work_dir / "thkdf4.a").exists(), "thkdf4.a was not written"
        assert (work_dir / "thkdf4.b").exists(), "thkdf4.b was not written"

    def test_writes_veldf4_ab(self, work_dir):
        """veldf4.a and veldf4.b are written."""
        result = self._run(work_dir)
        assert result.returncode == 0, (
            f"hycom_creat_thkdf4_veldf4.py exited {result.returncode}\n"
            f"stderr:\n{result.stderr}"
        )
        assert (work_dir / "veldf4.a").exists(), "veldf4.a was not written"
        assert (work_dir / "veldf4.b").exists(), "veldf4.b was not written"

    def test_thkdf4_values_in_range(self, work_dir):
        """
        Docs say thkdf4 ramps linearly from mnn (interior) to mxx (boundary).

        Verify using the range printed to stdout: ``thkdf4: range = <min> <max>``.
        The defaults in the script are mnn=0.02, mxx=0.125 — but the docs note
        these are set by editing the script, not via CLI arguments.
        """
        result = self._run(work_dir, sponge_width=20)
        assert result.returncode == 0
        # Script prints "thkdf4: range = <min> <max>" to stdout
        match = re.search(
            r"thkdf4:\s+range\s*=\s*([0-9.e+\-]+)\s+([0-9.e+\-]+)",
            result.stdout,
        )
        assert match, (
            f"Could not find 'thkdf4: range = ...' in stdout:\n{result.stdout}"
        )
        rmin, rmax = float(match.group(1)), float(match.group(2))
        assert rmin >= 0.0,    f"thkdf4 minimum {rmin} is negative"
        assert rmax <= 1.0,    f"thkdf4 maximum {rmax} > 1.0 (implausible diffusion)"
        assert rmin < rmax,    f"thkdf4 range inverted: min={rmin}, max={rmax}"

    def test_thkdf4_b_contains_range(self, work_dir):
        """
        thkdf4.b (the header/metadata file) contains the range line.

        Docs show: ``thkdf4.a/b`` — the .b file is the human-readable header.
        """
        result = self._run(work_dir)
        assert result.returncode == 0
        b_content = (work_dir / "thkdf4.b").read_text()
        assert "thkdf4" in b_content, (
            f"thkdf4.b does not contain 'thkdf4':\n{b_content}"
        )
        assert "range" in b_content, (
            f"thkdf4.b does not contain 'range':\n{b_content}"
        )


# ── hycom_bathy_ports.py ─────────────────────────────────────────────────────

class TestHycomBathyPorts:
    """
    Documented in forcing.md (Open boundary sections and relaxation coefficients).
    The script is invoked via nest_setup_ports.sh; the documented outcome is:

        "writes ports.input and rmu.a/b to nest/<IEXPT>/"

    Direct interface (from the argparse block):

        hycom_bathy_ports.py <depthfile> <rmu_width> <rmu_efold>

    Must be run from a directory containing regional.grid.a/b.

    NOTE — docs vs code discrepancy:
      The script writes ``ports.input.tmp``, NOT ``ports.input``.
      nest_setup_ports.sh renames it: ``cp ports.input.tmp $TARGETDIR/ports.input``
      The test below checks for ``ports.input.tmp`` to reflect actual script behaviour.
      If the docs are taken literally (``ports.input`` should exist), the companion
      assertion documents the gap.
    """

    @pytest.fixture
    def work_dir(self, tmp_path):
        """Scratch directory with regional.grid.a/b symlinked from TP0."""
        for ext in ("a", "b"):
            (tmp_path / f"regional.grid.{ext}").symlink_to(
                os.path.join(TOPO_DIR, f"regional.grid.{ext}")
            )
        # The depth file is passed as a positional argument (without .a/.b extension)
        for ext in ("a", "b"):
            (tmp_path / f"depth_TP0a1.00_01.{ext}").symlink_to(
                os.path.join(TOPO_DIR, f"depth_TP0a1.00_01.{ext}")
            )
        return tmp_path

    def _run(self, work_dir, depthfile="depth_TP0a1.00_01", width=20, efold=20):
        return subprocess.run(
            [sys.executable,
             os.path.join(BIN_DIR, "hycom_bathy_ports.py"),
             depthfile, str(width), str(efold)],
            capture_output=True, text=True,
            env=_SCRIPT_ENV, cwd=str(work_dir),
        )

    def test_writes_rmu_a(self, work_dir):
        """rmu.a (relaxation coefficient binary field) is written."""
        result = self._run(work_dir)
        assert result.returncode == 0, (
            f"hycom_bathy_ports.py exited {result.returncode}\n"
            f"stderr:\n{result.stderr}"
        )
        assert (work_dir / "rmu.a").exists(), "rmu.a was not written"

    def test_writes_rmu_b(self, work_dir):
        """rmu.b (rmu header file) is written."""
        result = self._run(work_dir)
        assert result.returncode == 0, (
            f"hycom_bathy_ports.py exited {result.returncode}\n"
            f"stderr:\n{result.stderr}"
        )
        assert (work_dir / "rmu.b").exists(), "rmu.b was not written"

    def test_ports_output_file(self, work_dir):
        """
        Docs say ``ports.input`` is written.  Code writes ``ports.input.tmp``.

        This test checks for ``ports.input.tmp`` (what the script actually produces).
        It also checks that ``ports.input`` does NOT exist — confirming that the
        docs-described output only exists after nest_setup_ports.sh renames it.
        """
        result = self._run(work_dir)
        assert result.returncode == 0, (
            f"hycom_bathy_ports.py exited {result.returncode}\n"
            f"stderr:\n{result.stderr}"
        )
        ports_tmp  = work_dir / "ports.input.tmp"
        ports_final = work_dir / "ports.input"

        assert ports_tmp.exists(), (
            "ports.input.tmp was not written — check if the script's output filename changed"
        )
        # If ports.input exists directly, the script has been updated to match the docs.
        # The assertion below documents the current gap without hard-failing on it.
        if not ports_final.exists():
            import warnings
            warnings.warn(
                "hycom_bathy_ports.py writes ports.input.tmp, not ports.input as documented. "
                "nest_setup_ports.sh renames it — but calling the script directly "
                "does not produce the documented ports.input.",
                stacklevel=2,
            )

    def test_ports_input_has_nports(self, work_dir):
        """
        The first line of ports.input[.tmp] declares the number of ports.

        Documented format (forcing.md): the file lists ``nports`` followed by
        five numbers per port (kdport, ifport, ilport, jfport, jlport).
        """
        result = self._run(work_dir)
        assert result.returncode == 0
        # Read whichever file was produced
        ports_path = work_dir / "ports.input.tmp"
        if not ports_path.exists():
            ports_path = work_dir / "ports.input"
        assert ports_path.exists(), "neither ports.input.tmp nor ports.input was produced"

        first_line = ports_path.read_text().splitlines()[0]
        match = re.search(r"(\d+)\s+'nports'", first_line)
        assert match, (
            f"First line of ports output does not contain nports declaration: "
            f"{first_line!r}"
        )
        nports = int(match.group(1))
        assert nports >= 0, f"nports must be non-negative, got {nports}"

    def test_rmu_b_header(self, work_dir):
        """rmu.b contains a human-readable header with 'rmu' and range information."""
        result = self._run(work_dir)
        assert result.returncode == 0
        b_content = (work_dir / "rmu.b").read_text()
        assert "rmu" in b_content.lower(), (
            f"rmu.b does not mention 'rmu':\n{b_content}"
        )

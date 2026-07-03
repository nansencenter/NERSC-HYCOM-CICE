"""
Shared fixtures for the NERSC-HYCOM-CICE test suite.

Environment detection
---------------------
* Native Linux with conda-forge ESMF (GitHub Actions, Betzy via native build):
  compile_model.sh is invoked directly.
* macOS with Podman + ``hycom-esmf-spike`` image (local developer Mac):
  compile_model.sh is invoked inside the container via ``podman run``.

To build the Podman image (first time only on Mac):
    CONTAINERS_MACHINE_PROVIDER=applehv podman machine init --now
    podman build -t hycom-esmf-spike docs/testing-esmf-spike-notes.md
    # or follow the Dockerfile in docs/testing-esmf-spike-notes.md
"""

import os
import platform
import shutil
import subprocess
import textwrap

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TP0_DIR   = os.path.join(REPO_ROOT, "TP0a1.00")
TOPO_DIR  = os.path.join(TP0_DIR, "topo")
EXPT_DIR  = os.path.join(TP0_DIR, "expt_01.0")

# ── build script run inside the container (Mac) ──────────────────────────────

_CONTAINER_BUILD_SCRIPT = textwrap.dedent("""\
    #!/bin/bash
    set -e
    set -x
    source /opt/conda/etc/profile.d/conda.sh 2>/dev/null

    mkdir -p /work
    for d in TP0a1.00 hycom cice bin input pythonlibs; do
      rsync -a --exclude='.git' /repo/$d/ /work/$d/
    done

    cat > /work/TP0a1.00/REGION.src <<'INNER'
    #!/bin/bash
    export R=TP0a1.00
    mydir=$(cd $(dirname ${BASH_SOURCE}) && pwd)
    tmp=$(basename $mydir)
    if [ "$tmp" != "${R}" ]; then
       echo "Error: Mismatch between path of region $R and path with REGION.src file:${tmp}"
       exit 1
    fi
    export NHCROOT=/work
    export HYCOM_ALL=$NHCROOT/hycom/hycom_ALL/hycom_2.2.72_ALL/
    export INPUTDIR=$NHCROOT/input/
    export BINDIR=$NHCROOT/bin/
    export HYCOM_PYTHON_ROUTINES=$BINDIR
    export MSCPROGS=$NHCROOT/hycom/MSCPROGS/
    INNER

    ln -sf /work/bin /work/TP0a1.00/bin || true

    export ESMF_DIR=/opt/conda
    export ESMF_MOD_DIR=/opt/conda/mod
    export ESMF_LIB_DIR=/opt/conda/lib

    cd /work/TP0a1.00/expt_01.0
    bash /work/bin/compile_model.sh -m openmpi gfortran
""")

# ── environment detection ─────────────────────────────────────────────────────

def _find_esmf_mk():
    """Return path to esmf.mk, or None if not found."""
    if "ESMFMKFILE" in os.environ:
        return os.environ["ESMFMKFILE"]
    candidates = [
        "/opt/conda/lib/esmf.mk",
    ]
    try:
        conda_base = subprocess.check_output(
            ["conda", "info", "--base"], text=True, stderr=subprocess.DEVNULL
        ).strip()
        candidates.insert(0, os.path.join(conda_base, "lib", "esmf.mk"))
    except Exception:
        pass
    return next((p for p in candidates if os.path.exists(p)), None)


def _podman_image_ok():
    if not shutil.which("podman"):
        return False
    r = subprocess.run(
        ["podman", "image", "exists", "hycom-esmf-spike"],
        capture_output=True,
    )
    return r.returncode == 0


def _build_env():
    """Return 'native', 'podman', or None."""
    if platform.system() == "Linux" and _find_esmf_mk():
        return "native"
    if platform.system() == "Darwin" and _podman_image_ok():
        return "podman"
    return None


# ── fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture(scope="session")
def tp0_topo():
    """Paths to committed TP0 topo files (always in the repo)."""
    return {
        "grid":      os.path.join(TOPO_DIR, "regional.grid"),
        "bathy":     os.path.join(TOPO_DIR, "depth_TP0a1.00_01"),
        "cice_grid": os.path.join(TOPO_DIR, "cice_grid.nc"),
        "idm": 100,
        "jdm": 110,
    }


@pytest.fixture(scope="session")
def tp0_blkdat():
    """Path to blkdat.input for TP0 expt_01.0."""
    return os.path.join(EXPT_DIR, "blkdat.input")


@pytest.fixture(scope="session")
def tp0_ice_in():
    """Path to CICE ice_in namelist for TP0 expt_01.0."""
    return os.path.join(EXPT_DIR, "ice_in")


@pytest.fixture(scope="session")
def compiled_hycom(tmp_path_factory):
    """
    Session fixture: compile hycom_cice for TP0 expt_01.0.

    Returns a dict ``{"success": bool, "stdout": str, "stderr": str}``.
    Skips if no build environment is available.
    """
    env_type = _build_env()
    if env_type is None:
        pytest.skip(
            "No build environment: need native Linux+conda-forge ESMF "
            "or Podman+hycom-esmf-spike image on macOS"
        )

    if env_type == "native":
        esmf_mk = _find_esmf_mk()
        esmf_prefix = os.path.dirname(os.path.dirname(esmf_mk))
        env = os.environ.copy()
        env.setdefault("ESMF_DIR",     esmf_prefix)
        env.setdefault("ESMF_MOD_DIR", os.path.join(esmf_prefix, "mod"))
        env.setdefault("ESMF_LIB_DIR", os.path.join(esmf_prefix, "lib"))
        result = subprocess.run(
            ["bash",
             os.path.join(REPO_ROOT, "bin", "compile_model.sh"),
             "-m", "openmpi", "gfortran"],
            cwd=EXPT_DIR,
            env=env,
            capture_output=True,
            text=True,
            timeout=1800,
        )
        return {
            "success": result.returncode == 0,
            "stdout":  result.stdout,
            "stderr":  result.stderr,
        }

    # podman
    script_dir = tmp_path_factory.mktemp("build_script")
    script_path = script_dir / "setup_and_build.sh"
    script_path.write_text(_CONTAINER_BUILD_SCRIPT)
    result = subprocess.run(
        [
            "podman", "run", "--rm",
            "-v", f"{REPO_ROOT}:/repo:ro",
            "-v", f"{script_path}:/setup_and_build.sh:ro",
            "hycom-esmf-spike",
            "bash", "/setup_and_build.sh",
        ],
        capture_output=True,
        text=True,
        timeout=1800,
    )
    return {
        "success": result.returncode == 0,
        "stdout":  result.stdout,
        "stderr":  result.stderr,
    }

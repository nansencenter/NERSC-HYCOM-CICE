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

import numpy as np
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


# ── GLORYS synthetic fixtures ─────────────────────────────────────────────────
#
# Three files matching the real GLORYS12/Mercator format:
#   coord file  — e3t (1-D), e1t, e2t, lon, lat, depth
#   mask file   — mask (byte, depth×lat×lon), deptho, deptho_lev
#   physics     — thetao, so, uo, vo, zos (float), time ("hours since 1950-01-01")
#
# Grid: 6 lon × 5 lat (all ≥ 30 °N) × 5 depth levels so the minlat=30 filter
# in nemo2archvz_regular.py keeps all points.

_G_NLO = 6
_G_NLA = 5
_G_NLV = 5
_G_LON = np.array([0.0, 0.0833, 0.1667, 0.25,   0.3333, 0.4167], dtype="f4")
_G_LAT = np.array([30.0, 30.0833, 30.1667, 30.25, 30.3333],       dtype="f4")
_G_DEP = np.array([1.0, 10.0, 50.0, 100.0, 200.0],                dtype="f4")
_G_E3T = np.array([1.0,  9.0, 40.0,  50.0, 100.0],                dtype="f4")
_G_FILL = np.float32(9.96921e+36)


def _write_glorys_coord(path: str) -> None:
    """Write a GLORYS-format coordinate file with lon/lat/depth/e3t/e1t/e2t."""
    import netCDF4
    with netCDF4.Dataset(path, "w") as ds:
        ds.Conventions = "CF-1.6"
        ds.source = "MERCATOR GLORYS12V1"
        ds.createDimension("longitude", _G_NLO)
        ds.createDimension("latitude",  _G_NLA)
        ds.createDimension("depth",     _G_NLV)
        for name, data, attrs in [
            ("longitude", _G_LON, {"units": "degrees_east",  "standard_name": "longitude", "axis": "X"}),
            ("latitude",  _G_LAT, {"units": "degrees_north", "standard_name": "latitude",  "axis": "Y"}),
            ("depth",     _G_DEP, {"units": "m", "positive": "down", "standard_name": "depth", "axis": "Z"}),
        ]:
            v = ds.createVariable(name, "f4", (name,), fill_value=_G_FILL)
            v[:] = data
            for k, val in attrs.items():
                setattr(v, k, val)
        v = ds.createVariable("e3t", "f4", ("depth",), fill_value=_G_FILL)
        v[:] = _G_E3T
        v.units = "m"; v.standard_name = "cell_thickness"; v.long_name = "Cell dimension along Z axis"
        for name in ("e1t", "e2t"):
            v = ds.createVariable(name, "f4", ("latitude", "longitude"), fill_value=_G_FILL)
            v[:] = np.full((_G_NLA, _G_NLO), 9260.0, dtype="f4")
            v.units = "m"


def _write_glorys_mask(path: str) -> None:
    """Write a GLORYS-format mask/bathy file with mask, deptho, and deptho_lev."""
    import netCDF4
    with netCDF4.Dataset(path, "w") as ds:
        ds.Conventions = "CF-1.6"
        ds.source = "GLORYS12V1"
        ds.createDimension("depth",     _G_NLV)
        ds.createDimension("latitude",  _G_NLA)
        ds.createDimension("longitude", _G_NLO)
        for name, data, dims in [
            ("longitude", _G_LON, ("longitude",)),
            ("latitude",  _G_LAT, ("latitude",)),
            ("depth",     _G_DEP, ("depth",)),
        ]:
            v = ds.createVariable(name, "f4", dims)
            v[:] = data
        v = ds.createVariable("mask", "i1", ("depth", "latitude", "longitude"))
        v[:] = np.ones((_G_NLV, _G_NLA, _G_NLO), dtype="i1")
        v.long_name = "Land-sea mask: 1 = sea ; 0 = land"
        v.standard_name = "sea_binary_mask"; v.units = "1"
        v = ds.createVariable("deptho", "f4", ("latitude", "longitude"), fill_value=_G_FILL)
        v[:] = np.full((_G_NLA, _G_NLO), 200.0, dtype="f4")
        v.units = "m"; v.standard_name = "sea_floor_depth_below_geoid"
        v = ds.createVariable("deptho_lev", "f4", ("latitude", "longitude"), fill_value=_G_FILL)
        v[:] = np.full((_G_NLA, _G_NLO), float(_G_NLV), dtype="f4")
        v.long_name = "Model level number at sea floor"; v.units = "1"


def _write_glorys_physics(path: str) -> None:
    """Write a GLORYS-format physics file with thetao, so, uo, vo, zos.

    Velocities have a vertical gradient so the baroclinic component is non-zero;
    nemo2archvz_regular.py raises RuntimeError if baroclinic velocity is all-zero.
    """
    import netCDF4
    # 1993-01-01 12:00 UTC in hours since 1950-01-01
    # Leap years 1950–1992: 1952,56,60,64,68,72,76,80,84,88,92 = 11
    # Days = 43*365 + 11 = 15706; +0.5 day = 15706.5 → ×24 = 376956 h
    t_val = np.float32(376956.0)
    with netCDF4.Dataset(path, "w") as ds:
        ds.Conventions = "CF-1.4"; ds.source = "MERCATOR GLORYS12V1"
        ds.createDimension("time",      1)
        ds.createDimension("depth",     _G_NLV)
        ds.createDimension("latitude",  _G_NLA)
        ds.createDimension("longitude", _G_NLO)
        v = ds.createVariable("time", "f4", ("time",))
        v[:] = [t_val]; v.units = "hours since 1950-01-01"; v.calendar = "gregorian"; v.axis = "T"
        for name, data, dims in [
            ("depth",     _G_DEP, ("depth",)),
            ("latitude",  _G_LAT, ("latitude",)),
            ("longitude", _G_LON, ("longitude",)),
        ]:
            v = ds.createVariable(name, "f4", (name,))
            v[:] = data
        s4 = (1, _G_NLV, _G_NLA, _G_NLO)
        for name, val, units, sname in [
            ("thetao", 10.0, "degrees_C", "sea_water_potential_temperature"),
            ("so",     35.0, "1e-3",      "sea_water_salinity"),
        ]:
            v = ds.createVariable(name, "f4", ("time","depth","latitude","longitude"), fill_value=_G_FILL)
            v[:] = np.full(s4, val, dtype="f4")
            v.units = units; v.standard_name = sname
        # Velocities must decrease with depth so the baroclinic component (layer minus
        # depth-average) is non-zero — the script raises RuntimeError if ul is all-zero.
        for name, base, units, sname in [
            ("uo", 0.10, "m s-1", "eastward_sea_water_velocity"),
            ("vo", 0.04, "m s-1", "northward_sea_water_velocity"),
        ]:
            data = np.zeros(s4, dtype="f4")
            for k in range(_G_NLV):
                data[0, k, :, :] = base * (1.0 - k / _G_NLV)
            v = ds.createVariable(name, "f4", ("time","depth","latitude","longitude"), fill_value=_G_FILL)
            v[:] = data
            v.units = units; v.standard_name = sname
        v = ds.createVariable("zos", "f4", ("time","latitude","longitude"), fill_value=_G_FILL)
        v[:] = np.zeros((1, _G_NLA, _G_NLO), dtype="f4")
        v.units = "m"; v.standard_name = "sea_surface_height_above_geoid"


@pytest.fixture(scope="session")
def glorys_fixtures(tmp_path_factory):
    """Synthetic GLORYS coord, mask/bathy, and physics NetCDF files (5×5×5, lat ≥ 30 °N)."""
    pytest.importorskip("netCDF4", reason="netCDF4 required for GLORYS preprocessing tests")
    d = tmp_path_factory.mktemp("glorys")
    coord = str(d / "GLO_MFC_001_MESH.nc")
    mask  = str(d / "GLO_MFC_001_mask_bathy.nc")
    # Filename must match MERCATOR-PHY-24-YYYY-MM-DD-12.nc
    # The script slices file[-16:-6] to extract the date "1993-01-01"
    phys  = str(d / "MERCATOR-PHY-24-1993-01-01-12.nc")
    _write_glorys_coord(coord)
    _write_glorys_mask(mask)
    _write_glorys_physics(phys)
    return {"coord": coord, "mask": mask, "physics": phys}


# ── ERA5 synthetic fixtures ───────────────────────────────────────────────────
#
# One file per variable, named 6h.VARNAME_%Y.nc, matching the newer (2025-style)
# format: float data, NaNf fill, double lat/lon/time, "hours since 1900-01-01".
# Grid: 1440 × 721 global 0.25° (needed for interpolation onto the TP0 HYCOM grid).
# Time: 4 steps covering 1993-01-01 00:00–18:00 UTC (6-hourly).
#
# Days from 1900-01-01 to 1993-01-01:
#   93 years; leap years 1904–1992 (not 1900) = 23
#   93×365 + 23 = 33968 days → 815232 h since 1900-01-01

_E5_NLO  = 1440
_E5_NLA  = 721
_E5_LON  = np.linspace(0.0, 359.75, _E5_NLO, dtype="f8")
_E5_LAT  = np.linspace(90.0, -90.0, _E5_NLA, dtype="f8")   # north→south
_E5_T0   = 815232.0  # hours since 1900-01-01 at 1993-01-01 00:00 UTC
_E5_TIMES = np.array([_E5_T0, _E5_T0+6, _E5_T0+12, _E5_T0+18], dtype="f8")

# (known_name, nc_varname, default_value, units, accumulated)
# Fields mirror input/era5.xml and era5+all: only names listed in
# atmosphere.py::_all_known_names are valid. Excluded:
#   blh, skt — not in _all_known_names (AtmosphericForcingError: Unknown field)
#   ci        — not in any production ERA5 config
#   ro        — assumed unit "1" is not convertible from "m s**-1" (after accumulation)
_ERA5_VARS = [
    ("10u",  "10U",  5.0,    "m s**-1",  False),
    ("10v",  "10V",  2.0,    "m s**-1",  False),
    ("2t",   "2T",   275.0,  "K",        False),
    ("2d",   "2D",   270.0,  "K",        False),
    ("msl",  "MSL",  101325.0, "Pa",     False),
    ("tcc",  "TCC",  0.5,    "1",        False),
    ("ssrd", "SSRD", 1.0e6,  "J m**-2",  True),
    ("strd", "STRD", 1.2e6,  "J m**-2",  True),
    ("ssr",  "SSR",  8.0e5,  "J m**-2",  True),
    ("str",  "STR",  4.0e5,  "J m**-2",  True),
    ("tp",   "TP",   0.001,  "m",        True),
]


def _write_era5_file(path: str, nc_varname: str, val: float, units: str) -> None:
    """Write a single ERA5 6-hourly file in the newer (2025-style) float format."""
    import netCDF4
    nt = len(_E5_TIMES)
    with netCDF4.Dataset(path, "w") as ds:
        ds.Conventions = "CF-1.7"
        ds.institution = "European Centre for Medium-Range Weather Forecasts"
        ds.createDimension("time",      nt)
        ds.createDimension("latitude",  _E5_NLA)
        ds.createDimension("longitude", _E5_NLO)
        v = ds.createVariable("time", "f8", ("time",))
        v[:] = _E5_TIMES
        v.units = "hours since 1900-01-01 00:00:00"; v.calendar = "gregorian"; v.axis = "T"
        v = ds.createVariable("latitude", "f8", ("latitude",))
        v[:] = _E5_LAT; v.units = "degrees_north"; v.axis = "Y"
        v = ds.createVariable("longitude", "f8", ("longitude",))
        v[:] = _E5_LON; v.units = "degrees_east"; v.axis = "X"
        v = ds.createVariable(nc_varname, "f4", ("time","latitude","longitude"),
                              fill_value=np.float32(float("nan")))
        v[:] = np.full((nt, _E5_NLA, _E5_NLO), val, dtype="f4")
        v.units = units; v.long_name = nc_varname


def _write_era5_xml(path: str, era5_dir: str) -> None:
    """Write an XML config for hycom_atmfor.py pointing to synthetic ERA5 files.

    atmosphere.py uses XPath 'forcing_datasets/forcing_dataset[@name="..."]' on an
    ElementTree, so <forcing_datasets> must be a child of the root element, not
    the root itself.
    """
    fields = []
    for known, ncvar, _, units, accum in _ERA5_VARS:
        acc_attr = ' accumulated="6h"' if accum else ''
        fields.append(
            f'    <field known_name="{known}" '
            f'file="{era5_dir}/6h.{ncvar}_%Y.nc" '
            f'varname="{ncvar}" '
            f'units="{units}"{acc_attr}/>'
        )
    xml = (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<config>\n'
        '<forcing_datasets>\n'
        '  <forcing_dataset name="ERA5_test" rootPath="" timestep="6h" format="netcdf">\n'
        + "\n".join(fields) + "\n"
        '  </forcing_dataset>\n'
        '</forcing_datasets>\n'
        '</config>\n'
    )
    with open(path, "w") as f:
        f.write(xml)


@pytest.fixture(scope="session")
def era5_fixtures(tmp_path_factory):
    """
    Synthetic ERA5 NetCDF files (one per variable) + XML config for hycom_atmfor.py.

    Returns dict with keys:
      "dir"  — directory containing 6h.VARNAME_1993.nc files
      "xml"  — path to the XML config file
    """
    pytest.importorskip("netCDF4", reason="netCDF4 required for ERA5 preprocessing tests")
    d = tmp_path_factory.mktemp("era5")
    era5_dir = str(d)
    for known, ncvar, val, units, _ in _ERA5_VARS:
        _write_era5_file(str(d / f"6h.{ncvar}_1993.nc"), ncvar, val, units)
    xml_path = str(d / "era5_forcing.xml")
    _write_era5_xml(xml_path, era5_dir)
    return {"dir": era5_dir, "xml": xml_path}

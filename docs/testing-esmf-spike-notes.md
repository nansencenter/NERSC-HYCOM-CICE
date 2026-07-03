# Coupled build (ESMF/CICE) test spike — handoff notes

Context: exploring whether the ESMF/CICE-coupled HYCOM-CICE build (`iceflg=2`,
TP0 `expt_01.0`) can be compiled and tested locally (Mac) and in GitHub Actions,
not just on Betzy/Olivia. This document records what was found so a follow-up
session (or a session on Betzy) can pick up without repeating the investigation.

## Goal

Design a `tests/` suite (modeled on `../xhycom/tests` and
[ucla-roms](https://github.com/CWorthy-ocean/ucla-roms)'s pytest-driven
compile+run+hash-check pattern) covering:

1. Pre-processing scripts (`bin/`, `pythonlibs/`)
2. The HYCOM/CICE Fortran source, run against real TP0 data
3. MSCPROGS post-processing

Priority: 1 and 2, always against real TP0 data (no synthetic fixtures, no
skip-on-missing-data), runnable both locally (Mac) and in GitHub Actions CI.

## Key findings

### `bin/compile_model.sh` is Linux-only
Line ~112 hard-exits (`exit 2`) if `uname -s` isn't `Linux`. Native macOS is a
non-starter without patching this script. **Decision: use Linux containers
instead** (Podman, since Docker Desktop has commercial licensing restrictions;
podman-py on PyPI is only a client SDK, not a container engine — the actual
Podman binary + a Linux VM is required on macOS either way).

### hycom_ALL and CICE are already portable to plain gfortran
- `hycom/hycom_ALL/hycom_2.2.72_ALL/config/gfortran_setup` and `macg5_setup`
  exist — no Intel compiler needed for pre/post-processing utilities.
- `cice/Release-5.1/bld/Macros.Linux.gfortran.openmpi` and
  `alt_macros/Macros.Darwin.CPOM.mi2` exist — CICE 5.1 builds with gfortran too.
- No existing HYCOM `RELO/config` file combines gfortran + OpenMPI + ESMF
  *without* OASIS/FABM baggage for a clean standalone-ESMF build; one may need
  to be added (see `Linux.gfortran.openmpi_cice`, which is close but should be
  checked/adapted).

### TP0 `expt_01.0` needs the coupled path, not standalone HYCOM
`TP0a1.00/expt_01.0/blkdat.input`: `iceflg=2` (coupled/ESMF), `ntracr=0` (no
BGC/FABM needed). So the default TP0 experiment genuinely exercises the
ESMF/CICE coupling — testing standalone HYCOM (`iceflg=0`) alone would not
cover what TP0 actually runs.

### ESMF is available prebuilt via conda-forge — no need to build it from source
`conda install -c conda-forge gfortran openmpi esmf netcdf-fortran` gives a
working `esmf.mk` (confirmed: ESMF 8.9.1, linux-aarch64 and presumably
linux-64/osx builds too). This is a huge simplification vs. HPC-style
`module load` + hand-built ESMF.

### Podman setup on macOS (for local testing)
- Docker Desktop has licensing restrictions for larger commercial use;
  **Podman** (Apache 2.0, `brew install podman`) was used instead as a
  license-free, Docker-CLI-compatible alternative.
- `podman machine init` defaults to the `libkrun` provider, which needs the
  `krunkit` binary from the `slp/krunkit` Homebrew tap — that tap is
  "untrusted" and requires explicit `brew trust`, which was avoided.
- **Fix: use the `applehv` provider instead** (native macOS virtualization,
  no third-party binary):
  ```bash
  CONTAINERS_MACHINE_PROVIDER=applehv podman machine init --now
  ```
  VM image pull is from `quay.io/podman/machine-os:6.0`, ~1.2GB — slow on this
  network (~20-30 min). Budget for this in CI-adjacent setup instructions but
  note GitHub Actions runners are native Linux and skip this VM step entirely.

### Working spike Dockerfile
```dockerfile
FROM condaforge/miniforge3:latest

RUN apt-get update && apt-get install -y csh gawk rsync && rm -rf /var/lib/apt/lists/*

RUN conda install -y -c conda-forge gfortran openmpi esmf netcdf-fortran make \
    && conda clean -afy

ENV CONDA_PREFIX=/opt/conda
ENV ESMFMKFILE=/opt/conda/lib/esmf.mk
ENV PATH=/opt/conda/bin:$PATH

WORKDIR /work
```
Gotcha: `csh` is not on conda-forge for linux-aarch64 — install it via `apt`
instead (works fine, it's just tcsh under Debian/Ubuntu).

### ESMF 5-era coupling code vs. ESMF 8.9 — probe results
Only two files in `hycom/RELO/src_2.2.98ZA-07Tsig0-i-sm-sse_relo_mpi/` use
ESMF directly: `mod_esmf_utils.F` and `mod_OICPL.F` (both guarded by
`#if (USE_ESMF)` / `#if (USE_ESMF_5)`).

Compiled each standalone against conda-forge ESMF 8.9 with:
```bash
cpp -P -DUSE_ESMF -DUSE_ESMF_5 -DIA32 -DREAL8 -DMPI -DSERIAL_IO -DRELO -DNERSC_HYCOM_CICE mod_X.F > /tmp/mod_X.f90
mpif90 -c /tmp/mod_X.f90 -o /tmp/mod_X.o -I$CONDA_PREFIX/mod -I$CONDA_PREFIX/include -ffixed-form -ffixed-line-length-none -J/tmp
```
(Note: these files are **fixed-form** Fortran despite the modernized
`use ESMF` syntax — must pass `-ffixed-form`, not `-ffree-form`, or gfortran
mis-parses the `&` continuations as garbage.)

- **Both files compile cleanly against ESMF 8.9 — no source changes needed.**

  **Correction to an earlier version of this note:** the first probe run
  reported `mod_esmf_utils.F` failing on two `ESMF_TimeIntervalGet(...,
  s_r8=...)` calls (lines 429, 440-441) with "no specific subroutine for the
  generic". That was a **probe bug, not a real incompatibility**: the
  production build configs (e.g. `Linux.gfortran.openmpi_cice`,
  `Linux.fram.gfortran_hycom`) always compile with
  `-fdefault-real-8 -fdefault-double-8`, which promotes plain `real`
  declarations to 8 bytes. The probe's first pass omitted those flags, so
  `ocnTimeStep_sr8`/`iceTimeStep_sr8` (declared plain `real`) were 4-byte
  reals — a kind mismatch against ESMF's `s_r8` (real*8) dummy argument,
  which is exactly the kind of thing that produces "no specific subroutine
  for the generic" in Fortran overload resolution. Re-running with
  `-fdefault-real-8 -fdefault-double-8 -fconvert=big-endian` added (matching
  the real build config) made both files compile clean on the first try.
  **Conclusion: the ESMF-5-era coupling code is source-compatible with ESMF
  8.9 as-is.** No porting work, no new CPP macro needed after all.

## Revised next steps

The main open blocker was resolved — no ESMF porting needed. What's still
unverified is a **full link+run**, not just these 2 files compiling in
isolation:

1. Attempt the full `compile_model.sh`-driven build of `hycom_cice` for TP0
   `expt_01.0` inside the Podman container (need to work around the script's
   Linux-only `uname` guard — trivial inside the container, since it reports
   `Linux` natively). This exercises the rest of the HYCOM/CICE source tree,
   the CICE build (`Macros.Linux.gfortran.openmpi`), and the final link step,
   none of which this spike has touched yet.
2. If it links, run a short TP0 integration and confirm sane output
   (`pythonlibs/abfile` can read the result).
3. Once compilation + a short run are confirmed working in the container,
   resume `tests/` suite design (session fixture running preprocessing →
   compile → run → hash-check, per the ucla-roms-inspired plan).
4. Validating on Betzy is still valuable as a second, HPC-native data point,
   but is no longer gating progress here — the user will do that separately
   later.

## Scratch artifacts from this spike (not committed)
- Dockerfile + probe script were written under this session's scratchpad dir
  (`/private/tmp/claude-502/<session-id>/scratchpad/esmf-spike/`), which is
  **session-specific and will not exist in a new session** — the exact
  contents needed to recreate it are inlined below.
- Podman machine `podman-machine-default` (applehv provider) and the built
  image `localhost/hycom-esmf-spike:latest` **do persist** across sessions
  (they live in Podman's actual VM/storage on disk, not the session
  scratchpad). A new session can go straight to `podman images` /
  `podman machine list` to check they're still there before rebuilding
  anything. If the machine was stopped, restart with
  `podman machine start`.

### STATUS AS OF SESSION END (2026-07-03): mid-way through a full build attempt

Confirmed so far (see above): both ESMF-dependent files compile clean in
isolation. **Not yet attempted: an actual full `compile_model.sh` build of
`hycom_cice` for TP0 `expt_01.0`.** A setup script for this was drafted but
**not executed** (interrupted by user to end the session cleanly). Resume
from here:

```bash
# Recreate Dockerfile if the image is gone (check `podman images` first):
cat > Dockerfile <<'EOF'
FROM condaforge/miniforge3:latest
RUN apt-get update && apt-get install -y csh gawk rsync && rm -rf /var/lib/apt/lists/*
RUN conda install -y -c conda-forge gfortran openmpi esmf netcdf-fortran make \
    && conda clean -afy
ENV CONDA_PREFIX=/opt/conda
ENV ESMFMKFILE=/opt/conda/lib/esmf.mk
ENV PATH=/opt/conda/bin:$PATH
WORKDIR /work
EOF
podman build -t hycom-esmf-spike .

# Draft build-attempt script (NOT YET RUN — this is the next step):
cat > setup_and_build.sh <<'EOF'
#!/bin/bash
set -e
set -x
source /opt/conda/etc/profile.d/conda.sh 2>/dev/null

# Copy needed source trees into a writable workspace (host repo is mounted read-only)
mkdir -p /work
for d in TP0a1.00 hycom cice bin input pythonlibs; do
  rsync -a --exclude='.git' /repo/$d/ /work/$d/
done

# REGION.src for TP0a1.00 (repo's input/REGION.src is a template for TP5;
# this is a TP0-specific copy with R and NHCROOT adjusted)
cat > /work/TP0a1.00/REGION.src <<'INNER'
#!/bin/bash
export R=TP0a1.00
mydir=$(cd $(dirname ${BASH_SOURCE}) && pwd)
tmp=$(basename $mydir)
if [ "$tmp" != "${R}" ] ;then
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

# Bypass compile_model.sh's hardcoded per-site ESMF path lookup by presetting
# these — the script uses them as-is if already set (see compile_model.sh
# ~line 163: "if ESMF_DIR && ESMF_MOD_DIR && ESMF_LIB_DIR already set, use them")
export ESMF_DIR=/opt/conda
export ESMF_MOD_DIR=/opt/conda/mod
export ESMF_LIB_DIR=/opt/conda/lib

cd /work/TP0a1.00/expt_01.0
bash /work/bin/compile_model.sh -m openmpi gfortran
EOF
chmod +x setup_and_build.sh

podman run --rm -v /Users/noroos/NERSC-HYCOM-CICE:/repo:ro \
  -v "$PWD/setup_and_build.sh:/setup_and_build.sh:ro" \
  hycom-esmf-spike bash /setup_and_build.sh
```

**Known unresolved question going into this step:** `compile_model.sh`
requires a `hycom/RELO/config/${MACROID}_cice` file to exist, where
`MACROID=Linux.gfortran.openmpi` for `-m openmpi gfortran` on a generic
(unrecognized-hostname) Linux box. The repo already has
`hycom/RELO/config/Linux.gfortran.openmpi_cice` (uses `-DUSE_ESMF -DUSE_ESMF_5`,
no OASIS/FABM) which should match, but this has **not been traced through the
CICE compile step (`comp_ice.esmf`) or the final `Make_cice.csh` link step** —
expect to hit and need to fix real errors there (missing config file variants,
CICE macro mismatches, etc.). Iterate on actual compiler/script output rather
than pre-tracing further.

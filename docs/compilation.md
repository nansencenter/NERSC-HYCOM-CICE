Before compiling, source the HPC environment file to load the correct modules and compilers.
See the [HPC environment](installation.md#hpc-environment) section for details.

::::{dropdown} Betzy (NRIS/Sigma2)

```{include} _snippets/betzy_hpc_env.md
```

::::

## BGC dependencies (MSCPROGS and FABM)

Only needed when running with biogeochemical modules (`ntracr` > 0 in `blkdat.input`).
These libraries only need to be rebuilt when the HPC modules are updated or the source changes.

### Compile MSCPROGS (libhycnersc.a)

Create a symlink for your machine and compiler, then build and install:

```bash
cd ${HOME}/NERSC-HYCOM-CICE/hycom/MSCPROGS/src/Make.Inc
ln -sf make.betzy.ifort make.inc
cd ${HOME}/NERSC-HYCOM-CICE/hycom/MSCPROGS/src
gmake clean && gmake all && gmake install
```

After installation, `libhycnersc.a` is at:

```
${HOME}/NERSC-HYCOM-CICE/hycom/MSCPROGS/lib/libhycnersc.a
```

For more detail, see the [MSCPROGS README](https://github.com/nansencenter/NERSC-HYCOM-CICE/tree/develop/hycom/MSCPROGS).

::::{dropdown} Missing build directories (NERSC-HYCOM-CICE_BIORANv2 only)

If using the `NERSC-HYCOM-CICE_BIORANv2` repository, create these missing directories
before running `gmake`:

```bash
mkdir -p ${HOME}/NERSC-HYCOM-CICE/hycom/MSCPROGS/src/Average/TMP
mkdir -p ${HOME}/NERSC-HYCOM-CICE/hycom/MSCPROGS/src/Hycom_mean/TMP
mkdir -p ${HOME}/NERSC-HYCOM-CICE/hycom/MSCPROGS/src/Perturb_Forcing-2.2.98_TP5/TMP
mkdir -p ${HOME}/NERSC-HYCOM-CICE/hycom/MSCPROGS/src/Perturb_Parameter/TMP
```

::::

### Compile FABM (libfabm.a)

Configure and build inside a `build/` directory within the FABM clone:

```bash
rm -rf ${HOME}/fabm/build
mkdir ${HOME}/fabm/build && cd ${HOME}/fabm/build
cmake ${HOME}/fabm \
    -DFABM_HOST=hycom \
    -DCMAKE_Fortran_COMPILER=ifort \
    -DFABM_INSTITUTES="ersem;nersc;gotm" \
    -DFABM_NERSC_BASE=${HOME}/nersc \
    -DFABM_ERSEM_BASE=${HOME}/ersem
make install
```

After installation, `libfabm.a` is at:

```
${HOME}/local/fabm/hycom/lib64/libfabm.a
```

## Compile HYCOM-CICE

Run the compile script from the experiment directory:

```bash
cd ${WORK}/<CONFIGNAME>/expt_<EXPT_ID>
bash ${HOME}/NERSC-HYCOM-CICE/bin/compile_model.sh ifort -u
```

The script manages the `build/` directory containing a per-experiment copy of the source:

- **First compile** (`build/` does not exist): source is automatically synced from the repository.
- **Recompile after repo changes** (use `-u`): source is resynced from the repository, overwriting any local modifications in `build/`.
- **Recompile without `-u`**: uses whatever is currently in `build/`, preserving any local edits.

To start completely fresh, delete `build/` before running the script:

```bash
rm -rf build
```

::::{dropdown} What compile_model.sh does behind the scenes

The script must be run from the experiment directory. It sources `REGION.src` and
`EXPT.src` to read configuration, then proceeds in the following order:

1. **Detect site and compiler**: derives `SITE` from the hostname and constructs a
   `MACROID` (e.g. `Linux.betzy.ifort`) used to select the correct Makefile macros.

2. **Resolve ESMF paths**: sets `ESMF_DIR`, `ESMF_MOD_DIR`, and `ESMF_LIB_DIR`.
   If these are already set in the environment they are used as-is; otherwise the script
   uses hardcoded paths for known sites (Betzy, Fram) or exits with an error on
   unknown systems.

3. **Select equation of state**: reads `SIGVER` from `EXPT.src` and links the
   matching `stmt_fns.h` include file in the build directory. This must be consistent
   with `thflag` in `blkdat.input` (the script checks this and exits if they differ).

4. **Sync source code**: rsyncs HYCOM and CICE source from `$NHCROOT/hycom/RELO/`
   into `build/` on first run, or when `-u` is passed. Existing `build/` directories
   are left untouched otherwise, preserving any local edits.

5. **Apply feature flags**: if a file named `hycom_feature_flags` exists in the
   experiment directory it is copied into the build directory and passed to the C
   preprocessor, allowing you to activate or deactivate specific HYCOM features.

6. **Compile CICE**: runs `comp_ice.esmf` with the grid dimensions read from
   `blkdat.input`. Only runs when `iceflg=2` in `blkdat.input`.

7. **Compile HYCOM-CICE**: builds the final `hycom_cice` executable and links it
   against the CICE objects.

::::

The executable `hycom_cice` is created at:

```
${WORK}/<CONFIGNAME>/expt_<EXPT_ID>/build/src_2.2.98ZA-07Tsig0-i-sm-sse_relo_mpi/hycom_cice
```

> **Note:** If compiling with the BGC module, ensure `ntracr` in `blkdat.input` is
> non-zero and copy the FABM configuration files and CICE namelist into the experiment
> directory before compiling:
>
> ```bash
> cd ${WORK}/<CONFIGNAME>/expt_<EXPT_ID>
> cp /nird/datalake/NS9481K/shuang/TP2_setup/exp02.6_seaclim_ref/fabm.yaml .
> cp /nird/datalake/NS9481K/shuang/TP2_setup/exp02.6_seaclim_ref/hycom_fabm.nml .
> cp /nird/datalake/NS9481K/shuang/TP2_setup/exp02.6_seaclim_ref/ice_in .
> ```

## Compile hycom_ALL

`hycom_ALL` provides domain-independent pre/post processing utilities used to prepare
external input files. Like MSCPROGS and FABM, it only needs to be compiled once per
machine/compiler setup, not per experiment.

::::{dropdown} What hycom_ALL compiles and what it is used for

`Make_all.com` builds executables into several subdirectories:

| Directory | Purpose |
|-----------|---------|
| `relax/` | Generate relaxation climatologies (`relaxi`, `relaxv`, `rmu`, …) |
| `topo/` | Grid and bathymetry processing |
| `archive/` | Archive file conversion (`archv2data2d`, `archt2archv`, …) |
| `force/` | Atmospheric forcing file preparation |
| `meanstd/` | Compute mean and standard deviation fields |
| `sample/` | Sample model output at points or sections |
| `subregion/` | Extract subregions from model output |
| `bin/` | General utilities, including `hycom_nest_dates` (used for nesting) |

All programs read `regional.grid.b` at runtime to get the grid dimensions for the
domain being processed. The same executables work for any HYCOM configuration.

These utilities are called by setup and pre-processing scripts such as
`create_ref_case.sh` (relaxation climatologies and MPI tile definition files) and
`nemo2hycom.sh` (boundary and initial conditions from NEMO output). Recompilation is
only needed when the HPC modules or compiler change.

::::


```bash
cd ${HOME}/NERSC-HYCOM-CICE/hycom/hycom_ALL/hycom_2.2.72_ALL
```

Open `Make_all.src` in a text editor. Find the `ARCH` line and set it to your compiler:

```
setenv ARCH intelIFC
```

Then compile:

```bash
csh ./Make_clean.com
csh ./Make_all.com
```

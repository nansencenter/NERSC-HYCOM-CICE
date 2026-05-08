## Model components

This setup couples several models and libraries to provide ocean, sea ice, and biogeochemical simulations.

1. **[NERSC-HYCOM-CICE](https://github.com/nansencenter/NERSC-HYCOM-CICE/tree/develop)** (develop branch)  
   Coupled ocean–sea ice–biogeochemistry framework developed at NERSC.

2. **[ECOSMO_operational / nersc](https://github.com/nansencenter/nersc)** (main branch)  
   Biogeochemical module optimized for operational use at NERSC. Cloned as the `nersc` package
   and used within the FABM framework.

3. **[GOTM](https://github.com/gotm-model/code)**  
   1D column ocean model. Only the light turbulence module is used, installed as part of FABM.

4. **[FABM](https://github.com/fabm-model/fabm)**  
   Fortran 2003 coupler for biogeochemical models. Provides the framework that links ECOSMO,
   ERSEM, and GOTM to HYCOM-CICE.

5. **[ERSEM](https://github.com/pmlmodelling/ersem)**  
   Marine biogeochemical and ecosystem model. Only the carbon chemistry module is used in
   this installation.

6. **NERSC Python library**  
   Included in the NERSC-HYCOM-CICE repository under `pythonlibs/`. Used for pre- and
   post-processing.

## Directory structure

The key convention is a separation between source and run directories: `HOME` holds
cloned repositories and compiled libraries, while `WORK` holds configuration files
and model output. This separation also matters for data safety: on many HPC systems the
work filesystem is subject to automatic purging, so source code, compiled binaries, and
configuration files should always live on home.

```
${HOME}/${USER}/
├── NERSC-HYCOM-CICE/
│   ├── pythonlibs/                  # NERSC Python libraries
│   ├── bin/                         # NERSC utility scripts
│   └── hycom/
│       ├── RELO/src_2.2.98.../      # HYCOM source files
│       └── MSCPROGS/
│           ├── bin/
│           └── lib/libhycnersc.a
├── ersem/                           # ERSEM source (carbon chemistry module)
├── fabm/                            # FABM source and compiled libraries
└── nersc/                           # NERSC-specific ECOSMO model source

${WORK}/${USER}/
└── <CONFIGNAME>/                    # e.g. TP2a0.10
    ├── bin -> (symlink to repo bin)
    ├── REGION.src
    ├── expt_<EXPT_ID>/              # e.g. expt_04.2
    │   ├── EXPT.src
    │   ├── SCRATCH/
    │   ├── data/
    │   └── build/src_.../hycom_cice
    ├── nest/<IEXPT>/                # e.g. nest/042
    ├── relax/<IEXPT>/
    └── force/synoptic/<IEXPT>/
```
> **Note:** This shows the finished directory structure for reference. The following sections build it up step by step.

`${HOME}` and `${USER}` are standard shell variables set automatically. `<CONFIGNAME>`,
`<EXPT_ID>`, and `<IEXPT>` are placeholders for user-defined values, with Betzy examples
in the table below.

| Variable | Description | Example (Betzy) |
|---|---|---|
| `HOME` | Root of the source/home filesystem | `/cluster/home` |
| `WORK` | Root of the run/work filesystem | `/cluster/work/users` |
| `USER` | Your username on the system | `nlo043` |
| `<CONFIGNAME>` | Grid/configuration name | `TP2a0.10` |
| `<EXPT_ID>` | Experiment identifier (dot notation) | `04.2` |
| `<IEXPT>` | Integer experiment identifier (used in directory names) | `042` |


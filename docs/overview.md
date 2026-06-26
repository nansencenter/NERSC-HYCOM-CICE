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

The roles (`HOME`, `WORK`) are the same on every machine, but the actual filesystem each
maps to is machine-specific. The examples below use Betzy; see the
[Olivia storage layout](#olivia-storage-layout) dropdown for how the same roles map onto
Olivia, where the storage areas differ in important ways.

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
:::{note}
This shows the finished directory structure for reference. The following sections build it up step by step.
:::

`${HOME}` and `${USER}` are standard shell variables set automatically. `<CONFIGNAME>`,
`<EXPT_ID>`, and `<IEXPT>` are placeholders for user-defined values, with Betzy examples
in the table below.

| Variable | Description | Example (Betzy) | Example (Olivia) |
|---|---|---|---|
| `HOME` | Root of the source/home filesystem | `/cluster/home` | `/cluster/home` |
| `WORK` | Root of the run/work filesystem | `/cluster/work/users` | `/cluster/work/projects/<PROJECT>` |
| `USER` | Your username on the system | `nlo043` | `nlo043` |
| `<PROJECT>` | Compute/storage allocation code | `nn9481k` | `nn9481k` |
| `<CONFIGNAME>` | Grid/configuration name | `TP2a0.10` | `TP2a0.10` |
| `<EXPT_ID>` | Experiment identifier (dot notation) | `04.2` | `04.2` |
| `<IEXPT>` | Integer experiment identifier (used in directory names) | `042` | `042` |

(<a name="olivia-storage-layout"></a>)
::::{dropdown} Olivia storage layout (NRIS/Sigma2)

Olivia's storage areas differ from Betzy's, so the `HOME`/`WORK` roles map onto different
filesystems and a few extra rules apply. See the
[Sigma2 storage documentation](https://documentation.sigma2.no/files_storage/clusters.html)
for the authoritative description.

| Area | Path | Role | Backup | Notes |
|---|---|---|---|---|
| Home | `/cluster/home/$USER` | small personal files, repos | yes | Only ~20 GiB — too small for the Python/container environment or build artifacts |
| Project | `/cluster/projects/<PROJECT>` | persistent: env, compiled binaries, archived output | no | Large quota, **not** purged — the long-term home for anything you want to keep |
| Work | `/cluster/work/projects/<PROJECT>` | `WORK` — active runs (`SCRATCH/`, `data/`) | no | Fast Lustre, shared per project, **auto-deleted after 21 days** |
| NIRD | `/nird/datalake/<NS-PROJECT>`, `/nird/datapeak/<NS-PROJECT>` | shared archive / collaborator data | — | Mounted on service (SVC) nodes (read-write) and compute nodes (**read-only**); **not** mounted on login nodes |

Two consequences for the layout above:

1. **There is no per-user work area.** Unlike Betzy's `/cluster/work/users/$USER`, Olivia
   provides only the per-project `/cluster/work/projects/<PROJECT>`, shared by everyone on
   the allocation. The run tree therefore lives at
   `/cluster/work/projects/<PROJECT>/$USER/<CONFIGNAME>/`.

2. **Work is purged after 21 days.** `SCRATCH/` can stay there, but move each experiment's
   `data/` output to `/cluster/projects/<PROJECT>` (or to NIRD) before it is deleted. The
   Python/container environment and compiled binaries must also live on
   `/cluster/projects/<PROJECT>`, not in the 20 GiB home (see
   [Python environment](installation.md#python-environment)).

:::{tip}
NIRD is not visible from the login nodes. To stage collaborator data from NIRD, copy it
on a service (SVC) node — where NIRD is read-write — or read it directly inside a batch
job, where compute nodes see it read-only. See
[Post-processing collaborator runs](mscprogs.md#post-processing-runs-from-nird) for the
recommended staging workflow.
:::

::::


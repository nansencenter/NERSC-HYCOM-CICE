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
| `WORK` | Root of the run/work filesystem | `/cluster/work/users` | `/cluster/projects/<PROJECT>` |
| `USER` | Your username on the system | `nlo043` | `nlo043` |
| `<PROJECT>` | Compute/storage allocation code | `nn9481k` | `nn9481k` |
| `<CONFIGNAME>` | Grid/configuration name | `TP2a0.10` | `TP2a0.10` |
| `<EXPT_ID>` | Experiment identifier (dot notation) | `04.2` | `04.2` |
| `<IEXPT>` | Integer experiment identifier (used in directory names) | `042` | `042` |

:::{note}
The Olivia `WORK` example is `/cluster/projects/<PROJECT>` (not the work filesystem) because
the recommended layout keeps the run tree on the non-purged project area and redirects only the
large `SCRATCH/` and `data/` directories to the purged scratch filesystem — and archives kept
output to NIRD. See [Olivia storage layout](#olivia-storage-layout) below.
:::

(<a name="olivia-storage-layout"></a>)
::::{dropdown} Olivia storage layout (NRIS/Sigma2)

Olivia's storage areas differ from Betzy's, so the `HOME`/`WORK` roles map onto different
filesystems and a few extra rules apply. See the
[Sigma2 storage documentation](https://documentation.sigma2.no/files_storage/clusters.html)
for the authoritative description.

`projects` and `work` are **different Lustre filesystems** with different sizes, striping,
and purge policies — the distinction drives the layout below:

| Area | Path | Capacity (shared by the allocation) | Striping | Backup | Purged | What it is for |
|---|---|---|---|---|---|---|
| Home | `/cluster/home/$USER` | ~20 GiB | — | yes | no | Dotfiles and small personal configs only — too small for the environment or build artifacts |
| Projects | `/cluster/projects/<PROJECT>` | 10 TiB hard cap | `stripe_count=1` (single OST) | no | **no** | **Small + permanent**: source, compiled binaries, container environment, experiment **configuration** and `build/` — *not* large output |
| Work | `/cluster/work/projects/<PROJECT>` | ~1.1 PB (hundreds of TiB free) | `stripe_count=4` | no | **yes, 21 days** | **Large + transient**: `SCRATCH/` and freshly-written `data/` output |
| NIRD | `/nird/datalake/<NS-PROJECT>`, `/nird/datapeak/<NS-PROJECT>` | large archive tier | — | — | no | **Large + keep long-term**: archived model output. Read-write on service (SVC) nodes, **read-only** on compute nodes, **absent** on login nodes |

**No per-user work area.** Unlike Betzy's `/cluster/work/users/$USER`, Olivia provides only
the per-project `/cluster/work/projects/<PROJECT>`, shared by everyone on the allocation. The
run tree therefore lives at `/cluster/work/projects/<PROJECT>/$USER/<CONFIGNAME>/`.

**Why the bulk run I/O (`SCRATCH/` and `data/`) goes on `work`, not `projects`.** It is
tempting to put the whole run on the non-purged `/cluster/projects` to dodge the 21-day clock,
but pointing the model's scratch and output I/O there causes real problems:

1. **Quota.** `projects` is a hard **10 TiB cap shared across the whole allocation**. A live
   run is the worst thing to point at it: `SCRATCH/` churns constantly and multi-year
   archive/restart output is large. Filling it doesn't just stop your job — it blocks
   **every collaborator on the allocation** from writing. `work` is a petascale scratch with
   hundreds of TiB free; the purge is the price for that capacity.
2. **Performance.** `projects` defaults to `stripe_count=1` (each file on a single object
   target), so large parallel writes serialise. `work` stripes 4-wide — which is what
   HYCOM's large `.a` archive and restart writes want. This is exactly why NRIS labels `work`
   "job data" and `projects` "user software".
3. **Shared-filesystem load.** Heavy job I/O on `projects` degrades it for everyone using it
   for storage.

:::{tip}
**When running in `projects` is fine:** a small or short one-off run whose footprint fits
comfortably in the free project quota — then you skip the purge bookkeeping. For production
or multi-year runs, use `work`.
:::

**So should everything sit in `/cluster/projects`? No — and neither should the large output.**
The 10 TiB project quota is shared by the whole allocation, and `data/` from a multi-year (and
especially BGC) run is far too big for it — which is exactly why shared model output lives on
NIRD. Split across **three tiers** by size and lifetime:

| Tier | What | Goes on |
|---|---|---|
| Small + permanent | source, container environment, experiment configuration, `build/` | `/cluster/projects/<PROJECT>` |
| Large + disposable | `SCRATCH/` | `/cluster/work/projects/<PROJECT>` |
| Large + keep | `data/` model output | written on `/cluster/work/projects`, then archived to **NIRD** |

**Recommended layout (hybrid).** Set `WORK=/cluster/projects/<PROJECT>` so the `<CONFIGNAME>/`
tree — configuration and `build/` — lives on the non-purged filesystem, then redirect the two
large directories to the work filesystem by overriding `S=` (scratch) and `D=` (data) in
`EXPT.src` (see [Configure EXPT.src](experiment-setup.md#configure-exptsrc)). Because `work` is
purged after 21 days — by file age, so output from early in a long run is at risk *while the run
is still going* — archive completed `data/` to NIRD on a **rolling** basis (e.g. per model year
as it finishes) from a service node, not just at the end. For a small run whose output
comfortably fits the project quota, leaving `data/` on `projects` is fine and skips the
archiving step.

```
/cluster/projects/<PROJECT>/$USER/          # small, persistent, NOT purged
├── NERSC-HYCOM-CICE/  fabm/  ersem/  nersc/    # source + compiled binaries (HOME role)
├── hycom-cice-env/                              # Python/container environment
└── <CONFIGNAME>/                                # experiment config (keep in git)
    ├── REGION.src
    └── expt_<EXPT_ID>/
        ├── EXPT.src, blkdat.input, srjob.sh     # configuration — must persist
        ├── build/.../hycom_cice                 # compiled executable
        ├── SCRATCH -> /cluster/work/projects/<PROJECT>/$USER/<CONFIGNAME>/expt_<EXPT_ID>/SCRATCH
        └── data    -> /cluster/work/projects/<PROJECT>/$USER/<CONFIGNAME>/expt_<EXPT_ID>/data

/cluster/work/projects/<PROJECT>/$USER/     # large, fast (4-wide stripe), PURGED after 21 days
└── <CONFIGNAME>/expt_<EXPT_ID>/
    ├── SCRATCH/                              # live run I/O — disposable
    └── data/                                 # model output — archive to NIRD before it ages out

/nird/datalake/<NS-PROJECT>/                 # large, long-term archive (push from a service node)
└── <run>/                                    # finished model output, kept
```

This redirect is safe because `expt_preprocess.sh` only creates (`mkdir -p`) and enters
(`cd`) the `$S` and `$D` directories — it never deletes them — so the overridden `S=`/`D=`
paths (or `SCRATCH`/`data` symlinks) are honoured transparently. The purge then only ever
touches `work`; the configuration and `build/` on `projects` are untouched, and the output is
yours to move to NIRD before it ages out.

:::{tip}
NIRD is not visible from the login nodes. To stage collaborator data from NIRD, copy it
on a service (SVC) node — where NIRD is read-write — or read it directly inside a batch
job, where compute nodes see it read-only. See
[Post-processing collaborator runs](mscprogs.md#post-processing-runs-from-nird) for the
recommended staging workflow.
:::

::::


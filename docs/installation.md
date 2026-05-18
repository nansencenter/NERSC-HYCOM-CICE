## Clone the repositories

Clone NERSC-HYCOM-CICE into `${HOME}` and switch to the `develop` branch:

```bash
cd ${HOME}
git clone https://github.com/nansencenter/NERSC-HYCOM-CICE.git
cd NERSC-HYCOM-CICE
git fetch --all
git checkout develop
```

Clone the biogeochemical components alongside it:

```bash
cd ${HOME}
git clone https://github.com/fabm-model/fabm.git
git clone https://github.com/pmlmodelling/ersem.git
git clone https://github.com/nansencenter/nersc.git
```

## Python environment

A dedicated conda environment keeps the model's Python dependencies isolated from other
projects and from the system Python, and makes the setup reproducible across machines.

### Set up conda

On HPC systems, conda may require machine-specific setup before creating the environment,
see the dropdown below for your machine. On a standard workstation, skip ahead to
[Create the environment](#create-the-environment).

::::{dropdown} Betzy (NRIS/Sigma2)

Betzy only provides Miniforge3 (not Anaconda or Miniconda). Full details are in the
[Sigma2 conda documentation](https://documentation.sigma2.no/software/userinstallsw/conda.html);
the essential steps are:

**1. Load and activate Miniforge3**

```bash
module load Miniforge3/24.1.2-0
source ${EBROOTMINIFORGE3}/bin/activate
```

After sourcing, your prompt should show `(base)`, confirming you are in the base conda environment. From here you can use `conda` to create and manage environments and install packages.

**2. Configure conda directories**

The default locations for package cache and environments are in `${HOME}`, which has
limited quota. Run the following once to redirect both to your project directory, replacing
`<PROJECT>` with your project code (e.g. `nn2993k`). Both settings are saved to `~/.condarc`.

```bash
conda config --append pkgs_dirs /cluster/projects/<PROJECT>/conda/${USER}/package-cache
conda config --append envs_dirs /cluster/projects/<PROJECT>/conda/${USER}
```

::::

### Create the environment

```bash
conda env create -f ${HOME}/NERSC-HYCOM-CICE/environment/python.yaml
conda activate hycom-cice
```

Then install the NERSC-specific libraries from the cloned repository:

```bash
pip install ${HOME}/NERSC-HYCOM-CICE/pythonlibs/modeltools
pip install ${HOME}/NERSC-HYCOM-CICE/pythonlibs/modelgrid
pip install ${HOME}/NERSC-HYCOM-CICE/pythonlibs/gridxsec
pip install ${HOME}/NERSC-HYCOM-CICE/pythonlibs/abfile
```

To upgrade the NERSC libraries, add `--upgrade` to each `pip install` command.

### Use the environment

The conda environment must be activated in any job or script that uses the model's Python
tools. On a standard workstation, `conda activate hycom-cice` is sufficient. On HPC
systems the conda installation itself must be loaded first.

::::{dropdown} Betzy (NRIS/Sigma2)

```{include} _snippets/betzy_python_activate.md
```

::::

## HPC environment

Each machine has a dedicated HPC environment file under `NERSC-HYCOM-CICE/environment/`.
Source it before compiling and include it in your job submission scripts before launching
the model.

```{include} _snippets/betzy_hpc_env.md
```

The dropdowns below render the actual environment files from `NERSC-HYCOM-CICE/environment/`
and are shown here for reference only. What you source is exactly what is displayed!

::::{dropdown} Betzy (NRIS/Sigma2)

:::{literalinclude} ../environment/betzy_env.sh
:language: bash
:::

::::

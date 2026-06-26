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

A dedicated Python environment keeps the model's Python dependencies isolated from other
projects and from the system Python, and makes the setup reproducible across machines.

### Set up conda

Betzy requires some one-time setup before conda can be used. On a standard workstation
or on Olivia, skip ahead to [Create the environment](#create-the-environment).

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

::::{dropdown} Betzy (NRIS/Sigma2) and workstation

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

::::

::::{dropdown} Olivia (NRIS/Sigma2)

Olivia uses [HPC-container-wrapper](https://documentation.sigma2.no/hpc_machines/olivia/software_stack.html#key-features-of-hpc-container-wrapper)
instead of conda directly. The wrapper encapsulates the environment inside a container,
reducing the number of visible files and improving file-system performance.

**1. Load the container wrapper**

Run the following in your current session. These are only needed during the build, so there is no need
to add them to `~/.bashrc`.

```bash
export http_proxy=http://10.63.2.48:3128/
export https_proxy=http://10.63.2.48:3128/
module load NRIS/CPU
module load hpc-container-wrapper
```

**2. Choose an installation path**

Set a variable pointing to where the environment will live permanently. This directory
must be kept after installation, as it contains the container and the executables:

```bash
export INSTALL_DIR=/cluster/projects/nn2993k/${USER}/hycom-cice-env
```

**3. Build the environment**

```bash
conda-containerize new --mamba \
    --prefix /cluster/projects/nn2993k/${USER}/hycom-cice-env \
    ${HOME}/NERSC-HYCOM-CICE/environment/python.yaml
```

**4. Install the NERSC-specific libraries**

Since the containerised environment cannot be modified directly, install the local
libraries via a post-installation script:

```bash
cat > /tmp/nersc_libs.sh << 'EOF'
pip install ${HOME}/NERSC-HYCOM-CICE/pythonlibs/modeltools
pip install ${HOME}/NERSC-HYCOM-CICE/pythonlibs/modelgrid
pip install ${HOME}/NERSC-HYCOM-CICE/pythonlibs/gridxsec
pip install ${HOME}/NERSC-HYCOM-CICE/pythonlibs/abfile
EOF
conda-containerize update /cluster/projects/nn2993k/${USER}/hycom-cice-env --post-install /tmp/nersc_libs.sh
```

**5. Activate the environment**

Run the following to make the environment available in your current session:

```bash
export PATH="/cluster/projects/nn2993k/${USER}/hycom-cice-env/bin:${PATH}"
```

Run this command at the start of each session or job script where you need the
environment. If you only use one Python environment, you can also add it to `~/.bashrc`
to activate it automatically.

:::{tip}
The steps above are also scripted under `environment/py-nhc/`: set `$PYLIBSPATH` in
`py-nhc-env-post.sh`, then run `create_python_env.sh` (with the `hpc-container-wrapper`
module loaded). It builds the container into `environment/py-nhc/py-nhc-env/`; add that
`bin/` directory to `$PATH` to activate.
:::

::::

### Use the environment

The conda environment must be activated in any job or script that uses the model's Python
tools. On a standard workstation, `conda activate hycom-cice` is sufficient. On HPC
systems the conda installation itself must be loaded first. On Olivia, the containerised
environment is activated by prepending its `bin` directory to `PATH`.

::::{dropdown} Betzy (NRIS/Sigma2)

```{include} _snippets/betzy_python_activate.md
```

::::

::::{dropdown} Olivia (NRIS/Sigma2)

```{include} _snippets/olivia_python_activate.md
```

::::

### Update the environment

::::{dropdown} Betzy (NRIS/Sigma2) and workstation

To sync the environment with changes to `python.yaml`:

```bash
conda activate hycom-cice
conda env update -f ${HOME}/NERSC-HYCOM-CICE/environment/python.yaml --prune
```

To add a package not in `python.yaml`:

```bash
conda activate hycom-cice
conda install <package>        # or: pip install <package>
```

::::

::::{dropdown} Olivia (NRIS/Sigma2)

The containerised environment cannot be modified directly. Use `conda-containerize update`
with a post-installation script listing the changes. First set the path to your existing
environment and reload the modules:

```bash
export http_proxy=http://10.63.2.48:3128/
export https_proxy=http://10.63.2.48:3128/
module load NRIS/CPU
module load hpc-container-wrapper
```

Then create a script with the packages to add or remove, and apply it:

```bash
cat > /tmp/update.sh << 'EOF'
conda install -y <package>
pip install <package>
EOF
conda-containerize update /cluster/projects/nn2993k/${USER}/hycom-cice-env --post-install /tmp/update.sh
```

If `python.yaml` has changed substantially (e.g. a new core dependency was added),
it is cleaner to rebuild from scratch by deleting the existing environment and repeating
the [Create the environment](#create-the-environment) steps:

```bash
rm -rf /cluster/projects/nn2993k/${USER}/hycom-cice-env
```

::::

## HPC environment

Each machine has a dedicated HPC environment file under `NERSC-HYCOM-CICE/environment/`.
Source it before compiling and include it in your job submission scripts before launching
the model.

::::{dropdown} Source HPC environment — Betzy (NRIS/Sigma2)

```{include} _snippets/betzy_hpc_env.md
```

::::

::::{dropdown} Source HPC environment — Olivia (NRIS/Sigma2)

```{include} _snippets/olivia_hpc_env.md
```

::::

The dropdowns below show the full content of each `*_env.sh` for reference — what you
source is exactly what is displayed.

::::{dropdown} View betzy_env.sh

:::{literalinclude} ../environment/betzy_env.sh
:language: bash
:::

::::

::::{dropdown} View olivia_env.sh

:::{literalinclude} ../environment/olivia_env.sh
:language: bash
:::

::::

## Environment setup

Each machine has a dedicated environment file under
`NERSC-HYCOM-CICE/environment/`. Source the appropriate file before compiling and include it in
your job submission scripts before launching the model.

```bash
source ${HOME}/NERSC-HYCOM-CICE/environment/betzy_env.sh   # adjust filename for your machine
```

The environment files load the required HPC modules, and configure the stack size.

## Modules and settings by machine

The dropdowns below render the actual environment files from `NERSC-HYCOM-CICE/environment/`
and are shown here for reference only. What you source is exactly what is displayed!

::::{dropdown} Betzy (NRIS/Sigma2)

:::{literalinclude} ../environment/betzy_env.sh
:language: bash
:::

::::


# TP2_setup
Instruction manual for installing TP2 (hycom 17km Arctic Ocean configuration) with ECOSMO (BGC module) on sigma2 HPC (Betzy etc.). This manual is developed based on 
[Shuang's TP2 setup note](https://docs.google.com/document/d/1z53p0V0knFjaSZTXhO1WUQvG1OMDk7nzLdp_zZKGwJk/edit?tab=t.0) and [REDAME.betzy](https://github.com/nansencenter/NERSC-HYCOM-CICE/blob/develop/README.betzy). 
Final goal of this repository is to provide bash installer and job submitter to complete neccessary tasks at command line (not there yet).

## Table of Contents
- [Rquirements](#requirements)
  - [*Model components*](#model-components)
  - [*Bash environment*](#bash-environment)
  - [*HPC modules on Betzy*](#hpc-modules-on-betzy)
- [Installation](#installation)
  - [*Environment variables for installation settings*](#environment-variables-for-installation-settings)
  - [*Cloning the model components*](#cloning-the-model-components)
- [Createing a new experiment](#createing-a-new-experiment)
  - [*Copy configuration to a work directory*](#copy-configuration-to-a-work-directory)
  - [*Modify experiment settings in REGION.src*](#modify-experiment-settings-in-region.src)
  - [*Create a new experiment folder*](#create-a-new-experiment-folder)
  - [*Copy and edit hycom_opt*](#copy-and-edit-hycom_opt)
  - [*Modify experiment settings in blkdat.input*](#modify-experiment-settings-in-blkdat.input)
  - [*Modify experiment settings in EXPT.src*](#modify-experiment-settings-in-expt.src)
- [Compiling hycom executable (hycom_cice)](#compiling-hycom-executable-(hycom_cice))
  - [*Compile libhycnersc.a/MSCPROGS*](#compile-libhycnersc.a/mscprogs)
  - [*Compile libfabm.a/FABM*](#compile-libfabm.a/fabm)
  - [*Compile HYCOM-CICE*](#compile-hycom-cice)
  - [*Compile hycom_ALL*](#compile-hycom_all)
- [Preparing external files](#preparing-external-files)
  - [*Prepare relaxation (climatology), river input and MPI tile definition files*](#prepare-relaxation-(climatology),-river-input-and-mpi-tile-definition-files)
  - [*Prepare nesting time scale definition files*](#prepare-nesting-time-scale-definition-files)
  - [*Prepare sponge layer settiing files*](#prepare-sponge-layer-settiing-filesc*)
- [Creating a job settings](#creating-a-job-settings)
  - [*Prepare restart files*](#prepare-restart-files)
  - [*Prepare off-line nesting files*](#prepare-off-line-nesting-files)
  - [*Prepare atmospheric forcing files*](#prepare-atmospheric-forcing-files)
  - [*Prepare SCRATCH folder*](#prepare-scratch-folder)
- [Job submission](#job-submission)
- [Visualization](#visualization)

### Rquirements

Installation of TOPAZ ocean model system with biogeochemical modules requires the following model components and environment setups. 

#### *Model components*

Each component is downloaded from its associated git repository:

1. [NERSC-HYCOM-CICE (develop branch)](https://github.com/nansencenter/NERSC-HYCOM-CICE/tree/develop)
   NERSC-HYCOM-CICE is an ocean-sea ice-biogeochemistry coupled model framework developed at NERSC.
2. [ECOSMO_operational (main branch)](https://github.com/nansencenter/nersc/blob/main/ecosmo/ecosmo_operational.F90)
   ECOSMO_operational is a biogeochemical modeling module optimized for operational tasks at NERSC. EOCSMO_operational package (nersc) is downloaded as a part of FABM package, but we use the latest version from NERSC git repository.
3. [GOTM](https://github.com/gotm-model/code)
   GOTM is a 1D column ocean model. In this installation we use a light module of GOTM installed as a part of FABM package.
5. [FABM](https://github.com/fabm-model/fabm)
   FABM is a Fortran 2003 programming framework (a coupler) for biogeochemical models of marine and freshwater systems. 
6. [ERSEM](https://github.com/pmlmodelling/ersem)
   ERSEM is a marine biogeochemical and ecosystem model. It describes the cycling of carbon, nitrogen, phosphorus, silicon, oxygen and iron through the lower trophic level pelagic and benthic ecosystems. In this installation we use a carbon chenistry module of ERSEM.
7. NERSC Python library (downloaded together with NERSC-HYCOM-CICE package).

#### *Bash environment*

Add the following line to `.bashrc`:

```bash
export LANG=en_US.UTF-8
export LC_ALL=en_US.utf8
```

#### *HPC modules on Betzy*

Add the following lines to `.bashrc`: 

```bash
ml purge
ml load CMake/3.23.1-GCCcore-11.3.0
ml load ESMF/8.3.0-iomkl-2022a
ml load FFTW/3.3.10-GCC-11.3.0
ml load UDUNITS/2.2.28-GCCcore-11.3.0
ml load Python/3.10.4-GCCcore-11.3.0
ml load GSL/2.7-intel-compilers-2022.1.0
```

After editing `.bashrc`, re-login or perform `source ~/.bashrc`. See [REDAME.betzy](https://github.com/nansencenter/NERSC-HYCOM-CICE/blob/develop/README.betzy) for  extra information.

### Installation

In this instruction, we use the following directory tree structure for hycom (version 2.2) installation on Betzy:
```bash
  /cluster
  ├── home
  │   └── ${USERNAME}                            # HOME
  │       ├── ${MODELNAME}                       # HOME_HYCOM (HYCOM-CICE home)
  │       │   └── ${HYCOM_REPO}
  │       │       ├── pythonlibs                 # LIB_PYTHON (NERSC HYCOM Python libs)
  │       │       ├── bin                        # BIN_HYCOM (NERSC HYCOM utility scripts)
  │       │       └── hycom
  │       │           ├── RELO
  │       │           │   └── src_2.2.98ZA-07Tsig0-i-sm-sse_relo_mpi # SRC_HYCOM (HYCOM source files)
  │       │           └── MSCPROGS
  │       │               ├── bin                # BIN_HYCOM_MSCPROGS (MSCPROGS bin folder)
  │       │               └── lib/libhycnersc.a  # LIB_HYCNERSC (MSCPROGS Nersclib library)
  │       ├── FABM                               # HOME_FABM (FABM/ECOSMO home)
  │       └── local/fabm/hycom/lib64/libfabm.a   # LIB_FABM (FABM library)    
  │
  └── work
      └── users
          └── ${USERNAME}
              └── ${MODELNAME}                   # WORK_HYCOM (HYCOM-CICE work)
                  └── ${CONFIGNAME}              # WORK_HYCOM_CNF (configuration folder)
                      ├── bin                    # 
                      ├── REGION.src
                      ├── expt_${NEWEXPERIMENT}  # WORK_HYCOM_EXP (experiment folder)
                      │   ├── EXPT.src
                      │   ├── SCRATCH            # WORK_HYCOM_SCR (SCRATCH folder)
                      │   ├── data               # WORK_HYCOM_DAT (data/restart files folder)
                      │   └── build
                      │       └── src_2.2.98ZA-07Tsig0-i-sm-sse_relo_mpi # HYCOM source files (copy from SRC_HYCOM)
                      │           └── hycom_cice # hycom executable
                      ├── nest
                      │   └── ${IEXPT}           # WORK_HYCOM_NST (nesting files folder)
                      ├── relax
                      │   └── ${IEXPT}           # WORK_HYCOM_RLX (relaxation files folder)
                      └── force
                          └── synoptic
                              └── ${IEXPT}       # WORK_HYCOM_FRC (forcing files folder)
```

#### *Environment variables for installation settings*

In order to make installation processes available through bash scripts, we set some key environment variables first:
```bash
export CONFIGNAME="TP2a0.10"   # configuration name 
export NEWEXPERIMENT=04.2      # new experiment number
export T="04"                  # topography version number
export IVERSN=22               # hycom version number x10
export MXBLCKS=9               # maximal number of ice blocks relates to MPI cores distribution
export ICECLIM=0               # prepare initial ice file from climatology (1) or not (otherwise)
export INITFLG=""              # start physics from climatology "--init" or from restart ""
export RSTRFQ=7                # restart file dumping frequency [day]
export MSLPRF=0                # msl pressure forcing flag (0=F,1=T)
export NRDFLG=0                # options for the radiation heat fluxes(0,1,...,7). Only to hycom 2.2 
export COMPILE_BIOMODEL="yes"  # turn on/off ("yes/no") BGC (FABM-ECOSMO) module
export NTRACR=1                # BGC on/off: 0-physics only; 1-biology restart; -1-biology initialized with climatology
export NMPI=504                # number of MPI tiles
```
for TP2 configuration with ECOSMO (BGC module) on Betzy. To setup hycom without ECOSMO choose:
```bash
export COMPILE_BIOMODEL="no"   # turn on/off ("yes/no") BGC (FABM-ECOSMO) module
export NTRACR=0                # BGC on/off: 0-physics only; 1-biology restart; -1-biology initialized with climatology
```
Note `NEWEXPERIMENT` is user defined variable with format `NN.M` that defines experiment folder name (see directory tree above), where `NN` is experiment folder ID (two degits integer) and `M` is its version ID  (sigle degit integer). 

Define HYCOM-CICE repository as
```bash
export HYCOM_REPO=NERSC-HYCOM-CICE                               # HYCOM-CICE repository name 
```
for hindcast/forecast run or as
```bash
export HYCOM_REPO=NERSC-HYCOM-CICE_BIORANv2                      # HYCOM-CICE repository name 
```
for BGC reanaysis run.

Next, we define `MODELNAME`. Choice of `MODELNAME` is up to user and for example, you can use
```bash
export MODELNAME='TP2_HINDCAST'                                              # common model folder name under HOME/WORK
```
or 
```bash
export MODELNAME='TP2_REANALYSIS'                                            # common model folder name under HOME/WORK
```

Another example is to use shortname of model configuration TP2a0.10:
```bash
set -u
export SHRTNAME="$(echo "$CONFIGNAME" | sed -E 's/[a-z].*$//')"              # eg. TP2a0.10 > TP2
export MODELNAME=$SHRTNAME                                                   # common model folder name under HOME/WORK
```

After setting bash variables, skelton of the directory tree can be initialized by bash scripts:

```bash
# setup directory tree
set -u

## define variables
export USERNAME=$(whoami)                                                      # your username
export IEXPT=$(echo "$NEWEXPERIMENT" | sed 's/\.//')                           # eg. "03.0" > "030" 
export HOME_HYCOM=/cluster/home/${USERNAME}/${MODELNAME}                       # HYCOM-CICE home directory
export HOME_FABM=/cluster/home/${USERNAME}/FABM                                # FABM(ECOSMO) home directory
export WORK_HYCOM=/cluster/work/users/${USERNAME}/${MODELNAME}                 # HYCOM-CICE work directory
export LIB_PYTHON=$HOME_HYCOM/$HYCOM_REPO/pythonlibs                           # NERSC hycom python libraries directory
export LIB_HYCNERSC=$HOME_HYCOM/$HYCOM_REPO/hycom/MSCPROGS/lib/libhycnersc.a   # Hycomlib library
export LIB_FABM=/cluster/home/${USERNAME}/local/fabm/hycom/lib64/libfabm.a     # FABM library
export BIN_HYCOM=$HOME_HYCOM/$HYCOM_REPO/bin                                   # hycom utility bin diectory
export BIN_HYCOM_MSCPROGS=$HOME_HYCOM/$HYCOM_REPO/hycom/MSCPROGS/bin           # hycom MSCPROGS bin diectory
export PATH="$BIN_HYCOM_MSCPROGS:$BIN_HYCOM:$HOME_HYCOM/$HYCOM_REPO/bin:$PATH" # add utilty bins to PATH

## create HOME_ and WORK_ directories
mkdir -p $HOME_HYCOM
mkdir -p $HOME_FABM        
mkdir -p $WORK_HYCOM
```
where `HOME_##` specifies a folder for storing model source files (*home*) and `WORK_##` specifies a folder for model configuration and execution (*build*).

#### *Cloning the model components*

Now you are ready to clone hycom under `HOME_HYCOM`:
```bash
set -u
# clone HYCOM-CICE
cd $HOME_HYCOM
git clone https://github.com/nansencenter/${HYCOM_REPO}.git

# switch to develop branch
cd $HOME_HYCOM/${HYCOM_REPO}
git fetch --all
git checkout develop
```
and to clone biogeochemical modules under `HOME_FABM`:
```bash
set -u
cd $HOME_FABM
# clone FABM
git clone https://github.com/fabm-model/fabm.git         # FABM

# clone ERSEM
git clone https://github.com/pmlmodelling/ersem.git      # ERSEM

# clone ECOSMO_operational
git clone git@github.com:nansencenter/nersc.git          # ECOSMO
#cp -r /cluster/projects/nn9481k/caglar/ECOSMO_code/nersc # ECOSMO
```

Nersc Python libraries can be installed by
```bash
pip install --user scipy
pip install --user cfunits
pip install --user netCDF4
pip install --user netcdftime
pip install --user numpy
pip install --user f90nml
pip install --user pyproj
```
and
```bash
pip install --user $LIB_PYTHON/modeltools
pip install --user $LIB_PYTHON/modelgrid
pip install --user $LIB_PYTHON/gridxsec
pip install --user $LIB_PYTHON/abfile
```
for the first installation and 
```bash
pip install --user --upgrade $LIB_PYTHON/modeltools
pip install --user --upgrade $LIB_PYTHON/modelgrid
pip install --user --upgrade $LIB_PYTHON/gridxsec
pip install --user --upgrade $LIB_PYTHON/abfile
```
for updating the exsisting libraries.

### Createing a new experiment

Before compilation of hycom executable, we need to create a new experiment folder ``.

#### *Copy configuration to a work directory*

First, copy hycom configuration and create symlink of hycom bin. For TP2, perform:
```bash
set -u
cd $WORK_HYCOM
cp -r $HOME_HYCOM/${HYCOM_REPO}/$CONFIGNAME .
cd $WORK_HYCOM/$CONFIGNAME
[ -f bin ] && rm bin
ln -sf $HOME_HYCOM/${HYCOM_REPO}/bin .
```

#### *Modify experiment settings in REGION.src*

Next, you need to prepare `REGION.src` for your experimental settings BEFORE making a new *experiment*. 
Normally, this step is done manually, but here is a bash script to achieve the step at command line:
```bash
set -u
cd $WORK_HYCOM/$CONFIGNAME
cp $HOME_HYCOM/${HYCOM_REPO}/input/REGION.src .
sed -i "/^export R=/c export R=$CONFIGNAME" REGION.src
sed -i "/^export NHCROOT=/c export NHCROOT=$HOME_HYCOM/${HYCOM_REPO}" REGION.src
```
where `R` is model configuration name and `NHCROOT` is where HYCOM-CICE source files are downloaded.

#### *Create a new experiment folder*

Now you create a new *experiment* from the existing *experiment* (01.0) under hycom work/configuration directory (which is the case if you clone HYCOM-CICE from git repository).
```bash
set -u
cd $WORK_HYCOM/$CONFIGNAME
bin/expt_new.sh 01.0 $NEWEXPERIMENT
```

After a new *experiment* folder `expt_$NEWEXPERIMENT` is created, some adjustments of `blkdat.input` and `EXPT.src` under the folder are required.

#### *Copy and edit hycom_opt* ⚠️**NEW**

```bash
set -u
cd $WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT
cp $HOME_HYCOM/$HYCOM_REPO/TP5a0.06/expt_01.0/hycom_opt .
sed -i 's/sssrmx_scalar=99\./sssrmx_scalar=.5/' hycom_opt
```

#### *Modify experiment settings in blkdat.input*

Experiment settings defined in the original `blkdat.input` are needed to be adjusted for specific hycom settings. 
You can copy `blkdat.input` from existing TP2 configuration with nesting boundary
```bash
BLKDATIN=/nird/datalake/NS9481K/shuang/TP2_setup/exp02.6_seaclim_ref/blkdat.input # TP2 with nesting boundary
```
or with climatology boundary
```bash
BLKDATIN=/nird/datalake/NS9481K/shuang/TP2_setup/exp02.6_seaclim_ref/blkdat.input_clim # TP2 with climatology boundary
```
then run
```bash
set -u
cd $WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT
mv blkdat.input blkdat.input_backup
cp $BLKDATIN blkdat.input
```

After importing new `blkdat.input`, you need to modify `iexpt`, `ntracr` and `rstefq` in `blkdat.input` depending on your experiment settings manually with file editor or run scripts below at command line:
```bash
set -u
vars=("$IEXPT" "$RSTRFQ" "$NTRACR")
keys=( "iexpt " "rstrfq"  "ntracr")
```
first and run:
```bash
set -u

cd $WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT
filename='blkdat.input'

# Define a dictionary (an associative array) of keyword and replacement.
declare -A replacements 
for i in "${!keys[@]}"; do
    replacements["${keys[$i]}"]="${vars[$i]}"
done

tempfile=$(mktemp)
while IFS= read -r line; do
    for keyword in "${!replacements[@]}"; do
        replacement="${replacements[$keyword]}"
        if [[ "$line" == *"$keyword"* ]]; then
            line=$(echo "$line" | sed -E "s/^([[:space:]]*)(-?[0-9]+)/\1$replacement/")
            break
        fi
    done
    echo "$line" >> "$tempfile"
done < "$filename"
mv "$tempfile" "$filename"

for key in "${keys[@]}"; do
    grep "$key" $filename
done
```

With hycom 2.2, we create new `blkdat.input` based on `blkdat.input_H2.298` by first setting environment variables:
```bash
set -u

cd $WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT

export VELDF4=-2               # diffusion velocity (m/s) for biharmonic momentum dissip. (use < 0 for reading exterbal file)
export THKDF4=-2               # diffusion velocity (m/s) for biharmonic thickness diffus. (use < 0 for reading exterbal file)
export TICEGR=2                # ENLN: temp. grad. inside ice (deg/m); =0 use surtmp
export MSLPRF=0                # msl pressure forcing flag (0=F,1=T)

vars=("$IEXPT" "$RSTRFQ" "$NTRACR" "$VELDF4" "$THKDF4" "$TICEGR" "$MSLPRF")
keys=( "iexpt " "rstrfq"  "ntracr"  "veldf4"  "thkdf4"  "ticegr"  "mslprf")

cp blkdat.input_H2.298 blkdat.input
```
then run:
```bash
set -u

cd $WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT
filename='blkdat.input'

# Define a dictionary (an associative array) of keyword and replacement.
declare -A replacements 
for i in "${!keys[@]}"; do
    replacements["${keys[$i]}"]="${vars[$i]}"
done

tempfile=$(mktemp)
while IFS= read -r line; do
    for keyword in "${!replacements[@]}"; do
        replacement="${replacements[$keyword]}"
        if [[ "$line" == *"$keyword"* ]]; then
            line=$(echo "$line" | sed -E "s/^([[:space:]]*)(-?[0-9]+)/\1$replacement/")
            break
        fi
    done
    echo "$line" >> "$tempfile"
done < "$filename"
mv "$tempfile" "$filename"

for key in "${keys[@]}"; do
    grep "$key" $filename
done
```
⚠️**Note**: This script is confirmed only with keywords in the above sample script. When you add other keyword to `replacements`, make sure first the string replacement is working as intended. 

#### *Modify experiment settings in EXPT.src*

Under `$WORK_HYCOM/$CONFIGNAME/topo` folder, you can find model topography data files `depth_TP2a0.10_01.[a,b]` where the last two digits represnt version number of the topography files. Here, we choose version `04` by setting `T="04"` in `EXPT.src`. Newly added variable `COMPILE_BIOMODEL` defines if HYCOM-CICE is coupled with ECOSMO. When `COMPILE_BIOMODEL="yes"`, the coupler is turned on and is turned off otherwise (eg, `COMPILE_BIOMODEL="no"`). The maximal number of ice blocks `MXBLCKS` relates to the MPI cores distribution, and could be recommened by the error message in a model run.
Normally, this step is done manually, but here is a bash script to achieve the step at command line:
```bash
set -u
cd $WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT
#echo ""                                                >> EXPT.src
#echo "# add new parameters for FABM and ICE"           >> EXPT.src
#echo "export MXBLCKS=${MXBLCKS}"                       >> EXPT.src # maximal number of ice blocks 
#echo "export COMPILE_BIOMODEL=\"${COMPILE_BIOMODEL}\"" >> EXPT.src # FABM coupler ON ("yes") or OFF ("no")
sed -i "/^X=/c X=\"$NEWEXPERIMENT\"" EXPT.src                       # X based on dir name
sed -i "/^E=/c E=\"$IEXPT\"" EXPT.src                               # E is X without "."
sed -i "/^T=/c T=\"$T\"" EXPT.src                                   # topography version
sed -i "/^export NMPI=/c export NMPI=$NMPI" EXPT.src
sed -i "/^export MXBLCKS=/c export MXBLCKS=$MXBLCKS" EXPT.src
sed -i "/^export COMPILE_BIOMODEL=/c export COMPILE_BIOMODEL=\"${COMPILE_BIOMODEL}\"" EXPT.src
```

- this command should be part of expt_new.sh with COMPILE_BIOMODEL, T as command line arguments.
- if you get error message about ice blocks exceed maximum, then change this number to the once recommend in the error message.

### Compiling hycom executable (hycom_cice)

Compilation of HYCOM_CICE with biogeochemical modules requires 2 pre-compiled libraries `libhycnersc.a` and `libfabm.a`.

#### *Compile libhycnersc.a/MSCPROGS*

First, setup proper Makefile include file `make.inc` for your environment:

```bash
set -u
cd $HOME_HYCOM/${HYCOM_REPO}/hycom/MSCPROGS/src/Make.Inc

compiler=ifort # FORTRAN compiler (eg. gfortran, gnu)
hostname=$(head -n 1 /etc/motd | awk -F' ' '{split($3, parts, "."); print parts[1]}') # extract hostname from motd
echo $hostname
[ -L make.inc ] && rm make.inc
ln -sf make.$hostname.$compiler make.inc
ls -l make.inc
```
⚠️**Note**: This script is confirmed only on sigma2 Fram and Betzy.

For example on Betzy, this script creates symlink `make.inc`:

```bash
make.inc -> make.betzy.ifort
```

If you set `HYCOM_REPO=NERSC-HYCOM-CICE_BIORANv2`, add the following patch before compilation (TODO: fix them on git repository):
```bash
set -u
mkdir -p $HOME_HYCOM/${HYCOM_REPO}/hycom/MSCPROGS/src/Average/TMP     # debug!
mkdir -p $HOME_HYCOM/${HYCOM_REPO}/hycom/MSCPROGS/src/Hycom_mean/TMP  # debug!
mkdir -p $HOME_HYCOM/${HYCOM_REPO}/hycom/MSCPROGS/src/Perturb_Forcing-2.2.98_TP5/TMP # debug!
mkdir -p $HOME_HYCOM/${HYCOM_REPO}/hycom/MSCPROGS/src/Perturb_Parameter/TMP # debug!
```

Now, compile MSCPROGS:
```bash
set -u
cd $HOME_HYCOM/${HYCOM_REPO}/hycom/MSCPROGS/src
gmake clean
gmake all
```

Once compilation successed, install executables under `MSCPROGS/bin` and `MSCPROGS/bin_setup`:

```bash
gmake install
```

After the compilation/installation, `libhycnersc.a` can be found under
```bash
$HOME_HYCOM/$HYCOM_REPO/hycom/MSCPROGS/lib
```

⚠️**Note**: The installation should be repeated everytime when HPC modules are renewed.

For further detail, see [Procedure for compiling/installing MSCPROGS](https://github.com/nansencenter/NERSC-HYCOM-CICE/tree/develop/hycom/MSCPROGS)

#### *Compile libfabm.a/FABM*

For clean install, remove `build` folder:
```bash
set -u
cd $HOME_FABM
[ -d build ] && rm -rf build
```
then compile fabm:
```bash
set -u
cd $HOME_FABM
mkdir build; cd build
cmake $HOME_FABM/fabm -DFABM_HOST=hycom -DCMAKE_Fortran_COMPILER=ifort -DFABM_INSTITUTES="ersem;nersc;gotm" -DFABM_NERSC_BASE=$HOME_FABM/nersc -DFABM_ERSEM_BASE=$HOME_FABM/ersem
make install
```

After the compilation, `libfabm.a` can be found under
```bash
${HOME}/local/fabm/hycom/lib64
```

#### *Compile HYCOM-CICE*

⚠️**Note**: BGC module should be turned on in `blkdat.input` for physics-BGC coupled model compilation by setting bash environment variable `NTRACR`.

Before the compilation, we need to copy FABM configuration files, `fabm.yaml` and `hycom_fabm.nml` and CICE namelist `ice_in`, to the `experiment` folder. In this installation manual, we use experiment setups for TP2 hindcast run:

```bash
set -u
cd $WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT
cp  /nird/datalake/NS9481K/shuang/TP2_setup/exp02.6_seaclim_ref/fabm.yaml .
cp  /nird/datalake/NS9481K/shuang/TP2_setup/exp02.6_seaclim_ref/hycom_fabm.nml .
cp  /nird/datalake/NS9481K/shuang/TP2_setup/exp02.6_seaclim_ref/ice_in .
```

For clean install, remove `build` folder:
```bash
set -u
cd $WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT
[ -d build ] && rm -rf build
```
then compile HYCOM-CICE:
```bash
set -u
cd $WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT
ln -sf $HOME_HYCOM/${HYCOM_REPO}/bin/compile_model.sh .    # temporal solution
ln -sf $HOME_HYCOM/${HYCOM_REPO}/bin/common_functions.sh . # temporal solution
bash compile_model.sh ifort -u
```

Once compilation finished, HYCOM-CICE executable `hycom_cice` is created under 
```bash
$WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT/build/src_2.2.98ZA-07Tsig0-i-sm-sse_relo_mpi`
```

For extra details, see [README.betzy](https://github.com/nansencenter/NERSC-HYCOM-CICE/blob/develop/README.betzy) or [README.fram](https://github.com/nansencenter/NERSC-HYCOM-CICE/blob/develop/README.fram).

#### *Compile hycom_ALL*

Before preparing external files to run hycom, hycom utilities `hycom_ALL` should be compiled. First, setup proper Makefile include file `Make_all.src` for your environment. Normally, this step is done manually, but here is a bash script to achieve the step at command line:
```bash
set -u
cd $HOME_HYCOM/${HYCOM_REPO}/hycom/hycom_ALL/hycom_2.2.72_ALL
compiler=intelIFC
if [ -f config/${compiler}_setup ]; then
   sed -i "/^setenv ARCH /c setenv ARCH $compiler" Make_all.src
else
   echo "Can not find setup file ${compiler}_setup, EXIT"
fi
```
⚠️**Note**: This script is confirmed only on sigma2 Fram and Betzy.

Then compile:
```bash
csh ./Make_clean.com
csh ./Make_all.com
csh ./Make_ncdf.com
```
Note compilation of `Make_ncdf.com` would be stuck due to error for the moment.
For further detail, see [Procedure for compiling programs in hycom_2.2.72_ALL](https://github.com/nansencenter/NERSC-HYCOM-CICE/blob/develop/hycom/hycom_ALL/README.md)

⚠️**Note [2026/01/14]**: These files failed to be compiled in ```Make_all.com```
```bash
hycom_profile_hybgen.f
hycom_profile_hybgen+.f
```

### Preparing external files

#### *Prepare relaxation (climatology), river input and MPI tile definition files*

First edit `create_ref_case.sh` under `$HOME_HYCOM/${HYCOM_REPO}/bin/`: 
```bash
set -u
cd $WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT
cp $HOME_HYCOM/${HYCOM_REPO}/bin/create_ref_case.sh .
ICORE=29    # number of tiles in i directrion
JCORE=26    # number of tiles in j directrion
sed -i "/^Icore=/c Icore=$ICORE" create_ref_case.sh
sed -i "/^Jcore=/c Jcore=$JCORE" create_ref_case.sh
sed -i "/^iceclim=/c iceclim=$ICECLIM" create_ref_case.sh
sed -i "/^iceclim=/c iceclim=$ICECLIM" create_ref_case.sh
sed -i 's/tile_grid\.sh/tile_grid.sh -s 1/g' create_ref_case.sh
chmod +x create_ref_case.sh
```
where `ICORE` and `JCORE` are for MPI tile defnition files which are created by `tile_grid.sh` which is called inside of `create_ref_case.sh`. ```-s``` option to `tile_grid.sh` adds extra freedom in the tiling operation (see Shuang's note for further detail).

⚠️**New** In order to optimize tiling (see Shuang's note below), one line in ```create_ref_case.sh``` should be edited from

```bash
tile_grid.sh -${Icore} -${Jcore} ${T} > $EDIR/log/ref_tiling.out 2>&1
```
to
```bash
tile_grid.sh -s 1 -${Icore} -${Jcore} ${T} > $EDIR/log/ref_tiling.out 2>&1
```
where ```-s``` option adds extra freedom in the tiling operation.

Next, execute:
```bash
set -u
cd $WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT
bash create_ref_case.sh
```
which generates the following files:

- phy relax z climatology
- phy relax hybrid climatology
- bio relax z climatology      (when ntracr=1/-1 in blkdat.in)
- bio relax hybrid climatology (when ntracr=1/-1 in blkdat.in)
- co2 relax z climatology      (when ntracr=1/-1 in blkdat.in)
- co2 relax hybrid climatology (when ntracr=1/-1 in blkdat.in)
- relaxation mask
- sea ice cover climatology    (when ICECLIM=0)

under `$WORK_HYCOM/$CONFIGNAME/relax/expt_$NEWEXPERIMENT` folder and

- river forcing

under `$WORK_HYCOM/$CONFIGNAME/force/rivers/$IEXPT` folder and

- kpar forcing

under `$WORK_HYCOM/$CONFIGNAME/force/seawifs` folder and

- MPI tile definition files:

```bash
depth_TP2a0.10_04.0504
depth_TP2a0.10_04.0504.ppm
```

under `$WORK_HYCOM/$CONFIGNAME/topo/partit` folder. Four degits number `0504` in the MPI tile definitio file name indicates parameter `NMPI` should be set to `504` on `EXPT.src`.

#### *Prepare nesting time scale definition files*

When nesting opetion is on ( in `blkdat,input`), nesting time scale definition files should be available under `$WORK_HYCOM/$CONFIGNAME/nest/$IEXPT`. Here we simply copy them from existing experiment settings:
```bash
DIR_NST=/nird/datalake/NS9481K/shuang/nest/TP2_expt023
```
for TP2 (where is location for TP5?) and run
```bash
set -u
cd $WORK_HYCOM/$CONFIGNAME/nest/$IEXPT
cp ${DIR_NST}/ports.input .
cp ${DIR_NST}/rmu.a .
cp ${DIR_NST}/rmu.b .
cp ${DIR_NST}/rmutr.a .
cp ${DIR_NST}/rmutr.b .
```

#### *Prepare sponge layer settiing files*

If sponge layer is turned on (specify negative values to ```thkdf4``` and ```veldf4``` in ```blkdat.input```), the following layer diffusion files should be available under `$WORK_HYCOM/$CONFIGNAME/relax/$IEXPT`:

> ```rmu.[a,b]```: 2D relaxation time scale definition files for physics [unit].
> 
> ```rmutr.[a,b]```: 2D relaxation time scale definition files are for biogeochemistry [unit].
> 
> ```thkdf4.[a,b]```: 2D biharmonic diffusion coeffient definition files for thickness [unit].
> 
> ```veldf4.[a,b]```: 2D biharmonic diffusion coeffient definition files for velocity [unit].

The layer diffusion files are generated by ```$BIN_HYCOM/hycom_creat_thkdf4_veldf4.py```. 

Before we run the script, we need to prepare model configuration files first:
```bash
set -u
cd $WORK_HYCOM/$CONFIGNAME/relax/$IEXPT/SCRATCH
ln -sf $WORK_HYCOM/$CONFIGNAME/topo/regional.grid.a .
ln -sf $WORK_HYCOM/$CONFIGNAME/topo/regional.grid.b .
ln -sf $WORK_HYCOM/$CONFIGNAME/topo/regional.depth.a .
ln -sf $WORK_HYCOM/$CONFIGNAME/topo/regional.depth.b .
```

Next, modify parameters ```mnn``` and ```mxx``` in ```hycom_creat_thkdf4_veldf4.py``` for a specific configuration:
```bash
MNN=0.05
MXX=0.125
rmuwidth=20
```
for TP2 and:
```bash
MNN=0.02
MXX=0.125
rmuwidth=20(?)
```
for TP5 and run
```bash
set -u
cd $WORK_HYCOM/$CONFIGNAME/relax/$IEXPT/SCRATCH
cat "$BIN_HYCOM/hycom_creat_thkdf4_veldf4.py" | \
  sed -e "s/^\([[:space:]]*\)mnn[[:space:]]*=.*/\1mnn=$MNN/" \
      -e "s/^\([[:space:]]*\)mxx[[:space:]]*=.*/\1mxx=$MXX/" \
  > hycom_creat_thkdf4_veldf4.py
```

Now you can run:
```bash
set -u
cd $WORK_HYCOM/$CONFIGNAME/relax/$IEXPT/SCRATCH
python hycom_creat_thkdf4_veldf4.py $rmuwidth

cp thkdf4.a ../
cp thkdf4.b ../
cp veldf4.a ../
cp veldf4.b ../
```

The diffusion files are ready under `$WORK_HYCOM/$CONFIGNAME/relax/$IEXPT` now.

### Creating job settings

In this guide, we assume we create a job setting for 1 week hycom run:
```bash
DATE_START="2016-08-27"
DATE_END="2016-09-03"
```

Before hycom job is submitted, we need to prepare external files specific to the given period of integration, reatrst files, atmospheric forcing files and nesting files.

#### *Prepare restart files*

Restart files for TP2 hindcast run can be found under:
```bash
DIR_RST=/nird/datalake/NS9481K/shuang/TP2_output/expt_02.6
```
then copy restart files to `$WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT/data` folder:
```bash
set -u

cd $WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT/
#[ -d data ] && rm -rf data
mkdir -p data/cice
#
DATE_RESTART=$(date -d "$DATE_START" +%Y_%j)
#
cd $WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT/data
echo "cd $WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT/data"
restart_afile=${DIR_RST}/restart/restart.${DATE_RESTART}_00_0000.a
restart_bfile=${DIR_RST}/restart/restart.${DATE_RESTART}_00_0000.b
[ -f "$restart_afile" ] && cp $restart_afile .
[ -f "$restart_afile" ] && echo "   cp $restart_afile ."
[ -f "$restart_bfile" ] && cp $restart_bfile .
[ -f "$restart_bfile" ] && echo "   cp $restart_bfile ."
#
cd $WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT/data/cice
echo "cd $WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT/data/cice"
restart_icefile=${DIR_RST}/cice/iced.${DATE_START}-00000.nc
[ -f "$restart_icefile" ] && cp ${restart_icefile} .
[ -f "$restart_icefile" ] && echo "   cp ${restart_icefile} ."
#
cd $WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT
echo "cd $WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT"
mkdir -p $WORK_HYCOM_DAT/cice
[ ! -f $WORK_HYCOM_DAT/$restart_afile ] && cp $restart_afile $WORK_HYCOM_DAT
[ ! -f $WORK_HYCOM_DAT/$restart_bfile ] && cp $restart_bfile $WORK_HYCOM_DAT
[ ! -f $WORK_HYCOM_DAT/cice/$restart_icefile ] && cp $restart_icefile $WORK_HYCOM_DAT/cice
```

#### *Prepare off-line nesting files*

With offline nesting option on ( in `blkdat.in`), we need to prepare nesting files. Off-line nesting files for TP2 hindcast run can be found under:
```bash
DIR_NST=/nird/datalake/NS9481K/shuang/nest/TP2_expt023
```

Since the nesting files for `TP2_expt023` are archived files (`*.tar.gz`), we use the following script to extract and copy neccessary files:
```bash
set -u

# create nesting files folder
WORK_HYCOM_NST=$WORK_HYCOM/$CONFIGNAME/nest/$IEXPT
mkdir -p $WORK_HYCOM_NST

# copy nesting files from archive
cd $WORK_HYCOM_NST
for (( d=$(date -u -d "$DATE_START" +%s); d<=$(date -u -d "$DATE_END" +%s); d+=86400 )); do
    DATE_NOW=$(date -u -d "@$d" +%Y-%m-%d)
    YYYY=$(date -d "$DATE_NOW" +%Y)
    DOY=$(date -d "$DATE_NOW" +%j)

    # list of nesting files
    afile=archv.${YYYY}_${DOY}_00.a
    bfile=archv.${YYYY}_${DOY}_00.b
    afile_fabm=archv_fabm.${YYYY}_${DOY}_00.a
    bfile_fabm=archv_fabm.${YYYY}_${DOY}_00.b
    #echo $DATE_NOW $afile

    if [ -f ${DIR_NST}/$afile ]; then
        cp ${DIR_NST}/$afile .
        cp ${DIR_NST}/$bfile .
        cp ${DIR_NST}/$afile_fabm .
        cp ${DIR_NST}/$bfile_fabm .
    fi

    # extract nesting files from archive

    if [ -d ${DIR_NST}/tar_files ]; then
        for f in ${DIR_NST}/tar_files/archv.${YYYY}_*.tar.gz; do
            range=$(echo "${f#${DIR_NST}/tar_files/archv.${YYYY}_}" | sed 's/\.tar\.gz//')
            #echo $range
            start=$((10#${range%_*}))
            end=$((10#${range#*_}))
            if (( DOY >= start && DOY <= end )); then
                tarfile=${DIR_NST}/tar_files/archv.${YYYY}_$range.tar.gz
            fi
        done
        if [ -f $tarfile ]; then
            echo "tar -xzf $tarfile $afile $bfile"
            tar -xzf $tarfile $afile $bfile
        fi

        for f in ${DIR_NST}/tar_files/archv_fabm.${YYYY}_*.tar.gz; do
            range=$(echo "${f#${DIR_NST}/tar_files/archv_fabm.${YYYY}_}" | sed 's/\.tar\.gz//')
            #echo $range
            start=$((10#${range%_*}))
            end=$((10#${range#*_}))
            if (( DOY >= start && DOY <= end )); then
                tarfile_fabm=${DIR_NST}/tar_files/archv_fabm.${YYYY}_$range.tar.gz
            fi
        done
        if [ -f $tarfile_fabm ]; then
            echo "tar -xzf $tarfile_fabm $afile_fabm $bfile_fabm"
            tar -xzf $tarfile_fabm $afile_fabm $bfile_fabm
        fi
    fi
done
```

After this, nesting files can be found under `$WORK_HYCOM/$CONFIGNAME/nest/$IEXPT`. How to create off-line nesting files, see Shuang's note.

#### *Prepare atmospheric forcing files (⚠️Skip when using srjob.sh to submit a job)*

Run the script:
```bash
set -u
cd $WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT
START=${DATE_START}"T00:00:00"
END=${DATE_END}"T00:00:00"
$BIN_HYCOM/atmo_synoptic.sh era5+lw $START $END
```
which will create a set of forcing files:

```bash
airtmp.[ab]
dewpt.[ab]
mslprs.[ab]
precip.[ab]
radflx.[ab]
shwflx.[ab]
vapmix.[ab]
wndewd.[ab]
wndnwd.[ab]
```
under ```$WORK_HYCOM/$CONFIGNAME/force/synoptic/$IEXPT```.

Since name of forcing files does not contain date marker, it is convenient to leave a stamp under the forcing directory:
```bash
set -u
cd $WORK_HYCOM/$CONFIGNAME/force/synoptic/$IEXPT
rm stamp_*
GDATE_START=$(date -d "$START" +"%Y%m%d")
GDATE_END=$(date -d "$END" +"%Y%m%d")
touch stamp_${GDATE_START}-${GDATE_END}
```
For example, this stamp file can be used in SLRUM job script to avoid duplication of forcing file generation.

#### *Prepare SCRATCH folder (⚠️Skip when using srjob.sh to submit a job)*

Then we are ready to prepare SCRATCH folder:

```bash
set -u
cd $WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT
../bin/expt_preprocess.sh $START $END $INITFLG
```
If successful, you see message:
```bash
No fatal errors. Ok to start model set up in ..
```
and you are ready to submit a hycom job.

### Job submission

We use default job script `srjob.sh` under `$WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT` to submit a hycom job. First, fix the original `srjob.sh` by
```bash
sed -i \
  -e 's|expt_postprocess.sh|bin/expt_postprocess.sh|g' \
  -e '/^done$/d' \
  srjob.sh
chmod +x srjob.sh
```
and, modify `srjob.sh` based on experiment settings (`START`, `END` and `INITFLG`) accordingly:
```bash
set -u

TIME="00:30:00"
START=${DATE_START}"T00:00:00"
END=${DATE_END}"T00:00:00"

cd $WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT
[ ! -r srjob.sh.backup ] && cp srjob.sh srjob.sh.backup

sed \
  -e "s|^START=.*|START=\"$START\"|" \
  -e "s|^END=.*|END=\"$END\"|" \
  -e "s|^INITFLG=.*|INITFLG=\"$INITFLG\"|" \
  -e "s|^#SBATCH --time=.*|#SBATCH --time=\"$TIME\"|" \
  srjob.sh.backup > srjob.sh
chmod +x srjob.sh
```
and submit `srjob.sh`:
```bash
set -u
cd $WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT
sbatch srjob.sh
```
Note TP2 7 days integration job takes about 10 min CPU time on Betzy.

After the submission, use command:
```bash
squeue -u $USER
```
to check status of the submitted job.

If job is successfuly finished, restart and daily mean files are moved to `$WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT/data` folder. The successful completion can be confirmed by
```bash
set -u
cd $WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT
cat log/hycom.stop
```
If it returns `GOODRUN`, job is successfuly finished. If not, check log records under `log` to see what went wrong.

There is useful tool to create a sample job script on Sigma2 HPC: https://open.pages.sigma2.no/job-script-generator/.

### Visualization

For plotting hycom output file under `$WORK_HYCOM/$CONFIGNAME/expt_$NEWEXPERIMENT/data` with Python, see sample Jupyter notebooks [`plot_2D_TP2_temp.ipynb`](https://github.com/nansencenter/TP2_setup/blob/main/plot_2D_TP2_temp.ipynb) for 2D mapping and [`plot_section_TP2.ipynb`](https://github.com/nansencenter/TP2_setup/blob/main/plot_section_TP2.ipynb) for vertical section mapping.

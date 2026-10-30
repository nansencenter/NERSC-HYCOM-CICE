### Example nesting file generation

Sample year: 2005, member: 4

```bash
# I haven't yet written a script to automate it. Whenever you create nesting files, make sure to copy them to the relevant folder in "/nird/datalake/NS9481K/www/NorCPM_nesting/YEAR/MEMBER/"
# After they are safely copied, remove archv* files in "TP2a0.10/nest/exp_number/"
```

```bash
# clone nesting file downloading scripts
git clone https://github.com/lilleannette/Script-tools.git
cd Script-tools/SEACLIM/
git checkout develop

# clone NERSC-HYCOM-CICE
git clone https://github.com/nansencenter/NERSC-HYCOM-CICE.git
git checkout test-tillnest2
```

```bash
cd Script-tools/SEACLIM/
# the following will download input year and member
# it will take some time, I suggest using "screen" to avoid disconnections
# ignore warnings
# potetially you can do this once with shell script
# e.g. loop over years and members
./Preproc_norcpm_ocn.sh 2005 4
# when that completes, you should have raw CPM files under "Script-tools/SEACLIM/noresm2-mm-seaclim_hindcast/noresm2-mm-seaclim_hindcast_20051101_mem004/"
```

```bash
# At this point CPMa1.00 and TP2a0.10 top folders should exist
cd CPMa1.00/expt_01.0
cp ../../Script-tools/SEACLIM/noresm2-mm-seaclim_hindcast/noresm2-mm-seaclim_hindcast_20051101_mem004/*merged* ../Nesting_files/

# The following script selects all merged files in Nesting_files. Unless you give them specific years and members, it will try to create nesting files for all norcpm merged files. You already have a copy of those merged files elsewhere, so it is safe to remove the files from previous nesting creation

# This one creates the nesting files
../bin/Nesting_noresm/cpm_to_hycom.sh ../../TP2a0.10/expt_08.0/ ../Nesting_files/*2005*mem004*merged_*.nc
```

```bash
# Assuming the nesting files are ready, the next step is to calculate Montgomery potential using a restart file. I have one for 1993, but maybe if you get all the start dates from Annette, make a list, and ask Vijith to send restart files for November 1st for each year that the ensemble will be started.
cd ../../TP2a0.10/expt_08.0/
# Assuming you have placed the restart file in the experiment folder
python ../bin/calc_montg1.py ../nest/080/archv.20[01]*.a  ./restart.2005_305_00_0000.a  ./
# Now you should a number of archive files in the experiment folder.
# The next step is to rename them. It starts from January, but the files you have start from November, so some first iterations will give an error.
./Renames.sh  2005 2010 # this code renames the archive files, and moves the original ones to Orig folder.

# The output archive files are your final product. You can copy them to Nird folder for sharing.
# As of writing this, the directory for sharing files are:
/nird/datalake/NS9481K/www/NorCPM_nesting/ 
```

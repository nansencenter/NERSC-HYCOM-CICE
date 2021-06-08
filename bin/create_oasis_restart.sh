# How to make restart files for the coupled HYCOM-neXtSIM.

# 1. Make sure NNOREST is set to tru in namecouple.
# 2. Make sure each variable exchange option is set to EXPOUT in namcouple.
# 3. Make sure ice.nc and ocean.nc are not copied to the SCRATCH directory during the run.
# 4. Run the model for a short time (2 time-steps is enough - it will probabaly crash).
# 5. Create ocean.nc using this script and make sure it is copied to SCRATCH in the next run.
# 5. Run the model again for a short time (2 time-steps is enough).
# 6. Create ice.nc using this script and make sure both ice.nc and ocean.nc are copied to SCRATCH during the run.
# 7. The restart files should be ok for initializing the model.  New ice.nc and ocean.nc are created at the end of each run.


#! /bin/bash
module load CDO/1.9.5-intel-2018b
module load NCO/4.7.9-intel-2018b
# make the sea ice restart.
# merge file
cdo -O merge I_taux_nxtsim_08.nc I_tauy_nxtsim_09.nc I_fwflux_nxtsim_10.nc I_rsnos_nxtsim_11.nc I_rsso_nxtsim_12.nc I_sfi_nxtsim_13.nc I_sic_nxtsim_14.nc I_psl_nxtsim_15.nc I_wspeed_nxtsim_16.nc  ice.tmp1.nc 

# select the second time step
cdo seltimestep,2 ice.tmp1.nc ice.tmp2.nc
#remove very high values
mv ice.tmp2.nc ice.nc

# make the ocean restart.
# remove very high values
cdo -expr,'O_SSTSST = ((O_SSTSST>10000.0)) ? (0.0) : O_SSTSST' O_SSTSST_hycomm_01.nc O_SSTSST.nc
cdo -expr,'O_SSSal = ((O_SSSal>10000.0)) ? (0.0) : O_SSSal' O_SSSal_hycomm_02.nc O_SSSal.nc
cdo -expr,'O_OCurx1 = ((O_OCurx1>10000.0)) ? (0.0) : O_OCurx1' O_OCurx1_hycomm_03.nc O_OCurx1.nc
cdo -expr,'O_OCury1 = ((O_OCury1>10000.0)) ? (0.0) : O_OCury1' O_OCury1_hycomm_04.nc O_OCury1.nc
cdo -expr,'O_E3T1st = ((O_E3T1st>10000.0)) ? (0.0) : O_E3T1st' O_E3T1st_hycomm_06.nc O_E3T1st.nc

# merge file
cdo -O merge O_SSTSST.nc O_SSSal.nc  O_OCurx1.nc O_OCury1.nc O_SSHght_hycomm_05.nc O_FraQsr_hycomm_07.nc O_E3T1st.nc ocean.tmp1.nc

# select the second time step
cdo seltimestep,2 ocean.tmp1.nc ocean.tmp2.nc
rm ocean.nc
ncwa -a time ocean.tmp2.nc ocean.nc



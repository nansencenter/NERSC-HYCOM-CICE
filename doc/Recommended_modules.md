###########################################################
# setting up the compiling software  environoment
  recommeneded on different HPC
#
# N.B.:
#     It would be updated with the time as the systemal compiler existed. 
#
###########################################################
#
# Fram
#
module restore system
module load NCL/6.6.2-intel-2018b
module load ESMF/7.1.0r-intel-2018b
module load FFTW/3.3.8-intel-2018b
module swap Python/3.6.6-intel-2018b Python/2.7.15-intel-2018b
module load ncview/2.1.7-intel-2018b
module load NCO/4.7.9-intel-2018b


#
# Surfsara
#
module purge
module load 2019
module load NCL
module load FFTW/3.3.8-intel-2018b
module load ncview/2.1.7-intel-2018b
module load CDO/1.9.5-intel-2018b
module swap Python/3.6.6-intel-2018b Python/2.7.15-intel-2018b
module swap Python/3.6.6-intel-2018b Python/2.7.15-intel-2018b
module swap GEOS/3.6.2-intel-2018b-Python-3.6.6 GEOS/3.6.2-intel-2018b-Python-2.7.15 
module swap GDAL/2.2.3-intel-2018b-Python-3.6.6 GDAL/2.2.3-intel-2018b-Python-2.7.15

For example:

to compiling hycom_cice: 

   compile_model.sh -m fram  -u ifort
.or.
   compile_model.sh -m surfsara  -u ifort


export SURFSARA_NETCDFF_ROOT=/sw/arch/RedHatEnterpriseServer7/EB_production/2019/software/netCDF/4.6.1-intel-2018b

N.B.: The shared project directory is /projects/0/bhao2 which should be difined in the related REGION.src



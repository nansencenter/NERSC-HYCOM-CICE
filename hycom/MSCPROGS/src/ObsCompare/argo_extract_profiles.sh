#!/bin/bash
# Routine extracts profiles for input files. It finds correct dates by checking
# dates of the files you specify. This involves checking locally (see ARGO_PATH below), 
# as well as looking in the FTP web site.
#
# Finally the vertical/horizontal interpolation is done, using the argocmp routine.
# 

# Check number of args - must be greater than or equal to 1
if [ $# -lt 1 ] ; then 
   echo "Program takes as input a model file of DAILY type"
   echo "It first retrieves the date from the input file name  and "
   echo "then tries to retrieve the argo profile data - first by "
   echo "looking at files in the local storage (env var ARGO_PATH), "
   echo "but if that fails it will look in the ifremer FTP site"
   echo "When the argo profile is downloaded, it is processed with "
   echo "The routine \"argocmp\" which extracts model data at locations"
   echo "of the argo profiles."
   echo 
   echo "The final  data is saved as individual profiles in text files"
   echo "under the directory \"ArgoCmp\". Also, if the optional argument"
   echo "-ncdump is specified, the model data, and levitus climatology "
   echo "will be dumped to the original netcdf file as well"
   echo
   echo "Usage:"
   echo "    $(basename $0) [-ncdump] files " 
   echo "where -ncdump is a optional argument"
   echo
   echo "NB: The main difference between this file and argo_extract_profiles_new.sh"
   echo "    is that this version only processes DAILY files"
   exit
fi


ncdump=""
if [ "$1" == "-ncdump" ] ;then
   ncdump="-ncdump"
fi

# This assumes routines are in the same location as this routine
basedir=$(dirname $0)


files=$@

#  Check for path of ARGO data
if [ "${ARGO_PATH}" == "" ] ; then
   echo "You should specify the path to ARGO profiles in environment variable"
   echo "ARGO_PATH. This is where the routine will for profiles first, and it"
   echo "is where downloaded profiles will be placed"
   exit 1
fi
## This location contins data cache - when not obtained from ifremer
#ARGO_PATH=/home/fimm/nersc/knutali/Data/Argo 

# Where argo data can be found - FTPDIR contains month/year as well (see below)
# NB: constrained to atlantic ocean for now. Change REGION below
REGION="atlantic_ocean"
FTPARGO=ftp.ifremer.fr
FTPDIRARGO="ifremer/argo/geo/$REGION/"


# For each file check that it is present in my Data dir
for i in $files; do

   # Check that this is a daily file. To be changed later on
   if echo $i | grep "DAILY" ; then

   rungen=$(echo $i | cut -c1-3)

   # Get date info from file name 
   thisrefyear=$(echo $i | sed "s/${rungen}DAILY_[0-9]\+_[0-9]\+_//" | sed "s/_.*//")
   thisjday=$(echo $i |    sed "s/${rungen}DAILY_[0-9]\+_[0-9]\+_[0-9]\+_//" |  sed  "s/\..*//")
   #echo $rungen  $thisrefyear $thisjday

   # Get actual date 
   act_date=$(jultodate $thisjday $thisrefyear 1 1)
   act_year=$(echo $act_date | cut -c1-4)
   act_month=$(echo $act_date | cut -c5-6)

   # This is the argo file name
   argofile=${act_date}_prof.nc
   argofilegz=${act_date}_prof.nc.gz

   echo $thisrefyear $thisjday $act_date $act_year $act_month
   #exit

   # Is the corresponding  file present in my data dir ?
   tmpgz=${ARGO_PATH}/geo/$REGION/${act_year}/${act_month}/${argofilegz}
   tmpdir=${ARGO_PATH}/geo/$REGION/${act_year}/${act_month}/
   if [ -f $tmpgz ] ; then
      gunzip < ${tmpgz} > $argofile
      echo "$argofile Retrieved from local backup"
   elif [ ! -f $tmpgz ] ; then

      # Try to retrieve from ifremer
ncftp << EOF
open $FTPARGO
cd $FTPDIRARGO/${act_year}/${act_month}/
bin
hash
prompt
get $argofile
bye
EOF

      # Copy to data to use later
      [ ! -d $tmpdir ] && mkdir -p $tmpdir
      [ ! -f $tmpgz  ] && gzip < $argofile > $tmpgz

      echo "$argofile:"
      if [ ! -f $argofile  ] ; then
         echo "Retrieve failed !"
      else
         echo "Retrieve successful"
      fi
   fi


   # The file SHOULD be here now
   if [ -f $argofile  ] ; then
      $basedir/argocmp $argofile $i $ncdump
   fi

   fi # Check for daily
done



# At this stage, we have a catalog Argocmp, with profiles from model and data.

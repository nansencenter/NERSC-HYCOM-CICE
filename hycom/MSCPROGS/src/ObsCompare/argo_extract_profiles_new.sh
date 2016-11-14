#!/bin/bash
# Routine extracts profiles for input files. it Finds correct dates by checking
# dates of the files you specify. This involves checking locally (see ARGO_PATH below), 
# as well as looking in the FTP web site.
#
# Finally the vertical/horizontal interpolation is done, using the argocmp routine.
# 

# Check number of args - must be greater than or equal to 1
if [ $# -lt 1 ] ; then 
   echo "Program takes as input a model file of DAILY or WEEKLY type"
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
   exit
fi

ncdump=""
if [ "$1" == "-ncdump" ] ;then
   ncdump="-ncdump"
   shift
fi

# This assumes routines are in the same location as this routine
basedir=$(dirname $0)

files=$@



# This location contins data cache - when not obtained from ifremer
#  Check for path of ARGO data
if [ "${ARGO_PATH}" == "" ] ; then
   echo "You should specify the path to ARGO profiles in environment variable"
   echo "ARGO_PATH. This is where the routine will for profiles first, and it"
   echo "is where downloaded profiles will be placed"
   exit 1
fi

# Where argo data can be found - FTPDIR contains month/year as well (see below)
REGION="atlantic_ocean"
FTPARGO=ftp.ifremer.fr
FTPDIRARGO=ifremer/argo/geo/$REGION/


# For each file check that it is present in my Data dir
for i in $files; do

   rungen=$(echo $i | cut -c1-3)


   # Check that this is a daily file. To be changed later on
   if echo $i | grep "DAILY" ; then

      rungen=$(echo $i | cut -c1-3)

      # Get date info from file name 
      thisrefyear=$(echo $i | sed "s/${rungen}DAILY_[0-9]\+_[0-9]\+_//" | sed "s/_.*//")
      thisjday=$(echo $i |    sed "s/${rungen}DAILY_[0-9]\+_[0-9]\+_[0-9]\+_//" |  sed  "s/\..*//")
      #echo $rungen  $thisrefyear $thisjday


      echo $thisrefyear $thisjday $act_date $act_year $act_month
      #exit

      let fstday=$thisjday
      let lstday=$thisjday

   elif echo $i | grep "${rungen}AVE_" ; then

      #echo "Not tested for AVE files"
      #exit

      thisrefyear=$(echo $i | cut -d "_" -f2 )
      thismonth=$(echo $i | cut -d "_" -f3 )
      thisweek=$(echo $i | cut -d "_" -f4 | sed "s/\..*//" )

      echo $thisyear
      echo $thismonth
      echo $thisweek

      numweeks=$( echo $thismonth $thisweek | perl -lane 'print 4*($F[0]-1)+$F[1];')

      fstday=$(echo $numweeks | perl -lane 'use POSIX qw(ceil floor); use List::Util qw(min max);  print max floor(365*($F[0]-1)/48),1')
      #lstday=$(echo $numweeks | perl -lane 'use POSIX qw(ceil floor); use List::Util qw(min max);  print min ceil(365*$F[0]/48),365')
      lstday=$(echo $numweeks | perl -lane 'use POSIX qw(ceil floor); use List::Util qw(min max);  print min floor(365*$F[0]/48),365')

      echo $numweeks $fstday $lstday
   fi # Check for daily

   let j=$fstday
   while [ $j -le $lstday ] ; do

      # Get actual date 
      act_date=$(jultodate $j $thisrefyear 1 1)
      act_year=$(echo $act_date | cut -c1-4)
      act_month=$(echo $act_date | cut -c5-6)

      # This is the argo file name
      argofile=${act_date}_prof.nc
      argofilegz=${act_date}_prof.nc.gz

      echo $argofile

      # Is the corresponding  file present in my data dir ?
      tmpgz=${ARGO_PATH}/geo/$REGION/${act_year}/${act_month}/${argofilegz}
      tmpdir=${ARGO_PATH}/geo/$REGION/${act_year}/${act_month}/
      if [ -f $tmpgz ] ; then
         gunzip < ${tmpgz} > $argofile
         echo "$argofile Retrieved from local backup"
      elif [ ! -f $tmpgz ] ; then

#      # Try to retrieve from ifremer
ncftp << EOF
open $FTPARGO
cd $FTPDIRARGO/${act_year}/${act_month}/
bin
hash
prompt
get $argofile
bye
EOF

         reslt=$?
         # Copy to data to use later
         [ ! -d $tmpdir ] && mkdir -p $tmpdir
         [ ! -f $tmpgz  ] && gzip < $argofile > $tmpgz

         echo "$argofile $?:"
         if [ ! -f $argofile  ] ; then
            echo "Retrieve failed !"
         else
            echo "Retrieve successful"
         fi
      fi


      # The file SHOULD be here now
      if [ -f $argofile  ] ; then
         # Copy to data to use later
         [ ! -d $tmpdir ] && mkdir -p $tmpdir
         [ ! -f $tmpgz  ] && gzip < $argofile > $tmpgz
         $basedir/argocmp $argofile $i $ncdump
      fi

      let j=$j+1
   done

done



# At this stage, we have a catalog Argocmp, with profiles from model and data.

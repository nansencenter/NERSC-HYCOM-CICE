#!/bin/bash
# 1)Retrieves DUACS file from FTP corresponding to hycom file date
# 2)Calls slacmp

usage=\
"This routine takes as input a daily or weekly (not yet) average file,
calculates the julian day relative to 1950-1-1, and retrieves the 
sla file corresponding to this date. It then calculates the  sla RMS
errors using the routine slacmp.

Usage: $(basename $0) hycom-daily-file
"

# This assumes slacmp is in same directory as this script
SLACMP=$(dirname $0)/slacmp

if [ $# -ne 1 ] ; then
   echo -e "$usage"
   exit 1
fi

# Check that this is a daily file. To be changed later on
fname=$1
rungen=$(echo $fname | cut -c1-3)
if echo $fname | grep "DAILY" > /dev/null ; then
   # Get date info from file name 
   thisrefyear=$(echo $fname | sed "s/${rungen}DAILY_[0-9]\+_[0-9]\+_//" | sed "s/_.*//")
   thisjday=$(echo $fname |    sed "s/${rungen}DAILY_[0-9]\+_[0-9]\+_[0-9]\+_//" |  sed  "s/\..*//")
   echo $thisrefyear $thisjday
   let fstday=$thisjday
   let lstday=$thisjday
elif echo $fname | grep "${rungen}AVE_" ; then
   thisrefyear=$(echo $fname | cut -d "_" -f2 )
   thismonth=$(echo $fname | cut -d "_" -f3 )
   thisweek=$(echo $fname | cut -d "_" -f4 | sed "s/\..*//" )
   echo $thisyear
   echo $thismonth
   echo $thisweek
   numweeks=$( echo $thismonth $thisweek | perl -lane 'print 4*($F[0]-1)+$F[1];')
   fstday=$(echo $numweeks | perl -lane 'use POSIX qw(ceil floor); use List::Util qw(min max);  print max floor(365*($F[0]-1)/48),1')
   lstday=$(echo $numweeks | perl -lane 'use POSIX qw(ceil floor); use List::Util qw(min max);  print min floor(365*$F[0]/48),365')
   echo $numweeks $fstday $lstday
   echo "Unsafe for weekly data - fix me (weekly data needs excact start/stop + a range of files)"
   exit 1
fi # Check for daily


# Convert days to rel day to 1950-1-1
jday0=$(datetojul $thisrefyear 1 1 1950 1 1)
let jday=$jday0+$fstday
jday=$(printf "%5.5d" $jday )
echo "Julian day rel 1950,1,1 = $jday (check me!)"

# Get sla data
duacsfile="msla_oer_merged_h_${jday}.nc.gz"
duacsfileuc=$(echo $duacsfile | sed "s/\.gz//")
FTPDUACS="ftp.cls.fr"
FTPPATH="/pub/oceano/AVISO/SSH/duacs/global/nrt/msla/merged/h/"

# Retrieve if not present
if [ ! -f $duacsfileuc ] ; then
ncftp <<EOF
open $FTPDUACS
cd $FTPPATH
bin
get $duacsfile
bye
EOF
gunzip $duacsfile
fi

# Run slacompare
echo "Running : $SLACMP $fname $duacsfileuc"
$SLACMP $fname $duacsfileuc

#!/bin/bash
# Sample script to generate ice drift using the "icedrift2" routine
#
# This is for an ensemble, using ICEDRIFT files and CERSAT file

usage=" Usage:
$(basename $0) [opional args] rungen bulletinyear bulletinday CERSATdriftfile
Optional arguments:
     -f number    --- First member to process (default=1  ) 
     -l number    --- Last  member to process (default=100) "


first=1
last=100
while getopts "f:l:" options; do
     case $options in
        f ) first=$OPTARG ;;
        l ) last=$OPTARG ;;
        * ) echo "$usage"
        exit 1;;
     esac
done
shift $(($OPTIND - 1))
if [ $# -ne 4 ] ;then 
   echo "$usage"
   exit 1
else
   rungen=$1
   byear=$2
   bday=$3
   CERSATfile=$4
fi





# Assumes ice drift routine is in same directory as this script
prog=$(dirname $0)
prog="$prog/icedrift2"

cnt=$first
while [ $cnt -le $last ] ; do


   # Sample -- the values here should be from CLI arguments
   #./icedrift2.sh RTH $byear $bday $CERSATfile $cnt
   $prog RTH $byear $bday $CERSATfile $cnt

   ret=$?

   if [ $ret -ne 0 ] ; then
      echo icedrift2.sh routine failed 
      exit 1
   fi

   let cnt=$cnt+1

done



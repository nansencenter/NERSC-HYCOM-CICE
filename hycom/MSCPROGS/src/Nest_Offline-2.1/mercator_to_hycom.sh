#!/bin/bash
progbase=/home/fimm/nersc/knutali/Progs/Setup/Nest_Offline
proghint=$progbase/hint_mercator
progvint=$progbase/vremap_mercator

usage=" Usage   $(basename $0)  -T restarttemplate  mercatorTfiles \n
Required arguments: -T restarttemplate \n
Optional arguments:  \n
        -T restarttemplate \n
        -R rungen \n
\n
restarttemplate is any valid restart file for the inner model\n
rungen is the three leter id used in restart files (default = xxx) \n
mercatorTfiles are the files for the T grid \n
Dates will be retrieved from the mercator files\n"

rungen="xxx"
while getopts "T:R:" options; do
     case $options in
        T ) template=$OPTARG
            [ -z $OPTARG ] && { echo -e $usage ; echo "-T needs a value" ; exit 1 ; };;
        R ) rungen=$OPTARG 
            [ -z $OPTARG ] && { echo -e $usage ; echo "-R needs a value" ; exit 1 ; };;
        * )  echo -e $usage
        exit 1;;
     esac
done

# Skip to files
shift $(($OPTIND - 1))

if [ ! $1 ] ; then
   echo -e $usage
   exit 1
fi

# Process mercator files
templatebase=$( echo $template | sed "s/\.[ab]//")
while  [ $1 ]  ; do 

   # Try to get dates from  mercator file
   datepart=$(echo $1 | sed "s/_gridT.*//"  | sed "s/.*_y//" )
   iyear=$( echo $datepart | cut -c1-4)
   imonth=$( echo $datepart | cut -c6-7)
   iday=$( echo $datepart | cut -c9-10)

   # Convert to julain day rel this year
   jday=$( datetojul $iyear $imonth $iday $iyear 1 1 | tail -c4)



   basename=$( echo $1 | sed "s/T\.nc//")
   echo $basename

   # This creates merchint.ab
   $proghint  $basename ||  exit 1

   # This creates mercvint.ab from template and merchint.ab
   $progvint $template ||  exit 1

   # NB -- assumes daily average
   newfilebase=${rungen}restart${iyear}_${jday}_12 

   # This finalizes things - creates restart files readable by hycom

   # First we must change the time header in the restart file
   line=$(head -n2 $templatebase.b |  tail -n1)
   fld1=$(echo $line |  cut -f1 -d "=")
   fld2=$(echo $line |  cut -f2 -d "=" | sed "s/^[ ]*//" )
   fld2=$(echo $fld2 | sed "s/[ ]*/ /" )
   fld21=$(echo $fld2 |  cut -f1 -d " ")
   fld22=$(echo $fld2 |  cut -f2 -d " ")



   hrefyear=1901
   jday2=$( datetojul $iyear $imonth $iday $hrefyear 1 1)
   jday2=$(echo $jday2 + .5 | bc )

   echo $line
   echo $fld1
   echo $fld2
   echo $fld21
   echo $fld22
   echo $jday2
   



   head -n1 $templatebase.b > $newfilebase.b 
   echo "$fld1 = $fld21   $jday2 " >>  $newfilebase.b 
   cat mercvint.b >> $newfilebase.b
   cp mercvint.a $newfilebase.a

   echo "New hycom files  $newfilebase.[ab] Created "
   echo "horizontally interpolated fields can be diagnosed from merchint.[ab]"



   shift
done

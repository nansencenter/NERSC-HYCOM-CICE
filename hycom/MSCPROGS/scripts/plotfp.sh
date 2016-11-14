#!/bin/bash
# Parse any .a file  and create postsript plots

# Get command line args
usage="
Generic routine for plotting fields from hycom-type .a files
Usage: $(basename $0) [-p plottype] [-s scale] time1 time2 file
Where:
   time1 is first index in file to plot
   time2 is last  index in file to plot
   file  is file to get array data from (ex forcing.airtmp.a)
  Optional:
   -p plottype : plottype is x11(screen), psl(landscape PS)  or psp (portrait PS)
   -s scale    : is a scaling factor

Example:
Plot  fields 1 to 100 from forcing.airtmp.a
   plotfp.sh -p x11 1 100 forcing.airtmp.a"

# Get optional argument
plottype=x11
scale=1
while getopts "p:s:" options; do
     case $options in
        p ) plottype=$OPTARG ;;
        s ) scale=$OPTARG ;;
        * ) echo -e "$usage"
        exit 1;;
     esac
done
shift $(($OPTIND - 1))

if [ $# -ne 3 ] ; then
   echo -e "$usage"
   exit 1 
fi
ptime=$1
ftime=$2
FILE=$3

[ "$plottype" != "x11" ] && [ "$plottype" != "psp" ] && [ "$plottype" != "psl" ] && \
echo "Supply plot type as x11 (screen), psp(postxcript portrait) or psl(postscript landscape" \
&& exit


# Parse hycom blkdat
[ ! -f  blkdat.input ]  && { echo "Input file blkdat.input does not exist "; exit 1 ; }
IDM=$(cat blkdat.input | grep idm | sed  "s/^[ ]*//" | cut -f1 )
JDM=$(cat blkdat.input | grep jdm | sed  "s/^[ ]*//" | cut -f1 )
KDM=$(cat blkdat.input | grep kdm | sed  "s/^[ ]*//" | cut -f1 )
echo IDM $IDM
echo JDM $JDM
echo KDM $KDM



# Do input files exist?
[ ! -f  $FILE ]  && { echo "Forcing file $FILE does not exist..." ;  exit 1 ; }

# Test for Routinedir
if [ "${HYCOM_ALL}" == "" ]; then
   echo "Environment variable HYCOM_ALL is empty"
   exit 1
fi

# Test for main routine
FP=${HYCOM_ALL}/plot/src/fp_$plottype 
if [ ! -f $FP ]; then
   echo "Cant find routine $FP"
   exit 1
fi

# Header for input file
TEMPLATE="FILE           
ATLb2.00
IDM     'idm   ' = longitudinal array size
JDM     'jdm   ' = latitudinal  array size
  1     'nperfr' = number of horizontal plots per frame
  0     'lalolb' = spacing of latitude/longitude labels
  0     'lalogr' = spacing of latitude/longitude grid over land (<0 land+sea)
  4     'loclab' = location of the contour label (1=upr,2=lowr,3=lowl,4=upl)
 13     'locbar' = location of the color bar     (1[0-4]=vert,2[0-4]=horiz)
  5     'kpalet' = palete (0=none,1=pastel,2=sst,3=gaudy,4=2tone,5=fc,6=ifc)
  1     'iorign' = i-origin of plotted subregion
  1     'jorign' = j-origin of plotted subregion
  0     'idmp  ' = i-extent of plotted subregion (<=idm; 0 implies idm)
  0     'jdmp  ' = j-extent of plotted subregion (<=jdm; 0 implies jdm)"
echo -e "$TEMPLATE" | sed "s/IDM/$IDM/" | sed "s/JDM/$JDM/" | sed "s/FILE/$FILE/" > plot.IN


# Skip to correct record
irec2=$ptime
while [ $irec2 -le $ftime ] ; do
   echo " $irec2	'nrec  '" >> plot.IN
   echo " Record $irec2 of  $FILE "        >> plot.IN
   echo "  $scale 	'qscale'" >> plot.IN
   echo " 0.0   'qq    '" >> plot.IN
   let irec2=$irec2+1
done
echo " -1	'nrec  '" >> plot.IN

$FP < plot.IN
if [ "$plottype" == "psp" ] || [ "$plottype" == "psl" ] ; then
   echo "File name is gmeta1.ps"
fi

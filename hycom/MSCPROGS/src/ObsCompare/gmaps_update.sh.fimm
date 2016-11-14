#!/bin/bash
# Routine processes argo profiles, and poroduces some plots, no stats yet
export PYTHONPATH=/home/fimm/nersc/knutali/opt/Python2.4/lib/python/
DATADIR=/work/laurentb/TOPAZ3_Diagnostics/
WORKDIR=/work/knutali/TOPAZ3/
RDIR=/home/fimm/nersc/knutali/Progs/TOOLS/ObsCompare/

TARGETDIR=/Network/Servers/sverdrup-e3.nersc.no/Data/topaz/Knut/html/Test/
TARGETUSER=topaz
TARGETHOST=nansen.nersc.no

SSH="ssh -i /home/fimm/nersc/knutali/.ssh/scriptkey_dsa"
SCP="scp -i /home/fimm/nersc/knutali/.ssh/scriptkey_dsa"

unique_bulletin_daily() {
   uni=$(ls -1 $@ | sed "s/_[0-9]\{4\}_[0-9]\{3\}.[ab]//g" | sed "s/.*...DAILY_//g" | uniq )
   echo $uni
}  
last_bulletin_daily() {
   tmp=$(unique_bulletin_daily $@)
   lastbulletin=$( echo $tmp | tr ' ' '\n' | tail -n1 )
   echo $lastbulletin
}  

#############################################################################################

cd $WORKDIR




# List of vital files
datafilelist="regional.grid.a regional.grid.b regional.depth.a regional.depth.b grid.info depths800x880.uf Data/latlon.dat" 

# Copy vital files here
for i in $datafilelist ; do
   [ ! -f $i ] && cp $DATADIR/$i $i
done

# Clean up argocmp dir
[ -d ArgoCmp ] && rm -rf ArgoCmp/*

lastbulletin=$( last_bulletin_daily $DATADIR/DAILY_F07/*DAILY* )
year=$(echo $lastbulletin | cut -f1 -d "_")  
jday=$(echo $lastbulletin | cut -f2 -d "_") 
bdate=$(jultodate $jday $year 1 1)
#echo $lastbulletin $year $jday $bdate

# Test for existing bulletindate file
echo "Last bulletin was at $bdate "
[ -f lastbulletin.gmaps ] && [ "$(cat lastbulletin.gmaps)" == "$bdate" ] && \
   { echo "This bulletin has already been processed (delete lastbulletin.gmaps if you disagree)" ; exit 1 ; }




# Copy Forecast07 files here & process them
dailyfiles=$DATADIR/DAILY_F07/*DAILY_${lastbulletin}*.[ab]
for i in $dailyfiles ; do
   echo  $i
   [ ! -f $(basename $i) ] && cp $i .

   # Start processing files
   $RDIR/argo_extract_profiles_new.sh  $(basename $i)
done

## copy argomp files to ArgoCmp/tmp
mkdir ArgoCmp/tmp
find ArgoCmp/[0-9]* -type f -exec cp {} ArgoCmp/tmp/ \;

# Descend into that dir and produce line plots
cd ArgoCmp/tmp/
python2.4 $RDIR/plotargocmp.py -dAgg

# Convert all those pngs to thumbnails - uses imagemagick
for i in *N.png ; do
   convert -geometry 200x150 $i $(echo $i | sed "s/\.png/\_thumb.png/")
done

# Now copy to topaz@nansen :-)
basedir=TOPAZ3_${bdate}
echo base dir is  $basedir

$SSH $TARGETUSER@$TARGETHOST "mkdir $TARGETDIR/$basedir ; chmod 755 $TARGETDIR/$basedir"
$SCP *.png $TARGETUSER@$TARGETHOST:$TARGETDIR/$basedir
$SSH $TARGETUSER@$TARGETHOST "find $TARGETDIR/$basedir/ -type f -exec chmod 644 {} \;"

cd $WORKDIR
echo $bdate > lastbulletin.gmaps

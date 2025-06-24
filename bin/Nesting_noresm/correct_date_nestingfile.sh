#! /bin/bash                                                                                                                     
# Beacuse the 15th of each month will not the the exact date HYCOM wants,  
# nesting files need to be renamed with that date.
# 
# To use from the experiment folder
# bin/correct_date_nestingfile.sh start_year end_year yrflag nesting_frequency

start_year=$1
end_year=$2
yrflag=$3  #from blkdat
nstfq=$4   #from blkdat
source EXPT.src
source ../REGION.src
cd ../nest/$E
syear=${start_year}#2095
eyear=${end_year} #2101
sjday=$(((($syear-1901)*1461)/4 + 30))
ejday=$(((($eyear-1901)*1461)/4 - 60))
echo $sjday $ejday
#echo $yrlag $nestfq $nestfq $sjday $ejday | $NHCROOT/hycom/hycom_ALL/hycom_2.2.72_ALL/bin/hycom_nest_dates > correct_dates.txt


linecount=0
for ((year=syear; year<eyear; year+=1)); do
   for sday in 016 046 075 106 136 167 197 228 259 289 320 350 ; do
       linecount=$(($linecount + 1))
       echo $linecount
       dstr=`sed -n ${linecount}p correct_dates.txt`
       echo $sday $dstr
       echo cp archv.${year}_${sday}_00.a archv.$dstr.a
       cp archv.${year}_${sday}_00.a archv.$dstr.a                                                                              
       cp archv.${year}_${sday}_00.b archv.$dstr.b                                                                              
       cp archv_fabm.${year}_${sday}_00.a archv_fabm.$dstr.a
       cp archv_fabm.${year}_${sday}_00.b archv_fabm.$dstr.b
       echo "mv"                                                                                                                  
       mv archv.${year}_${sday}_00.* Orig/                                                                                     
       mv archv_fabm.${year}_${sday}_00.* Orig/
   done
done

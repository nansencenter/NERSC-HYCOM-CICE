#! /bin/bash

syear=2034
eyear=2101
sjday=$(((($syear-1901)*1461)/4 + 30))
ejday=$(((($eyear-1901)*1461)/4 - 60))
echo $sjday $ejday
echo 3 30.4375 30.4375 $sjday $ejday | ~/HYCOM-CICE/NERSC-HYCOM-CICE/hycom/hycom_ALL/hycom_2.2.72_ALL/bin/hycom_nest_dates > correct_dates.txt

linecount=0
for ((year=syear; year<eyear; year+=1)); do
   for sday in 016 046 075 106 136 167 197 228 259 289 320 350 ; do
       linecount=$(($linecount + 1))
       echo $linecount
       dstr=`sed -n ${linecount}p correct_dates.txt`
       echo $sday $dstr
       echo cp archv.${year}_${sday}_00.a archv.$dstr.a
#       cp archv.${year}_${sday}_00.a archv.$dstr.a
#       cp archv.${year}_${sday}_00.b archv.$dstr.b
       cp archv_fabm.${year}_${sday}_00.a archv_fabm.$dstr.a
       cp archv_fabm.${year}_${sday}_00.b archv_fabm.$dstr.b       
#       echo mv
#       mv archv.${year}_${sday}_00.* Montg/
       mv archv_fabm.${year}_${sday}_00.* Orig/
   done
done

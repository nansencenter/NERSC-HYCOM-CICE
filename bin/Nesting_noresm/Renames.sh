#! /bin/bash
mkdir -p Orig
syear=$1
eyear=$2
sjday=$(((($syear-1901)*1461)/4 + 30))
ejday=$(((($eyear-1901)*1461)/4 - 60))
echo $sjday $ejday
echo 3 30.4375 30.4375 $sjday $ejday | /cluster/work/users/cagyum/NERSC-HYCOM-CICE/hycom/hycom_ALL/hycom_2.2.72_ALL/bin/hycom_nest_dates > correct_dates.txt

linecount=0
for ((year=syear; year<eyear; year+=1)); do
   for sday in 015 046 074 105 135 166 196 227 258 288 319 349 ; do
       
       existing=$sday
       if (( (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0) )); then
           # If it's a leap year AND we are past February, shift the day by 1
           if [ "$sday" -gt 046 ]; then
               existing=$(printf "%03d" $((10#$sday + 1)))
           fi
       fi
       
       linecount=$(($linecount + 1))
       echo $linecount
       dstr=`sed -n ${linecount}p correct_dates.txt`
       echo $sday $dstr
       echo "cp archv.${year}_${existing}_00.a archv.$dstr.a"
       cp archv.${year}_${existing}_00.a archv.$dstr.a
       cp archv.${year}_${existing}_00.b archv.$dstr.b
#       cp archv_fabm.${year}_${sday}_00.a archv_fabm.$dstr.a
#       cp archv_fabm.${year}_${sday}_00.b archv_fabm.$dstr.b       
       mv archv.${year}_${existing}_00.* Orig/
   done
done


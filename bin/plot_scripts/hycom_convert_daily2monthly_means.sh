#!/bin/bash
#SBATCH --job-name=TP6_deaily_2_monthlymean
#SBATCH --time=29:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH -o  ./Hycmean_out_.out 

# used for convert the daily archm file into monthly or yearly
#
# Based on the input year to search the related archm files 
# To create the basic information file as the inputs for hycom_mean
#
# Ex.: hycommean.sh <????>
# EX.: sbatch ./hycom_convert_archm_2monthly_means.sh year
# Output: TMP9999_999.a[b]
#
exec 1> ./Hycom_convert_archm2monthly_out.log 2>&1
# Give notice to SMS that job has started
echo "hyc2proj job STARTED" > ./Hycmean.start_msg
#
if [ $# -ne 1 ]
then
  	echo "${0}:Usage: ${0} refyear"
   exit 1
fi
# Give notice to SMS that job has started
echo "hyc2proj job STARTED" > ./Hycmean.start_msg
rid="TP6"
#rid="NM8"
#iy=$1
#export mnth=10
rfyy=$1
nxyy=$(( ${rfyy} + 1 ))
echo "---> ${rfyy}"
echo "---> ${nxyy}"

#---------------------First create links to the daily files--------------------------------
for ((year=${rfyy}; year<${nxyy} ; year+=1)); do
   #
 totday=366
 for (( day=1; day<=${totday}; day+=1 )); do
       Jsday=`echo -n 00${day} | tail -3c`
       Sdate0=`/home/sm_alfal/sea/TOPAZ6/NERSC-HYCOM-CICE/hycom/MSCPROGS/bin/jultodate ${Jsday} ${year}  1 0`
       yy=`echo ${Sdate0:0:4}`
       mm=`echo ${Sdate0:4:2}`
       dd=`echo ${Sdate0:6:2}`
       #
       echo " renaming ... ${rid}archm.${year}_${Jsday}_12.a  ${rid}archm${yy}_${mm}_${dd}.a "
       echo " renaming ... ${rid}archm.${year}_${Jsday}_12.b  ${rid}archm${yy}_${mm}_${dd}.b "
       #
       ln -sf ${rid}archm.${year}_${Jsday}_12.a  ${rid}archm${yy}_${mm}_${dd}.a 
       ln -sf ${rid}archm.${year}_${Jsday}_12.b  ${rid}archm${yy}_${mm}_${dd}.b 
 done
done


#---------------------Now compute monthly means and std--------------------------------
Progdir=/home/sm_alfal/sea/TOPAZ6/NERSC-HYCOM-CICE/hycom/hycom_ALL/hycom_2.2.72_ALL/meanstd/src
#
for (( iy=${rfyy}; iy<${nxyy} ; iy+=1)); do
#for (( iy=2009; iy<2010; iy+=1 )); do
   #
   rm dy2mn.txt
   export kk=`grep "'kdm   ' =" blkdat.input | awk '{printf("%2d", $1)}'`
   #
   for (( mm=1; mm<13; mm+=1 )); do
      smm=`echo -n 0${mm} | tail -2c`
      echo "compute hycmean: monthly mean for month -----------> " ${smm}
      rm list.txt
      #
      if [ ! -r list.txt ]; then
        #ls NM8archm*${iy}_${smm}_*.a >list.txt
        ls ${rid}archm*${iy}_${smm}_*.a >list.txt
      fi
      Fn=`cat list.txt | wc -l`
      echo "No. of days =====--> " $Fn
      # 
      #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      echo "compute hycmean: monthly mean        : month -----------> " ${smm}
      if [ -f dy2mn.txt ]; then
        echo "File dy2mn.txt found remove it!"
        rm dy2mn.txt
      fi 
      echo ${kk} ' kk' > dy2mn.txt
      echo '1  single' >> dy2mn.txt
      echo '0  meansq' >> dy2mn.txt
      echo "${Fn} narchs" >> dy2mn.txt
      cat list.txt >> dy2mn.txt
      echo '0 narchs' >> dy2mn.txt
      echo "${rid}monthly_${iy}_${smm}" >> dy2mn.txt
      cat dy2mn.txt
      if [ -r dy2mn.txt ]; then
        ${Progdir}/hycom_mean < dy2mn.txt
      fi
      #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      echo "compute hycmean: monthly mean squire  : month -----------> " ${smm}
      if [ -f dy2mn.txt ]; then
        echo "File dy2mn.txt found remove it!"
        rm dy2mn.txt
      fi 
      echo ${kk} ' kk' > dy2mn.txt
      echo '1  single' >> dy2mn.txt
      echo '1  meansq' >> dy2mn.txt
      echo "${Fn} narchs" >> dy2mn.txt
      cat list.txt >> dy2mn.txt
      echo '0 narchs' >> dy2mn.txt
      echo "${rid}monthly_meansq_${iy}_${smm}" >> dy2mn.txt
      # 
      if [ -r dy2mn.txt ]; then
        ${Progdir}/hycom_mean < dy2mn.txt
      fi
      #---------------------------------------------------------------------
   done
done

echo "hyc2proj job STARTED" > ./Hycmean.end_msg
#


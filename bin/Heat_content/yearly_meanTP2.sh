

Sourdrt='/nird/datalake/NS9481K/shuang/seaclim/output/old/3D/PHY/'

# archv_20080108_12.nc

vars="thetao so"

for ii in ${vars}; do
   echo $ii
   for jj in `seq 1993 2024`; do
      Fclim=TP2${ii}_archv${jj}.nc
      N0=$(ls ${Sourdrt}3DT-${ii}/TOPAZ2_1m-m_3DT-${ii}_${jj}??-${jj}??.nc 2>/dev/null | sed -n '$=')
      if [ ! -s ${Fclim} -a $N0 -eq 12 ]; then
         ncea -Oh ${Sourdrt}3DT-${ii}/TOPAZ2_1m-m_3DT-${ii}_${jj}??-${jj}??.nc ${Fclim}
      fi
   done
done



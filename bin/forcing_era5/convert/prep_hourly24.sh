
# Ensure the procesed NC file has a time variable with hourly unit

if [ $# -eq 1 ]; then
   File0=$1
   echo ${File0}
   # extract the time steps into hourly .nc
   J0=$(datetojul 1970 1 1 1900 1 1)
   (( Jh = $J0 * 24 ))
   FileT=time_${File0#*hourly.}
   ncks -Oh -v valid_time ${File0} ${FileT} 
   ncap2 -Oh -s "time=${Jh}+valid_time/3600." -s "time@calendar=\"gregorian\"" -s "time@units=\"hours since 1900-01-01 00:00:00\"" ${FileT} tmp.${FileT}
   ncks -Oh -C -x -v valid_time tmp.${FileT} ${FileT}
   ncks -Oh -3 ${FileT} tmp.${FileT}
   ncrename -Oh -d valid_time,time tmp.${FileT}
   ncks -Oh -4 -L 4 tmp.${FileT} ${FileT}
   rm tmp.${FileT}

   # append the time steps into hourly .nc
   File1=pre_${File0#*hourly.}   
   if [ ! -s ${File1} ]; then
      ncks -Oh -C -x -v valid_time,expver ${File0} ${File1}
   fi
   ncrename -Oh -d valid_time,time ${File1}
   [ ! -s ${File1} ] && { echo "No ${File1}";  exit 1; }
   if [ -r ./hourly ]; then
      ncks -Ah -v time ${FileT} ${File1}
      ncks -4 -L 4 ${File1} ./hourly/hourly.${File0#*hourly.}
      [ -s ${File1} ] && rm ${File1}
   fi

else
   echo "Missing the input NC file like *hourly.????_yyyy.nc"
   exit 0
fi

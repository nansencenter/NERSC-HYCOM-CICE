
# Concatenate two netcdf files without missing time steps
# Usage:

if [ ! -s env.drt ]; then
   echo $(pwd)
   echo "Missing the directory definition from env.drt"
   exit 1
fi

if [ $# -gt 0 ]; then
   yy=$1
else
   yy=2026
fi

source env.drt

cd ${Rundir}

# diretory check and clean
echo "diretory check and clean ... "
dirts="6h_00 new"
for ii in ${dirt}; do
  if [ ! -s $ii ]; then
     mkdir $ii
     if [ ${#ii} -gt 3 ]; then
        cp ${ERA5proj}/6h.*_${yy}.nc ${ii}/. 
     fi
  else
     rm ${ii}/*
  fi
done

#Sour2='/cluster/projects/nn9481k/ERA5_6h/'
Sour1='./6h_00/'
Sour2='./6h/'
Outdrt='./new'

for i in $(ls ${Sour2}6h.*.nc); do
   j=${Sour1}6h.${i#*6h.}
   k=${Outdrt}/6h.${i#*6h.}
   if [ -s $i -a -s $j -a ! -s $k ]; then
      echo ""
      echo " preparing $k ... "
      ncrcat -Oh $j $i $k
   fi
done

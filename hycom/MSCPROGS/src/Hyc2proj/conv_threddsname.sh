for i in WOA2005*.nc; do 
   secnum=$(ncdump $i | grep station_group_number | sed "s/.*=[ ]*//" |\
      sed "s/[ ]*;//");

   secname=$(ncdump $i | grep station_group_name | sed "s/.*=[ ]*//" |\
      sed "s/[ ]*;//");

   secnum=$(printf "%2.2d" $secnum)
   echo "$secnum --> $secname  -->  $i"

   #copy file to new name

   if echo $secname | grep MOORINGS ; then
      cp $i moorings.nc
   else
      cp $i section${secnum}.nc
   fi
done


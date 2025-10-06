
# dealing with the predefined-year's field (1-15)
#
if [ ! $# -eq  2 ]; then
   echo "Wrong input for $0"
   echo "Example: "
   echo "       do_6hourly.sh yyyy nn"
   exit 0
fi

# need the default input direectory 
[ ! -r ./hourly ] && { echo missing the hourly directory under $pwd; exit 0; } 
[ ! -r ./6h ] && mkdir 6h     # for output 

Indir=./hourly

iyy=$1
ii=$2

case ${ii} in
    1) 
      varname='10U'  
      oldname='u10'  
      par0=128
      par1=165 ;;
    2)
      varname='10V' 
      oldname='v10'  
      par0=128
      par1=166 ;;
    3)
      varname='2D'  
      oldname='d2m'  
      par0=128
      par1=168 ;;
    4)
      varname='2T'  
      oldname='t2m'  
      par0=128
      par1=167 ;;
    5)
      varname='BLH'  
      oldname='blh'  
      par0=128
      par1=159 ;;
    6)
      varname='MSL'  
      oldname='msl'  
      par0=128
      par1=151 ;;
    7)
      varname='TCC'  
      oldname='tcc'  
      par0=128
      par1=165 ;;
    8)
      varname='CI'  
      oldname='siconc'  
      par0=128
      par1=31 ;;
    9)
      varname='SKT'  
      oldname='skt'  
      par0=128
      par1=235 ;;
   10)
      varname='SSR'  
      oldname='ssr'  
      par0=128
      par1=176 ;;
   11)
      varname='STR'  
      oldname='str'  
      par0=128
      par1=177 ;;
   12)
      varname='SSRD'  
      oldname='ssrd'  
      par0=128
      par1=169 ;;
   13)
      varname='STRD'  
      oldname='strd'  
      par0=128
      par1=175 ;;
   14)
      varname='RO'  
      oldname='ro'  
      par0=128
      par1=205 ;;
   15)
      varname='TP'  
      oldname='tp'  
      par0=128
      par1=164 ;;
esac 

Fini=6h.${varname}_${iyy}.nc
if [ -s ${Indir}/hourly.${varname}_${iyy}.nc -a ! -s ./6h/${Fini} ]; then
   echo $ii ' + ' ${ivarnam} ' - ' ${varname} 
   if [ ! -s ${Fini} ]; then
      # correct the attribute for oldname
      # ncatted -a att_nm, var_nm, mode, att_typ, att_val in.nc [out.nc]
      ncatted -Oh -a coordinates,${oldname},m,c,"time latitude longitude" ${Indir}/hourly.${varname}_${iyy}.nc 
      cdo splithour ${Indir}/hourly.${varname}_${iyy}.nc ${varname}_${iyy}.steps
      if [ $ii -le 9 ]; then
         cdo mergetime ${varname}_${iyy}.steps06.nc ${varname}_${iyy}.steps12.nc ${varname}_${iyy}.steps18.nc ${varname}_${iyy}.steps00.nc ${Fini}

      #cdo mergetime ${varname}_${iyy}.steps03.nc ${varname}_${iyy}.steps06.nc ${varname}_${iyy}.steps09.nc ${varname}_${iyy}.steps12.nc ${varname}_${iyy}.steps15.nc ${varname}_${iyy}.steps18.nc ${varname}_${iyy}.steps21.nc ${varname}_${iyy}.steps00.nc 3h.${varname}_${iyy}.nc 
      else
         if [ -s ${varname}_${iyy}.steps00.nc ]; then
            cdo -b F32 enssum ${varname}_${iyy}.steps00.nc ${varname}_${iyy}.steps01.nc ${varname}_${iyy}.steps02.nc ${varname}_${iyy}.steps03.nc ${varname}_${iyy}.steps04.nc ${varname}_${iyy}.steps05.nc 6h01.${varname}_${iyy}.nc 
            cdo -b F32 enssum ${varname}_${iyy}.steps06.nc ${varname}_${iyy}.steps07.nc ${varname}_${iyy}.steps08.nc ${varname}_${iyy}.steps09.nc ${varname}_${iyy}.steps10.nc ${varname}_${iyy}.steps11.nc 6h02.${varname}_${iyy}.nc 
            cdo -b F32 enssum ${varname}_${iyy}.steps12.nc ${varname}_${iyy}.steps13.nc ${varname}_${iyy}.steps14.nc ${varname}_${iyy}.steps15.nc ${varname}_${iyy}.steps16.nc ${varname}_${iyy}.steps17.nc 6h03.${varname}_${iyy}.nc 
            cdo -b F32 enssum ${varname}_${iyy}.steps18.nc ${varname}_${iyy}.steps19.nc ${varname}_${iyy}.steps20.nc ${varname}_${iyy}.steps21.nc ${varname}_${iyy}.steps22.nc ${varname}_${iyy}.steps23.nc 6h04.${varname}_${iyy}.nc 

	    if [ -s 6h04.${varname}_${iyy}.nc ]; then
	       cdo mergetime 6h0?.${varname}_${iyy}.nc tmp.${Fini}
               echo "cdo  expr,${oldname}=(${oldname} < 0.) ? 0. : ${oldname}"
               cdo expr,"${oldname}=( ${oldname} < 0. ) ? 0. : ${oldname}" tmp.${Fini} ${Fini}
	       echo "cdo shifttime,+6hours ${Fini} tmp.${Fini}"
               [ -s ${Fini} ] && rm tmp.${Fini}
               cdo shifttime,+6hours ${Fini} tmp.${Fini}
	       if [ -s tmp.${Fini} ] ; then
                  ncks -Oh -4 -L 4 tmp.${Fini} ${Fini} 
	          rm 6h0?.${varname}_${iyy}.nc 
                  rm tmp.${Fini}
	       fi
	    fi
         fi
      fi
   fi
   echo "last change for 6h.${varname} .."
   if [ -s ${Fini} -a -s ./6h ]; then
      # processing of changes for the varialbes name and add some attributes
      ncatted -Oh -a table,${oldname},a,s,${par0} ${Fini}
      ncatted -Oh -a code,${oldname},a,s,${par1} ${Fini}
      ncrename -Oh -v ${oldname},${varname} ${Fini}
      ncks -Oh -4 -L 4 ${Fini} ./6h/${Fini}
    #  Ftmp=tmp.${varname}_${iyy}.nc
    #  if [ $ii -gt 9 ]; then
    #     #ncap2 -Oh -s 'where(${varname}<0.) ${varname}=0;' ${Fini} ${Ftmp}
#	 #ncks -4 -L 4 ${Ftmp} ${Fini}
#	 #ncpdq ${Ftmp} ${Fini}
#         [ -s ${Ftmp} ] && rm ${Ftmp}
#      else
#	 ncks -Oh -4 -L 4 ${Fini}
#         [ -s ${Ftmp} ] && mv ${Ftmp}
#      fi
      if [ -s ./6h/${Fini} ]; then
         rm ${varname}_${iyy}.steps*.nc
         [ -s ${Fini} ] && rm ${Fini}
      fi
   fi
else
   echo " The file could be ready: ${Fini}"
   echo ""
fi

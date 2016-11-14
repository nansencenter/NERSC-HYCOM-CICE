#!/bin/bash

[ $# -ne 1 ] && { echo "usage: $(basename $0) <rungen>"  ; exit ; }
rungen=$1

# Clim_weekly program
prog=/home/fimm/nersc/knutali/Progs/Setup/Nest_Offline/clim_weekly


[ -f clim_weekly.in ] && rm clim_weekly.in 

# Prepare infile for fortran program clim_weekly
echo "Preparing infile for clim_weekly - rungen $rungen"
imonth=1
while [ $imonth -le 12 ] ; do
   iweek=1
   imonth2=$(echo 0$imonth | tail -c 3)

   while [ $iweek -le 4 ] ; do

      echo "# month=$imonth2 , week=$iweek" >> clim_weekly.in

      # "9999" is the climate field ;-) I was thinking 6666 but
      # its too obvious
      ls -1 ${rungen}AVE_[0-8]???_${imonth2}_${iweek}.a >> clim_weekly.in

      let iweek=$iweek+1
   done

   let imonth=$imonth+1
done

nice $prog

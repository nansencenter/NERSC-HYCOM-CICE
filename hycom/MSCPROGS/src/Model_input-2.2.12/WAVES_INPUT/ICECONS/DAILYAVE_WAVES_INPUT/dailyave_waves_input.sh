# #! /bin/bash

# dailyav_as_input.sh
# Author: Timothy Williams
# Date:   20130531, 15:45:32 CEST

# input directory;
in_dir=$1
# in_dir="/work/timill/TP4a0.12/expt_01.0/data"

# input directory;
# outdir="out"
outdir=.


for f in ${in_dir}/*DAILY_*.a
   # eg TP4DAILY_2010_169_2010_173.a
   do
      # input files
      nf=${#f}
      dn=2
      jj=$(($nf - $dn))
      g=${f:0:${jj}}.b
      # echo $f
      # echo $g
   
      # convert dump date to calendar date
      dn=10
      jj=$(($nf - $dn))
      dump_year=${f:$jj:4}
      # echo $dump_year

      dn=5
      jj=$(($nf - $dn))
      dump_day=${f:$jj:3}
      # echo $dump_day

      cdate=`~fanf/bin/jultodate $dump_day $dump_year 1 1`

   
      # output files
      # (just links to input files, but with names
      #  that are independent of the model region
      #  and the start date)
      fo=$outdir/dailyave_waves_input_${cdate}.a
      go=$outdir/dailyave_waves_input_${cdate}.b
      ln -sf $f $fo
      ln -sf $g $go
   done

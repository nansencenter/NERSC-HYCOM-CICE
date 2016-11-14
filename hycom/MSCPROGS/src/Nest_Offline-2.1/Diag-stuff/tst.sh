prog=/home/fimm/nersc/knutali/Progs/Setup/Nest_Offline/checknest

export LD_LIBRARY_PATH=/local/Matlab-R14sp3/bin/glnxa64/:$LD_LIBRARY_PATH


[ $# -ne 1 ] && exit
bound=$1

for i in nest_????_???_$bound ; do

   echo $i

   iyear=$( echo $i | sed "s/nest_//" | sed "s/_.*//")
   iday=$(  echo $i | sed "s/nest_...._//" | sed "s/_.*//")

#   echo $iyear $iday $bound | $prog
$prog << EOF
$iyear $iday
$bound
EOF


done

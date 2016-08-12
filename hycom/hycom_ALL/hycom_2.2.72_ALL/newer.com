#
# --- newer source file than the latest tar bundle
#
set echo
#
setenv N /u/data/wallcraf/TAR/hycom/ALL
ll $N
#
#cd ~/hycom/ALL
find . -name "*.[Ffch]" -newer $N -print

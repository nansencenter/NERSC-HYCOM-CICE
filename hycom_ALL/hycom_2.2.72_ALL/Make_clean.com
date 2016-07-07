#!/bin/csh
#
#set echo
#
# --- Usage:  ./Make_clean.com
#
# --- clean all HYCOM pre/post processing source directories
#
source Make_all.src
#
printenv ARCH
#
if (! -e ./config/${ARCH}_setup) then
  echo "ARCH = " $ARCH "  is not supported"
  exit 1
endif
#
setenv A $cwd
#
#foreach d ( archive force meanstd plot relax sample subregion topo)
foreach d ( archive force meanstd plot relax sample subregion topo cice ncom roms)
  echo "CLEANING ${d}/src:"
  cd ${A}/${d}/src
  make clean
end
#
# ./bin/Make_clean.com does not use Make_clean.src, but
#  should not need editing for any supported machine type.
#
echo "PROCESSING bin:"
cd ${A}/bin
csh Make_clean.com >&! Make_clean.log
cat Make_clean.log

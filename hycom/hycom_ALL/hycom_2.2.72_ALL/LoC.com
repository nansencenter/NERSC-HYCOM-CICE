#!/bin/csh
#
# --- count lines of code (simplest version)
#
setenv A $cwd
#
foreach d ( archive force meanstd plot relax sample subregion topo )
  cd ${A}/${d}/src
  cat *.F *.f *.c *.h >! ${A}/LoC_${d}.tmp
end
#
cd ${A}/bin
cat *.F *.f *.c *.h  >! ${A}/LoC_bin.tmp
#
cd $A
cat     LoC_*.tmp | wc -l
/bin/rm LoC_*.tmp

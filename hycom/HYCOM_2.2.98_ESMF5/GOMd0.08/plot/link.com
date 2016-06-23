#!/bin/csh
#
set echo
#
# --- plot directory softlinks for topography and source.
#
setenv R GOMd0.08
setenv N 02
#
touch   regional.grid.a
/bin/rm *.[ab]
#
/bin/ln -s ../topo/regional.grid.a   regional.grid.a
/bin/ln -s ../topo/regional.grid.b   regional.grid.b
#
/bin/ln -s ../topo/depth_${R}_${N}.a regional.depth.a
/bin/ln -s ../topo/depth_${R}_${N}.b regional.depth.b
#
if (-e ../topo/landsea_${R}.a) then
  /bin/ln -s ../topo/landsea_${R}.a regional.mask.a
endif
#
if (! -e src) then
  /bin/ln -s ../../ALL/plot/src_2.1.00 src
endif


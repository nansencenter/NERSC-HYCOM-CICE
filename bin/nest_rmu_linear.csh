#!/bin/csh
# Create rmu_nest files with variant e-folding time 
set echo
set time = 1
#

setenv TOOLS /home/sm_alfal/sea/TOPAZ5/NERSC-HYCOM-CICE/hycom/hycom_ALL/hycom_2.2.72_ALL
setenv TOOLS ~/HYCOM-tools
set path = ( ${path} ${TOOLS}/bin )
C
C --- Generate HYCOM nest rmu e-folding time, nest_rmu.a
C --- Script has idm and jdm values hardwired.
C --- S,E,N boundaries: 20 grid pts with 3-15 day e-folding time

C --- set copy command
if (-e ~${user}/bin/pget) then
C --- remote copy
  setenv pget ~${user}/bin/pget
  setenv pput ~${user}/bin/pput
else
C --- local copy
  setenv pget cp
  setenv pput cp
endif
setenv E `echo $cwd | awk '{print substr($0,length($0)-2,length($0))}'`
setenv X `echo ${E} | awk '{printf("%04.1f", $1*0.1)}'`

setenv DS /home/sm_alfal/sea/TOPAZ5/NERSC-HYCOM-CICE/TP5a0.06/

setenv R TP5a0.06
setenv T 05
C
C --- E is experiment number, from EXPT.src
C --- R is region identifier, from EXPT.src
C --- T is topog. identifier, from EXPT.src
C
C --- P is primary path,
C --- C is scratch directory,
C --- D is permanent directory,
C
#source EXPT.src
C
setenv P relax/${E}
setenv S ${DS}/${P}/
setenv IDM    `grep idm ${DS}/topo/regional.grid.b | awk '{print $1}'`
setenv JDM    `grep jdm ${DS}/topo/regional.grid.b | awk '{print $1}'`

mkdir -p $S
cd       $S

C
touch   rmu_linear fort.51 fort.51A regional.grid.a regional.grid.b
/bin/rm rmu_linear fort.51 fort.51A regional.grid.a regional.grid.b
C
C --- Input.
C
touch fort.51 fort.51A
if (-z fort.51) then
  ${pget} ${DS}/topo/depth_${R}_${T}.b fort.51 &
endif
if (-z fort.51A) then
  ${pget} ${DS}/topo/depth_${R}_${T}.a fort.51A &
endif
C
touch regional.grid.a regional.grid.b
if (-z regional.grid.b) then
  ${pget} ${DS}/topo/regional.grid.b regional.grid.b &
endif
if (-z regional.grid.a) then
  ${pget} ${DS}/topo/regional.grid.a regional.grid.a &
endif
C
touch rmu_linear
if (-z rmu_linear) then
  ${pget} ${TOOLS}/relax/src/rmu_linear . &
endif
wait
chmod a+rx rmu_linear
C
/bin/rm -f core 
touch core
/bin/rm -f fort.21 fort.21A
C
setenv FOR021A fort.21A
setenv FOR051A fort.51A
C
#234567890123456789012345678901234567890123456789012345678901234567890123456789
cat <<'E-o-D' >! fort.99
S,N boundaries: 20 grid pts with 3-15 day e-folding time
 345      'if    ' = first i point of sub-region (<=0 to end)
 694      'il    ' = last  i point of sub-region
   2      'jf    ' = first j point of sub-region
   5      'jl    ' = last  j point of sub-region
   3.0    'efoldA' = bottom left  e-folding time in days
   3.0    'efoldB' = bottom right e-folding time in days
   9.0    'efoldC' = top    right e-folding time in days
   9.0    'efoldD' = top    left  e-folding time in days
 345      'if    ' = first i point of sub-region (<=0 to end)
 694      'il    ' = last  i point of sub-region
   5      'jf    ' = first j point of sub-region
  11      'jl    ' = last  j point of sub-region
   9.0    'efoldA' = bottom left  e-folding time in days
   9.0    'efoldB' = bottom right e-folding time in days
  12.0    'efoldC' = top    right e-folding time in days
  12.0    'efoldD' = top    left  e-folding time in days
 345      'if    ' = first i point of sub-region (<=0 to end)
 694      'il    ' = last  i point of sub-region
  11      'jf    ' = first j point of sub-region
  20      'jl    ' = last  j point of sub-region
  12.0    'efoldA' = bottom left  e-folding time in days
  12.0    'efoldB' = bottom right e-folding time in days
  15.0    'efoldC' = top    right e-folding time in days
  15.0    'efoldD' = top    left  e-folding time in days
 17       'if    ' = first i point of sub-region (<=0 to end)
 206      'il    ' = last  i point of sub-region
 756      'jf    ' = first j point of sub-region
 759      'jl    ' = last  j point of sub-region
   9.0    'efoldA' = bottom left  e-folding time in days
   9.0    'efoldB' = bottom right e-folding time in days
   3.0    'efoldC' = top    right e-folding time in days
   3.0    'efoldD' = top    left  e-folding time in days
 17      'if    ' = first i point of sub-region (<=0 to end)
 206      'il    ' = last  i point of sub-region
 750       'jf    ' = first j point of sub-region
 756      'jl    ' = last  j point of sub-region
  12.0    'efoldA' = bottom left  e-folding time in days
  12.0    'efoldB' = bottom right e-folding time in days
   9.0    'efoldC' = top    right e-folding time in days
   9.0    'efoldD' = top    left  e-folding time in days
 17      'if    ' = first i point of sub-region (<=0 to end)
 206      'il    ' = last  i point of sub-region
 740       'jf    ' = first j point of sub-region
 750      'jl    ' = last  j point of sub-region
  15.0    'efoldA' = bottom left  e-folding time in days
  15.0    'efoldB' = bottom right e-folding time in days
  12.0    'efoldC' = top    right e-folding time in days
  12.0    'efoldD' = top    left  e-folding time in days
  -1      'if    ' = first i point of sub-region (<=0 to end)
'E-o-D'
./rmu_linear
C
C
C --- Output.
C
/bin/mv fort.21  ./nest_rmu.b
/bin/mv fort.21A ./nest_rmu.a
#
#rcp ./nest_rmu.b newton:${D}/nest_rmu.b
#rcp ./nest_rmu.a newton:${D}/nest_rmu.a
C
C --- generate invday for plotting
C
/bin/rm -f     nest_rmu_mask.a
hycom_bandmask nest_rmu.a 800 760 0.0 0.0 nest_rmu_mask.a >! nest_rmu_mask.b
cat            nest_rmu_mask.b
C
/bin/rm -f     nest_rmu_invday.a
hycom_expr     nest_rmu_mask.a INV 800 760 86400.0 1.0 nest_rmu_invday.a >! nest_rmu_invday.b
cat            nest_rmu_invday.b
C
hycom_zonal    nest_rmu_invday.a   800 760 1 >! nest_rmu_invday.zonal
C
C --- Delete all scratch directory files.
C
/bin/rm -f core fort* rmu_linear regional.grid.*

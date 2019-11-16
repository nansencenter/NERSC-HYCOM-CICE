#
set echo
#
# --- form the annual mean sst correction, ${S} * SSTme
# --- SSTme = model_annual_mean_SST - clim_annual_mean_SST
#
# --- Hardwired for GLBa0.72 (idm jdm = 800 760)
#
#cd /scr/${user}/hycom/GLBa0.72/force/offset
#
# --- S = degC to W/m^2; B = min allowed SSTme; T = max allowed SSTme;
# --- typically, T-B is 4 degC centered about basin-wide SSTme.
#
#setenv S -45.0
setenv S -60.0
#setenv S -60.0
#setenv B  -1.25
#setenv T   2.75
#
#setenv B  -0.25
#setenv T   6.75
#
#setenv B   0.0
#setenv T   4.0
#run8
#setenv B   0.0
#setenv T   3.0
#run10
#setenv B   0.0
#setenv T   2.75
#run11
setenv B   -0.5
setenv T   2.75
#
setenv E 010
setenv X expt_01.1
setenv Y 0005
#
touch   ${E}_sst_me.a offlux_${E}.a
/bin/rm ${E}_s*       offlux_${E}.?
#
#hycom_extract ../../${X}/data/meanstd/RSvs${E}_${Y}_sst.a 800 760 1 1 1 1 ${E}_sst_me.a >! ${E}_sst_me.b
hycom_extract TP5exp${E}_${Y}_sst.a 800 760 1 1 1 1 ${E}_sst_me.a >! ${E}_sst_me.b
#
hycom_clip   ${E}_sst_me.a    800 760 ${B} ${T} ${E}_sst_me_cl.a  >! ${E}_sst_me_cl.b
hycom_smooth ${E}_sst_me_cl.a 800 760  2        ${E}_sst_me_cls.a >! ${E}_sst_me_cls.b
#
echo "Heat Flux Offset (W/m^2)"           >! offlux_${E}.b
echo "${S}*(RS_SSTann-${X}_SSTann)"       >> offlux_${E}.b
echo "model mean from year ${Y}"          >> offlux_${E}.b
echo "anomaly between ${B} and ${T} degC" >> offlux_${E}.b
echo "i/jdm =  800  760"                  >> offlux_${E}.b
hycom_expr   ${E}_sst_me_cls.a ONE 800 760 ${S}  0.0 offlux_${E}.a | head -1 >> offlux_${E}.b
#
hycom_zonal_lat offlux_${E}.a 1 0.75 ./regional.grid.a >! offlux_${E}.txt

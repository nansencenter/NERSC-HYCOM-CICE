#
# --- create an update .tar bundle, based on newer.log
#
cd ~
tar cvf hycom_2.1.01_ALL_update1.tar \
hycom/ALL/plot/src_2.1.00/blkin.f \
hycom/ALL/plot/src_2.1.00/getdepth.f \
hycom/ALL/plot/src_2.1.00/xsecij.f \
hycom/ALL/plot/src_2.1.00/micomproc.f \
hycom/ALL/plot/src_2.1.00/hycomproc.f \
hycom/ALL/plot/src_2.1.00/hycomtest.f \
hycom/ALL/force/src_2.1.00/aphf_100_co.f \
hycom/ALL/force/src_2.1.00/aphf_1125_ec.f \
hycom/ALL/force/src_2.1.00/aphf_100_fn.f \
hycom/ALL/force/src_2.1.00/aphf_125_fn.f \
hycom/ALL/force/src_2.1.00/aphf_1875_nc.f \
hycom/ALL/force/src_2.1.00/aphf_250_ec.f \
hycom/ALL/force/src_2.1.00/aphf_100_soc.f \
hycom/ALL/subregion/src_2.1.00/sub_topog.f \
hycom/ALL/archive/src_2.1.00/horout_nc.f >&! hycom_2.1.01_ALL_update1.tar.lis


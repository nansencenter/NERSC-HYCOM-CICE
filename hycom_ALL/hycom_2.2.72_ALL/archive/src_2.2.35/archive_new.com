#
set echo
#
/bin/rm archive_new.tar
tar cvf archive_new.tar archv2data2d.f archv2data2t.f archv2data3z.f archv2restart.f horout.f horout_nc.f layer2z.f mod_restart.F restart2archv.f trim_archv.f >! archive_new.tar.lis

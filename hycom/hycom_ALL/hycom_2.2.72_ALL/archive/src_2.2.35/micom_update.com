#!/bin/csh
#
set echo
#
# --- Use *archv*.f as templates to update micom versions.
#
sed -e 's?lhycom/.true. /,?lhycom/.false./,?g' hycomarchv.f    >! micomarchv.f
sed -e 's?lhycom/.true. /,?lhycom/.false./,?g' archv2data2d.f  >! archm2data2d.f
sed -e 's?lhycom/.true. /,?lhycom/.false./,?g' archv2data3z.f  >! archm2data3z.f
sed -e 's?lhycom/.true. /,?lhycom/.false./,?g' archv2restart.f >! archm2restart.f

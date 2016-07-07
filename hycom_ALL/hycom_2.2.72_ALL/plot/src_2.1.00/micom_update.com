#
set echo
#
# --- Use hycom*.f as templates to update micom versions.
#
sed -e 's?lhycom/.true. /?lhycom/.false./?g' hycomproc.f  >! micomproc.f

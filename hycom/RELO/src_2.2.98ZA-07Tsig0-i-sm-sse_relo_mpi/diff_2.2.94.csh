#
set echo
#
#setenv C ~/hycom/RELO/src_2.2.94-17Tsig2-i-sm-sse_relo_mpi
setenv C ~/hycom/RELO_GLB/src_2.2.94-17Tsig2-i_relo_mpi
#
foreach f ( Makefile Make.csh )
  echo "*****     *****     *****     *****     *****     *****     *****"
  diff -bw $f $C
end
foreach f ( *.h *.c )
  echo "*****     *****     *****     *****     *****     *****     *****"
  diff -ibw $f $C
end
#allow for possible switch from .f to .F or .F to .f
foreach f ( *.f *.F )
  echo "*****     *****     *****     *****     *****     *****     *****"
  diff -ibw $f $C/$f:r.[Ff]
end

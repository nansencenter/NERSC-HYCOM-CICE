#
set echo
#
setenv C ~/hycom/RELO/src_2.2.97-17Tsig2-i-sm-sse_relo_one
#setenv C ~/hycom/RELO/src_2.2.97-17Tsig2-p_relo_one
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

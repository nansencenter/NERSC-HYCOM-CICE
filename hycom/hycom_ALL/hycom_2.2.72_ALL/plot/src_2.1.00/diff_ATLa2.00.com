#
set echo
#
foreach f ( README* Makefile *.h *.f *.F )
  echo "*****     *****     *****     *****     *****     *****     *****"
  diff -ibw $f ~/hycom/ATLa2.00/plot/src_2.0.03
end

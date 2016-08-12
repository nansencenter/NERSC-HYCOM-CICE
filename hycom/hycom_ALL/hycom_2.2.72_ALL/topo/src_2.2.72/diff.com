#
set echo
#
foreach f ( READ* Makefile *.f *.F )
  echo "*****     *****     *****     *****     *****     *****     *****"
  diff -ibw $f ~/hycom/ATLd0.32/topo/src_2.0.00
end

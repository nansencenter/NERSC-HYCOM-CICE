#
set echo
#
foreach f ( READ* Makefile *.h *.f *.F )
  echo "*****     *****     *****     *****     *****     *****     *****"
  diff -ibw $f ../../../ATLa2.00/topo/src_2.0.00
end

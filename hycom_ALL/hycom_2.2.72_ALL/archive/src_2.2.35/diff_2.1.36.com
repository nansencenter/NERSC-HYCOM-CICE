#
set echo
#
foreach f ( README* Makefile *.h *.f *.F )
  echo "*****     *****     *****     *****     *****     *****     *****"
  diff -ibw $f ../src_2.1.36
end

#
set echo
#
foreach f ( * )
  echo "*****     *****     *****     *****     *****     *****     *****"
  diff -ibw $f ../../ATLb2.00/config
end

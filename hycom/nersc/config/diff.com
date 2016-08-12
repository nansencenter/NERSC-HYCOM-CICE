#
set echo
#
foreach f ( RE* *_* )
  echo "*****     *****     *****     *****     *****     *****     *****"
  diff -ibw $f ../../ATLb2.00/config
end

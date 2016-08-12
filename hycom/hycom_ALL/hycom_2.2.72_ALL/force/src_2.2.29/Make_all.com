#
#set echo
#
# --- Usage:  ./Make_all.com >& Make_all.log
#
# --- make all force setup executables
#
# --- set ARCH to the correct value for this machine.
#
source Make_all.src
#
printenv ARCH
#
if (! -e ../../config/${ARCH}_setup) then
  echo "ARCH = " $ARCH "  is not supported"
  exit 1
endif
#
foreach m ( ap diff kp kpc kphfc pzero riv_mon riv_hf time_interp tp wi wi_curl wi_magstress wi_meanfit wc zero )
  make ${m} ARCH=${ARCH} >&! Make_${m}
  if ($status) then
    echo "Make failed:" ${m}
  else
    echo "Make worked:" ${m}
  endif
end
#
foreach m ( aphf_add aphf_climo aphf_diurnal aphf_extend aphf_flcorr aphf_meanfit aphf_monthly aphf_offset aphf_scale aphf_tacorr kphf_table )
  make ${m} ARCH=${ARCH} >&! Make_${m}
  if ($status) then
    echo "Make failed:" ${m}
  else
    echo "Make worked:" ${m}
  endif
end

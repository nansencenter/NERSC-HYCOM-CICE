#!/bin/csh
#set echo
#
# --- hardlink all Make_all.src files to the one in this directory.
#
foreach f ( */*/Make_all.src )
  /bin/rm -f $f
  /bin/ln Make_all.src $f
end

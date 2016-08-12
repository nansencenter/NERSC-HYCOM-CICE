#
#set echo
#
foreach f ( *.[Ffc] )
  foreach ff ( ../*/src_*/$f )
    /bin/rm -f $ff
    /bin/ln $f $ff
  end
end

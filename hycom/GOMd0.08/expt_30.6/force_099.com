#
set echo
date
#
# --- build all the forcing files for 1999
#
foreach y ( 099 )
# foreach m ( a b c d e f g h i j k l )
  foreach m ( i j k l )
    foreach t ( G L )
      awk -f 306.awk y01=${y} ab=${m} 306${t}.com >! 306${t}${y}${m}.com
      csh 306${t}${y}${m}.com >&! 306${t}${y}${m}.log &
    end 
    wait
    date
  end
end

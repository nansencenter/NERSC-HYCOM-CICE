      PROGRAM HYCOM_TIDAL_FOREMAN
      IMPLICIT NONE
C
C   hycom_tidal_foreman - Usage: hycom_tidal_foreman fssh.a f_day e_day i_hour constit [chunk [grid.a]] ftide.a
C
C                 Foreman tidal analysis of a SSH time series
C
C  fssh.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C   Each array represents a snapshot of sea surface height, and the
C   time between snapshots must be constant and a whole number of hours.
C
C  ftide.a will contain two fields for the real and imaginary parts of
C   each consituent in constit, in the order specified by constit.
C
C  f_day and e_day are the start and end of the time series, allowed
C   formats are: yyyy_ddd_hh, or yyyy_mm_dd_hh, or the hycom model day
C   as a floating point number.  A hycom model day will be rounded to
C   the nearest hour.
C
C  i_hour is the (integer) time increment of the time series in hours.
C   The time interval from f_day to e_day in hours must be an integer
C   multiple of i_hour.
C
C  constist is a string of tidal consistuent names separated by commas
C   with no spaces.  The order is arbitrary, but sets the order of
C   ftide.a.  Allowed values: m2,s2,n2,k2,k1,o1,p1,q1,mf,mm,m4,ms4,mn4.
C
C  chunk is the chunk size, it must be a power of 2 and the largest
C   allowed value of 4096 is also the default.  A smaller chunk size
C   requires less memory at the expense of more passes through fssh.a.
C
C  grid.a is a hycom grid file, default regional.grid.a.  Note that 
C   the corresponding grid.b must also exist.  Also, chunk must be
C   explicitly specified if grid.a is included on the command line.
C
C  Based on the "new" IOS Tidal Package by  Mike Foreman, see:
C  http://www.pac.dfo-mpo.gc.ca/science/oceans/tidal-marees/index-eng.htm
C
C  this version for "serial" Unix systems.
C
C  Danial R. Moore (QNA) and Alan J. Wallcraft (NRL), August 2012.
C

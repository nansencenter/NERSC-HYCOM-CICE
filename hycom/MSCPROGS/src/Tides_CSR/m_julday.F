      module m_julday
      contains
      subroutine julday (itimes,stime,utc)
c
c...purpose:  to convert between calendar day and julian date (utc).
c
c   output:
c       utc      double precision julian date
c
c   input:
c      itimes    integer array containing month, day, year, hour, minute
c      stime     floating point seconds
c
cGE      implicit double precision (a-h,o-z)
      dimension itimes(5)
c
      jd = itimes(2)-32075+1461*(itimes(3)+4800+(itimes(1)-14)/12)/4+367
     *   *(itimes(1)-2-(itimes(1)-14)/12*12)/12-3*((itimes(3)+4900+
     *   (itimes(1)-14)/12)/100)/4
c
      utc = dble(float(jd))-0.5d+0+dble(float(itimes(4)))/24.0d+0+
     *   dble(float(itimes(5)))/1440.0d+0+dble(stime)/86400.0d+0
c
      end subroutine
      end module

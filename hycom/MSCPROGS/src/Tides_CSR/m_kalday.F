      module m_kalday
      contains
      subroutine kalday (utc,mdyhms)
c
cGE      implicit double precision (a-h,o-z)
      dimension mdyhms(6),iptime(5)
c
c...compute the hour, minute, and second.
c
      t = utc
    5 jd = t
      xtime = mod((t-jd)*24.0d0+12.0d0,24.0d0)
      iptime(4) = xtime
      xtime = (xtime-iptime(4))*60.0d0
      iptime(5) = xtime
      ptimes = (xtime-iptime(5))*60.0d0
c
c...compute the year, month, and day.
c
      if (iptime(4).lt.12) jd = jd+1
      lx = jd+68569
      nx = 4*lx/146097
      lx = lx-(146097*nx+3)/4
      iptime(3) = 4000*(lx+1)/1461001
      lx = lx-1461*iptime(3)/4+31
      iptime(1) = 80*lx/2447
      iptime(2) = lx-2447*iptime(1)/80
      lx = iptime(1)/11
      iptime(1) = iptime(1)+2-12*lx
      iptime(3) = 100*(nx-49)+iptime(3)+lx
      do 10 i = 1,5
         mdyhms(i) = iptime(i)
   10 continue
      mdyhms(3) = mdyhms(3)-1900
      mdyhms(6) = ptimes+0.5d0
c
      if (mdyhms(4).eq.23 .and. mdyhms(5).eq.59 .and.
     .  mdyhms(6).eq.60) then
        t = t + 1.0d-6
        go to 5
      endif
c
      end subroutine
      end module

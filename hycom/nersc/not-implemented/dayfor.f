      logical function leapyear(iyr)
      implicit none
      integer iyr
      if ( (mod(iyr,4)==0 .and. mod(1901+iyr,400)==0).or. 
     &     mod(iyr,4)/=0 ) then
        leapyear=.false.
      else
        leapyear=.true.
      end if
      end function
c
C      subroutine forday(dtime,yrflag, iyear,iday,ihour)
C      implicit none
Cc
C      real*8  dtime
C      integer yrflag, iyear,iday,ihour
Cc
Cc --- converts model day to "calendar" date (year,ordinal-day,hour).
Cc
C      real*8  dtim1,day
C      integer iyr,nleap
Cc
C      if     (yrflag.eq.0) then
Cc ---   360 days per model year, starting Jan 16
C        iyear =  int((dtime+15.001d0)/360.d0) + 1
C        iday  =  mod( dtime+15.001d0 ,360.d0) + 1
C        ihour = (mod( dtime+15.001d0 ,360.d0) + 1.d0 - iday)*24.d0
Cc
C      elseif (yrflag.eq.1) then
Cc ---   366 days per model year, starting Jan 16
C        iyear =  int((dtime+15.001d0)/366.d0) + 1
C        iday  =  mod( dtime+15.001d0 ,366.d0) + 1
C        ihour = (mod( dtime+15.001d0 ,366.d0) + 1.d0 - iday)*24.d0
Cc
C      elseif (yrflag.eq.2) then
Cc ---   366 days per model year, starting Jan 01
C        iyear =  int((dtime+ 0.001d0)/366.d0) + 1
C        iday  =  mod( dtime+ 0.001d0 ,366.d0) + 1
C        ihour = (mod( dtime+ 0.001d0 ,366.d0) + 1.d0 - iday)*24.d0
Cc
C      elseif (yrflag.eq.3) then
Cc ---   model day is calendar days since 01/01/1901
C        iyr   = (dtime-1.d0)/365.25d0
C        nleap = iyr/4
C        dtim1 = 365.d0*iyr + nleap + 1.d0
C        day   = dtime - dtim1 + 1.d0
C        if     (dtim1.gt.dtime) then
C          iyr = iyr - 1
C        elseif (day.ge.367.d0) then
C          iyr = iyr + 1
C        elseif (day.ge.366.d0 .and. mod(iyr,4).ne.3) then
C          iyr = iyr + 1
C        endif
C        nleap = iyr/4
C        dtim1 = 365.d0*iyr + nleap + 1.d0
Cc
C        iyear =  1901 + iyr
C        iday  =  dtime - dtim1 + 1.001d0
C        ihour = (dtime - dtim1 + 1.001d0 - iday)*24.d0
Cc
C      else
C        write(6,*)
C        write(6,*) 'error in forday - unsupported yrflag value'
C        write(6,*)
C        stop '(forday)'
C      endif
C      return
C      end
Cc
CKAL  Basically does the opposite of forday
      subroutine dayfor(dtime,yrflag, iyear,iday,ihour)
      implicit none
c
      real*8  dtime
      integer yrflag, iyear,iday,ihour
c
c --- converts "calendar" date (year,ordinal-day,hour) to model day
c
      real*8  dtim1,day
      integer iyr,nleap
      integer k, diy, iday2
      logical, external :: leapyear
c
      if     (yrflag.eq.0) then
c ---   360 days per model year, starting Jan 16
        dtime =  (iyear-1) * 360.d0 +
     &           iday               +
     &           ihour/24.d0    -
     &           16.d0
c
      elseif (yrflag.eq.1) then
c ---   366 days per model year, starting Jan 16
        dtime =  (iyear-1) * 366.d0 +
     &           iday               +
     &           ihour/24.d0        -
     &           16.d0
c
      elseif (yrflag.eq.2) then
c ---   366 days per model year, starting Jan 01
        dtime =  (iyear-1) * 366.d0 +
     &           iday               -
     &           ihour/24.d0        -
     &           1.d0
c
      elseif (yrflag.eq.3) then
c ---   model day is calendar days since 01/01/1901
        iyr=iyear-1901
        dtime=0.
        do k=1,iyr
           if ( leapyear(k)) then
              dtime=dtime+366.d0
           else
              dtime=dtime+365.d0
           end if
        end do
c
        iday2=iday
        do while (iday2+ihour/24.>0.d0)
           iyr=iyr+1
           if ( leapyear(k)) then
              diy=366.d0
           else
              diy=365.d0
           end if
           if ((iday2+(ihour/24.))/diy>=1.d0) then
              dtime=dtime+diy
           else
              dtime=dtime+iday2+ihour/24.d0
           end if
           iday2=iday2-diy
        end do
      else
        write(6,*)
        write(6,*) 'error in forday - unsupported yrflag value'
        write(6,*)
        stop '(forday)'
      endif
      return
      end

program datetojul
use mod_year_info
implicit none 

!Which day was xxx, knowing that 1.1.1950 was a sunday.

#if defined(IARGC)
   integer*4, external :: iargc
#endif
   integer :: year, month,day,julday
   integer :: ryear, rmonth,rday
   character(len=10) :: tmpchar
   character(len=6) :: cjd
   integer :: fmt_

   
   if (iargc()==4) then ! 

      ! Actual days
      call getarg(1,tmpchar) ; read(tmpchar,'(i6)') julday

      ! Reference days
      call getarg(2,tmpchar) ; read(tmpchar,'(i4)') ryear
      call getarg(3,tmpchar) ; read(tmpchar,'(i2)') rmonth
      call getarg(4,tmpchar) ; read(tmpchar,'(i3)') rday

   else
      print '(a)','*****************************************************************'
      print '(a)','Program converts from julian day relative to reference date'
      print '(a)','and returns the corresponding  date in the format YYYYMMDD  where '
      print '(a)','YYYY=year, MM=month, DD=day of month'
      print '(a)',''
      print '(a)','Usage: jultodate  julianday refyear refmonth refday_of_month'
      print '(a)','Example - the command  \"jultodate 455 2008 1 1 \" returns 20090331'
      print '(a)','*****************************************************************'
      stop 
   end if

   fmt_=1

   ! Convert to date
   call juliantodate(julday,year,month,day,ryear,rmonth,rday)

   if (fmt_==1) then
      write(*,'(i4.4,i2.2,i2.2)') year,month,day
   end if



end program datetojul
   



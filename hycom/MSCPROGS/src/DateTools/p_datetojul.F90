program datetojul
use mod_year_info
implicit none 

#if defined(IARGC)
   integer*4, external :: iargc
#endif
   integer :: year, month,day,julday
   integer :: ryear, rmonth,rday
   character(len=10) :: tmpchar
   character(len=6) :: cjd

   
   if (iargc()==6) then ! 

      ! Actual days
      call getarg(1,tmpchar) ; read(tmpchar,'(i4)') year 
      call getarg(2,tmpchar) ; read(tmpchar,'(i2)') month
      call getarg(3,tmpchar) ; read(tmpchar,'(i3)') day

      ! Reference days
      call getarg(4,tmpchar) ; read(tmpchar,'(i4)') ryear 
      call getarg(5,tmpchar) ; read(tmpchar,'(i2)') rmonth
      call getarg(6,tmpchar) ; read(tmpchar,'(i3)') rday

   else
      print '(a)','*****************************************************************'
      print '(a)','Program converts from a specified date and returns julian day    '
      print '(a)','relative to the specified reference date'
      print '(a)',''
      print '(a)','Usage: datetojul year month day refyear refmonth refday'
      print '(a)','Example - the command  \"datetojul 2008 1 15 2008 1 1\" returns 14'
      print '(a)','*****************************************************************'
      stop
   end if

   ! Convert to Julian day
   julday =datetojulian(year,month,day,ryear,rmonth,rday)

   write(*,'(i5.5)') julday
end program datetojul
   



program datetojul
use mod_year_info
implicit none 


   integer*4, external :: iargc
   integer :: year, month,day,julday
   integer :: ryear, rmonth,rday
   character(len=10) :: tmpchar
   character(len=6) :: cjd

   
   if (iargc()==6) then ! 

      ! Actual days
      call getarg(1,tmpchar) ; read(tmpchar,'(i)') year 
      call getarg(2,tmpchar) ; read(tmpchar,'(i)') month
      call getarg(3,tmpchar) ; read(tmpchar,'(i)') day

      ! Reference days
      call getarg(4,tmpchar) ; read(tmpchar,'(i)') ryear 
      call getarg(5,tmpchar) ; read(tmpchar,'(i)') rmonth
      call getarg(6,tmpchar) ; read(tmpchar,'(i)') rday

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
   



module m_datetojulian
integer, dimension(12),parameter :: months_standard = &
   (/31,28,31,30,31,30,31,31,30,31,30,31/)
integer, dimension(12),parameter :: months_leapyear = &
   (/31,29,31,30,31,30,31,31,30,31,30,31/)
contains

integer function datetojulian(year,month,day,ryear,rmonth,rday)
   implicit none
   integer, intent(in) :: year,month,day,ryear,rmonth,rday
   integer :: iyear,sum_days,months(12)

   sum_days=0
   do iyear=ryear,year
      sum_days=sum_days+dayinyear(iyear)
      !print *,sum_days
   enddo

   ! Subtract from start of ref year to reference date
   months=monthsofyear(ryear)
   sum_days=sum_days          &
      - sum(months(1:rmonth)) &
      + months(rmonth) - rday + 1
   !print *,sum_days



   ! Subtract from end date in last year to end of year
   months=monthsofyear(year)
   sum_days=sum_days          &
      - sum(months(month:12)) &
      + day -1
   !print *,sum_days

   datetojulian=sum_days
end  function datetojulian





subroutine juliantodate(jday,year,month,day,ryear,rmonth,rday)
   implicit none
   integer, intent(in) :: jday,ryear,rmonth,rday
   integer, intent(out):: year,month,day
   integer :: iyear,sum_days,months(12),imonth,iday

   sum_days=0

   ! Subtract from start of ref year to reference date
   sum_days=sum_days+dayinyear(ryear)
   !print *,sum_days
   months=monthsofyear(ryear)
   sum_days=sum_days          &
      - sum(months(1:rmonth)) &
      + months(rmonth) - rday + 1
   !print *,sum_days


   ! Add years until beyond julian day
   iyear=ryear+1
   do while (sum_days<jday)
      sum_days=sum_days+dayinyear(iyear)
      iyear=iyear+1
      !print *,sum_days
   enddo
   if (sum_days>jday) then
      iyear=iyear-1
   end if

   imonth=12
   months=monthsofyear(iyear)
   do while (sum_days>jday)
      sum_days=sum_days-months(imonth)
      !print *,sum_days
      imonth=imonth-1
   enddo
   imonth=mod(imonth,12)+1

   iday=1
   do while (sum_days<jday)
      sum_days=sum_days+1
      iday=iday+1
   end do
      
   year=iyear
   month=imonth
   day=iday
end  subroutine juliantodate






integer function dayinyear(iyear)
   implicit none
   integer, intent(in) :: iyear
   !print *,iyear,mod(iyear,4)
   if (mod(iyear,4)==0 ) then
      if (mod(iyear,400)==0) then
         dayinyear=366
      else if (mod(iyear,100)==0) then
         dayinyear=365
      else
         dayinyear=366
      end if
   else
      dayinyear=365
   end if
end function dayinyear


function monthsofyear(iyear)
   implicit none
   integer :: monthsofyear(12)
   integer, intent(in) :: iyear
   if (dayinyear(iyear)==366) then
      monthsofyear=months_leapyear
   else
      monthsofyear=months_standard
   end if
end function

end module m_datetojulian

module mod_year_info
   type year_info
      integer iyy             ! year
      integer imm             ! month
      integer idd             ! day in year
      integer ihh             ! hours
      integer iss             ! seconds
      integer idm             ! current day in month
      character(len=4) cyy    ! year nr  (char)
      character(len=2) cmm    ! month nr (char)
      character(len=3) cdd    ! day in year  (char)
      character(len=2) chh    ! hour nr (char)
      character(len=4) css    ! second nr (char)
      character(len=3) month  ! month 'JAN' etc
      character(len=2) cdm    ! current day in month (char)
      integer totdim(12)      ! total Days In Months 
      integer daysinyear      ! total Days In year
   end type year_info

integer, dimension(12),parameter :: months_standard = &
   (/31,28,31,30,31,30,31,31,30,31,30,31/)
integer, dimension(12),parameter :: months_leapyear = &
   (/31,29,31,30,31,30,31,31,30,31,30,31/)

contains


subroutine year_day(time,refyear,tt,rforce)
   implicit none

! Calculates current day and year knowing that
! 1992 was skuddaar,  using model day (time) and 
! model refyear.

   real, intent(in)    :: time
   integer, intent(in) :: refyear
   character(len=5), intent(in) :: rforce
   type(year_info), intent(out) :: tt
   integer, parameter  :: extyear=1992

   integer days,days_in_year,i,itime
   real hr
   integer year,day,hour,second
   integer days_in_month(12)
   integer ttday





   itime=int(time)
   days=0
   do i=refyear,refyear+100000
      if (mod(i-extyear,4) == 0) then
         days_in_year=366
      else
         days_in_year=365
      endif

      if (rforce == 'month') days_in_year=360

      days=days+days_in_year
      if (itime < days) then
         year=i
         days=days-days_in_year
         day=itime-days
         hour=int( (time-float(itime))*24.0 )
         second=nint( ((time-float(itime))*24.0 - float(hour))*3600.0 )
         if (second == 3600) then
            hour=hour+1
            second=0
         endif
         if (hour == 24) then
            hour=0
            day=day+1
         endif
         if (day > days_in_year) then
            year=year+1
            day=0
         endif
         exit
      endif
   enddo
   tt%iyy=year
   tt%idd=day
   tt%ihh=hour
   tt%iss=second
   write(tt%cyy,'(i4.4)')tt%iyy
   write(tt%cdd,'(i3.3)')tt%idd
   write(tt%chh,'(i2.2)')tt%ihh
   write(tt%css,'(i4.4)')tt%iss

   if (days_in_year == 360) then
      do i=1,12
         days_in_month(i)=30
      enddo
   else
      days_in_month(1)=31
      if (days_in_year == 366) then
        days_in_month(2)=29
      else
        days_in_month(2)=28
      endif
      days_in_month(3)=31
      days_in_month(4)=30
      days_in_month(5)=31
      days_in_month(6)=30
      days_in_month(7)=31
      days_in_month(8)=31
      days_in_month(9)=30
      days_in_month(10)=31
      days_in_month(11)=30 
      days_in_month(12)=31
   endif
   tt%totdim=days_in_month

   ttday=tt%idd
   do i=1,12
      if (ttday >= days_in_month(i)) then
         ttday=ttday-days_in_month(i)
      else
         tt%imm=i
         write(tt%cmm,'(i2.2)')tt%imm
         tt%idm=ttday
         write(tt%cdm,'(i2.2)')ttday+1
         exit
      endif
   enddo
         
   select case (tt%imm)
   case (1)
      tt%month='JAN'
   case (2)
      tt%month='FEB'
   case (3)
      tt%month='MAR'
   case (4)
      tt%month='APR'
   case (5)
      tt%month='MAY'
   case (6)
      tt%month='JUN'
   case (7)
      tt%month='JUL'
   case (8)
      tt%month='AUG'
   case (9)
      tt%month='SEP'
   case (10)
      tt%month='OCT'
   case (11)
      tt%month='NOV'
   case (12)
      tt%month='DEC'
   end select

   tt%daysinyear=days_in_year
   

end subroutine year_day




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


end module mod_year_info

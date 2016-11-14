module m_year_day
contains
subroutine year_day(time,refyear,tt,rforce)
   use mod_year_info
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
         write(tt%cdm,'(i2.2)')ttday
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
end module m_year_day

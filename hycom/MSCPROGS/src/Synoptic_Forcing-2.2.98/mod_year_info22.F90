module mod_year_info22
implicit none
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
   integer, dimension(12),parameter :: months_360 = &
      (/30,30,30,30,30,30,30,30,30,30,30,30/)
   integer, dimension(12),parameter :: months_365 = &
      (/31,28,31,30,31,30,31,31,30,31,30,31/)
   integer, dimension(12),parameter :: months_366 = &
      (/31,28,31,30,31,30,31,31,30,31,30,31/)

contains



subroutine previousmonth(tt,year,month)
implicit none
type(year_info), intent(in) :: tt
integer, intent(out) :: month,year
month=mod(12+tt%imm-2,12)+1
year = tt%iyy - month/12
end subroutine


subroutine nextmonth(tt,year,month)
implicit none
type(year_info), intent(in) :: tt
integer, intent(out) :: month,year
month=mod(tt%imm,12)+1
year = tt%iyy + tt%imm/12
end subroutine



subroutine hycomday(dtime,iyear,imonth,idom,ihour,iss,tt,yrflag)
   implicit none
! Calculates hycom time based on a given _date_, also sets tt
   real*8,  intent(out) :: dtime
   integer, intent(in)  :: iyear, idom, imonth,ihour, iss
   type(year_info), intent(out) :: tt
   integer, intent(in)  :: yrflag

   integer days,days_in_year,k,idd
   integer months(12)
   real*8 dtime2

   if (yrflag==1 .or. yrflag==2) then
      print *,'hycomday not reliable for  yrflag=1,2'
      print *,'fix me'
      stop
   end if

   ! Get julian date of year based on month  and day of month
   ! NB - day of month starts from 1 (!)
   months=monthsofyear(iyear,yrflag)

   idd=0
   do k=1,imonth-1
      idd=idd+months(k)
   end do

   if (idom>months(imonth)) then
      print *,'impossible date : ',iyear, imonth, idom
      stop 
   else
      idd=idd+idom
   end if
   
   ! hycom time
   !print *,'calling dayfor in hycomday:',iyear,idd
   call dayfor(dtime,yrflag,iyear,idd,0)
   dtime=dtime+ihour/24.d0+iss/86400.d0
   !print *,'end calling dayfor in hycomday:',dtime

   ! set year_day
   call year_day(dtime,iyear,tt,yrflag)
end subroutine

   

subroutine year_day(dtime,refyear,tt,yrflag,refflag)
   implicit none
! Calculates current day and year knowing that
! 1992 was skuddaar,  using model day (dtime) and 
! model refyear.
   real*8, intent(in)    :: dtime
   integer, intent(in) :: refyear
   !character(len=5), intent(in) :: rforce
   type(year_info), intent(out) :: tt
   integer, intent(in) :: yrflag
   logical, optional   :: refflag ! if true, dtime is rel 1.1 of refyear
   integer days,days_in_year,i
   integer hour,second
   integer days_in_month(12)
   integer ttday,iyear, iday, ihour
   logical refflag2
   real*8 dtime2

   refflag2=.false.
   if (present(refflag)) refflag2=refflag

! Follow hycoms time convention
   if (refflag2) then
      call dayfor(dtime2,yrflag,refyear,1,0) ! NB - 1 jan is day 1
      dtime2=dtime+dtime2
   else
      dtime2=dtime
   end if
   call forday(dtime2,yrflag,iyear,iday,ihour)
   iday = iday - 1  ! hycom starts at 1 - we start at 0

   days_in_year =daysinyear  (iyear,yrflag)
   days_in_month=monthsofyear(iyear,yrflag)
   
   tt%totdim=days_in_month
   hour=int((dtime2-floor(dtime2))*24)
   second=nint(   ((dtime2-floor(dtime2))*24.0 - float(hour))   *3600.0 )
   tt%iyy=iyear
   tt%idd=iday
   tt%ihh=hour
   tt%iss=second
   write(tt%cyy,'(i4.4)')tt%iyy
   write(tt%cdd,'(i3.3)')tt%idd
   write(tt%chh,'(i2.2)')tt%ihh
   write(tt%css,'(i4.4)')tt%iss

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


!   integer function daysinyear(year,rforce)
!      implicit none
!      integer,intent(in) :: year
!      character(len=5),intent(in) :: rforce
!      ! If year is undivisible by 4, or if year is divisible by 100 and is
!      ! undivisible by 400 (e.g. 1900),  THEN year is not a leapyear
!      if (mod(year,4)/=0 .or. (mod(year,100)==0.and.mod(year,400)/=0)) then
!         daysinyear=365
!      else
!         daysinyear=366
!      end if
!      if (trim(rforce) == 'month') daysinyear = 360
!   end function daysinyear

   integer function daysinyear(year,yrflag)
   implicit none
   integer,intent(in) :: year
   integer,intent(in) :: yrflag
   ! If year is undivisible by 4, or if year is divisible by 100 and is
   ! undivisible by 400 (e.g. 1900),  THEN year is not a leapyear
   !if (mod(year,4)/=0 .or. (mod(year,100)==0.and.mod(year,400)/=0)) then
   !   daysinyear=365
   !else
   !   daysinyear=366
   !end if

   if (yrflag==0) then
      daysinyear=360
   elseif (yrflag==1) then
      daysinyear=365
   elseif (yrflag==2) then
      daysinyear=366
   elseif (yrflag==3) then
      if (mod(year,4)/=0 .or. (mod(year,100)==0.and.mod(year,400)/=0)) then
         daysinyear=365
      else
         daysinyear=366
      end if
   end if
   end function daysinyear




integer function datetojulian(year,month,day,ryear,rmonth,rday)
   implicit none
   integer, intent(in) :: year,month,day,ryear,rmonth,rday
   integer :: iyear,sum_days,months(12)

   sum_days=0
   do iyear=ryear,year
      !sum_days=sum_days+dayinyear(iyear)
      sum_days=sum_days+daysinyear(iyear,3)
      !print *,sum_days
   enddo

   ! Subtract from start of ref year to reference date
   months=monthsofyear(ryear,3)
   sum_days=sum_days          &
      - sum(months(1:rmonth)) &
      + months(rmonth) - rday + 1
   !print *,sum_days



   ! Subtract from end date in last year to end of year
   months=monthsofyear(year,3)
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
   !sum_days=sum_days+dayinyear(ryear)
   sum_days=sum_days+daysinyear(ryear,3)
   !print *,sum_days
   months=monthsofyear(ryear,3)
   sum_days=sum_days          &
      - sum(months(1:rmonth)) &
      + months(rmonth) - rday + 1
   !print *,sum_days


   ! Add years until beyond julian day
   iyear=ryear+1
   do while (sum_days<jday)
      sum_days=sum_days+daysinyear(iyear,3)
      iyear=iyear+1
      !print *,sum_days
   enddo
   if (sum_days>jday) then
      iyear=iyear-1
   end if

   imonth=12
   months=monthsofyear(iyear,3)
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


!integer function dayinyear(iyear)
!   implicit none
!   integer, intent(in) :: iyear
!   !print *,iyear,mod(iyear,4)
!   if (mod(iyear,4)==0 ) then
!      if (mod(iyear,400)==0) then
!         dayinyear=366
!      else if (mod(iyear,100)==0) then
!         dayinyear=365
!      else
!         dayinyear=366
!      end if
!   else
!      dayinyear=365
!   end if
!end function dayinyear


function monthsofyear(iyear,yrflag)
   implicit none
   integer :: monthsofyear(12)
   integer, intent(in) :: iyear,yrflag
   if (yrflag==0) then
      monthsofyear=months_360
   elseif (yrflag==1) then
      monthsofyear=months_365
   elseif (yrflag==2) then
      monthsofyear=months_366
   elseif (yrflag==3) then
      if (daysinyear(iyear,3)==366) then
         monthsofyear=months_leapyear
      else
         monthsofyear=months_standard
      end if
   end if
end function


end module mod_year_info22

module mod_diagnostics
implicit none


! Information regarding diagnostic/assimilation/restart times
   type diass
      integer nday                   ! Days for diagnostic or assimilation
      integer nhour                  ! Days for diagnostic or assimilation
      real*8 day                     ! Days for diagnostic or assimilation (hycom time format)
      logical ass                    ! true for assimilation time
      logical dia                    ! true for diagnosis time
      logical res                    ! true for saving restart file
   end type diass

!! Information about what to print at diagnostic times
!   type printing
!      logical cf                     ! true for saving central forecast
!      logical ice                    ! true for saving ice variables
!      logical eco                    ! true for saving ecosystem variables
!      logical ave                    ! true for saving ensemble mean
!      logical var                    ! true for saving ensemble variance
!   end type printing

   integer, parameter :: nrdday=1000
   integer,save       :: iday1, iday2
   type (diass),save  :: dday(0:nrdday) 

   logical, save :: loutput
   integer, save :: outday


contains


   logical function ioday(rt,dday,refyear,yrflag,baclin)
      use mod_xc
      use mod_year_info, only: year_info,daysinyear
      implicit none
      type (year_info),intent(in) :: rt
      type (diass),    intent(in) :: dday
      integer,         intent(in) :: refyear
      real,            intent(in) :: baclin
      integer         ,intent(in) :: yrflag 

      integer :: year, theday, theyear

      ! rt is days in year, dday is days from reference year. Fix it.
      theday  = dday%nday
      theyear = refyear
      do year = refyear,rt%iyy-1
         theday  = theday - daysinyear(year,yrflag)
         theyear = year
      end do

      ! Set status
      if (theday < 0 .or. theday > daysinyear(rt%iyy,yrflag) - 1 ) then
         ioday = .false.
      else

         ! For timesteps <= 3600.
         if (baclin <= 3600) then

            ioday = rt%idd==theday .and. rt%ihh==dday%nhour.and.rt%iss<baclin
         else if (baclin>3600 .and. baclin <= 86400 ) then

            ioday = rt%idd==dday%nday.and. (                                          &
                   ((dday%nhour - rt%ihh)*3600+rt%iss<baclin.and.dday%nhour-rt%ihh>=0 )   &
                   .or.                                                                   &
                   (rt%ihh<dday%nhour.and.rt%ihh*3600+rt%iss+baclin>86400 )               &
                   )
         else ! Not likely ...
            print *,'Timestep is not supported by this routine'
            call xcstop('(ioday)')
            stop '(ioday)'
         end if
      end if
      !print *,theyear,rt%iyy,theday,rt%idd,ioday
   end function ioday


   ! Set output flags and outday based on what is specified in
   ! ddays
   subroutine set_outputday(rt,refyear,yrflag,baclin,restart)
   use mod_year_info, only: year_info
   implicit none
   type(year_info), intent(in) :: rt
   integer, intent(in) :: refyear, yrflag
   real   , intent(in) :: baclin
   logical, intent(in) :: restart

   integer :: i

   ! Output days
   loutput=.false.
   do i=iday1,iday2+1
      if (ioday(rt,dday(i),refyear,yrflag,baclin)) then
         outday=i
         loutput=.true.
         exit
      endif
   enddo

   !Special case for restart
   if (restart) then
      outday=iday1
      loutput=.true.
      !Do not dump restart files initially
      dday(iday1)%res=.false.
   end if
   end subroutine


end module mod_diagnostics

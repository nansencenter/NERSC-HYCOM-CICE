module m_icestate_solar

implicit none

   !Thermodynamic parameters
   !REAL, PARAMETER :: absw1_cl = .6071    ! First water short wave absorption coef. for clear skies []
   !REAL, PARAMETER :: absw2_cl = .1119    ! Second water short wave absorption coef. for clear skies []
   !REAL, PARAMETER :: absw1_ov = .4212    ! First water short wave absorption coef. for overcast skies []
   !REAL, PARAMETER :: absw2_ov = .1292    ! Second water short wave absorption coef. for overcast skies []

   !REAL, PARAMETER :: s0       = 1365.    ! Solar constant [Wm^-2]
   !REAL, PARAMETER :: absh2o   = .09      ! Absorption of water and ozone, -

   real, parameter :: snwlim_trans=1e-4   ! Limiting snow thickness retaining sw at surface

   !integer, parameter :: ifrac=24             ! Hours in a day 


contains




! -----------------------------------------------------------------------
! ------------------------- FUNCTION SUN_RADFL --------------------------
! -----------------------------------------------------------------------
! Calculates day averaged incident solar radiative flux at the marine 
! interface, as a function of latitude and timeref (number of seconds
! since the beginning of the experiment, on jan. 1st, 00:00:00 TU.  
!subroutine sun_radfl(cawdir,radfl,ccx,rlatx,timerefx,rforce) !,ldebug)
!  use mod_icestate_thermo
!  implicit none
!  REAL, intent(out):: cawdir,radfl
!  REAL, INTENT(in) :: ccx,rlatx
!  REAL*8, INTENT(in) :: timerefx
!  character*5, intent(in) :: rforce
!  real dangle, decli, sundv, sin2, cos2, scosz, stot, &
!       bioday, biohr, hangle, cosz, srad, sdir, sdif, altdeg, cfac, &
!       ssurf
!  integer npart
!  !logical ldebug
!  real*8 daymodulus,day
!
!  if (rforce=='month') then
!     daymodulus=360.     ! Yearflag=0, each year has 360 days
!  else
!     daymodulus=365.2425 ! Average number of days in a gregorian calendar year
!  end if
!
!
!  !dtime  = timerefx
!  !day    = aint(dtime)                      ! accumulated day number 
!  !day    = modulo(day,365.)                 ! 0 < day < 364
!  !dangle = pi2*day/365                     ! day-number-angle, in rad.
!
!
!  day=dmod(timerefx,daymodulus)
!  dangle = pi2*day/daymodulus                     ! day-number-angle, in rad.
!
!  !if(ldebug) print *,'KAL solar 1',day,timerefx,daymodulus
!
!!
!! compute astronomic quantities
!  decli = .006918+.070257*sin(dangle)   -.399912*cos(dangle)   &
!                 +.000907*sin(2.*dangle)-.006758*cos(2.*dangle)  &
!                 +.001480*sin(3.*dangle)-.002697*cos(3.*dangle)
!
!  sundv = 1.00011+.001280*sin(dangle)   +.034221*cos(dangle)   &
!                 +.000077*sin(2.*dangle)+.000719*cos(2.*dangle)
!
!  sin2 = sin(rlatx)*sin(decli)
!  cos2 = cos(rlatx)*cos(decli)
!
!! split each day into ifrac parts, and compute the solar radiance for
!! each part. By assuming symmetry of the irradiance about noon, it
!! is sufficient to compute the irradiance for the first 12 hrs of
!! the (24 hrs) day (mean for the first 12 hrs equals then the mean
!! for the last 12 hrs)
!  scosz = 0.
!  stot  = 0.
!
!  !if (ldebug)  print *,'KAL solar 2',decli,sundv,sin2,cos2,day,rlatx,decli
!
!  DO npart = 1,ifrac
!    bioday = day+0.5*(npart-.5)/FLOAT(ifrac)
!    biohr  = bioday*86400.                   ! hour of day in seconds
!    biohr  = modulo(biohr+43200.,86400.)     ! hour of day;  biohr=0  at noon
!
!    hangle = pi2*biohr/86400.                ! hour angle, in radians
!    cosz   = max(0.,sin2+cos2*cos(hangle))   ! cosine of the zenith angle
!    scosz  = scosz+cosz                      !  ..accumulated..
!
!    srad = s0*sundv*cosz                     ! extraterrestrial radiation
!
!    !KAL obs: .7^100 = 3x10^-16 , the min function avoids underflow
!    !KAL obs: .7^40  = 6x10^-7 , the min function avoids underflow
!    !sdir = srad*0.7**(1./(cosz+epsil1))      ! direct radiation component
!    !sdir=srad*0.7**(min(100.,1./(cosz+epsil1)))    !direct radiation component
!    sdir=srad*0.7**(min(40.,1./(cosz+epsil1)))    !direct radiation component
!
!    sdif = ((1.-absh2o)*srad-sdir)*.5        ! diffusive radiation component
!
!    ! solar noon altitude in degrees
!    !altdeg = max(0.,asin(sin2+cos2))*360./pi2
!    ! KAL - keep asin argument within allowed bounds - especially important when
!    ! default real is real*4.....
!    altdeg = max(0.,asin(max(-1.,min(sin2+cos2,1.))))*360./pi2
!
!
!    cfac  = (1.- 0.62 * ccx + 0.0019*altdeg) ! cloud correction (Reed, 77)
!    !cfac  = 1.
!    ssurf = (sdir+sdif)*cfac
!
!    stot = stot+ssurf
!
!     !if (ldebug)  print *,'KAL solar 3',npart,cfac,altdeg,sin2+cos2
!
!  END DO 
!
!  scosz=scosz/float(ifrac)                      !24-hrs mean of  cosz
!  radfl = stot/FLOAT(ifrac)                     !24-hrs mean shortw rad in w/m^2
!
!  cawdir    = 1.-amax1(0.15,0.05/(scosz+0.15))  ! Co-albedo over water for dir. light
!  !if (ldebug)  print *,'KAL solar 4',scosz,radfl,cawdir
!
!
!END subroutine sun_radfl
!










! -----------------------------------------------------------------------
! ------------------------ SUBROUTINE CALC_TRANS ------------------------
! -----------------------------------------------------------------------
! If qswx is the solar flux at the ice-air interface, only part
! of it, (1 - I0) * qswx is absorbed by the top of the ice slab.    
!
!   - If the ice is snow covered, or if the heat reservoir Qstore is 
! full, then I0 = 0. Else :
!   - If hice > 0.1 m, then I0 = 0.17 
!   - If hice < 0.1 m, then I0 = 1 - 0.83 * hice / 0.1 
SUBROUTINE calc_trans(icem,inodx,qmaxx,transx)
  USE mod_icestate
  IMPLICIT NONE
  type(t_ice), intent(inout) :: icem(nthick)
  REAL, DIMENSION(nthick), INTENT(out)   :: inodx,qmaxx,transx
  integer hk

  qmaxx   = 0.
  inodx   = 0.
  transx  = 0.

  !qmaxx = qst_frac * hofusn0 * icem%hice

  ! Restrict to top 1m
  qmaxx = qst_frac * hofusn0 * min(icem%hice,1.)

  where ( (icem%hsnw<snwlim_trans.and.(qmaxx-icem%qstore)>epsil1.and.icem%fice>epsil1) )
     transx = exp(-1.5 * (icem%hice - .1))
     where (icem%hice>.1) 
       inodx  = .17
     elsewhere
       inodx  = 1. - .83 * icem%hice / .1
     endwhere
   endwhere

END SUBROUTINE calc_trans


end module m_icestate_solar



   

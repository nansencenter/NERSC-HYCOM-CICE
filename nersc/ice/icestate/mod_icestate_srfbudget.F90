module mod_icestate_srfbudget
   implicit none


! ====================================================================
! Routines for calculating the surface energy budget in the icestate
! model. This is done by expressing the budget as 4th order polynomials.
! This might seem awkward if we consider calculating heat fluxes, but is
! useful when we want to solve for surface temperature of the ice
! ====================================================================






contains





! ====================================================================
! Routine for calculating ice surface temperature using atmospheric
! heat fluxes and conductive and thermal inertia heat fluxes.
! The routine uses function calls for different flux components
! and sets the heat balance up as a function of surface temperature. 
! The routine sets  up and solves:
!
! a(0) + a(1)*T^1 + a(2)*T^2 + a(3)*T^3 + a(4)*T^4 =0 
!
! For now we onlys olve a linear equation
! --------------------------------------------------------------------
! Note: In this manner it is easy to iterate and/or use Newton solver
! ====================================================================
   elemental function srftemp_ice(icem,t_old,relhum,tair,wind,qsw,slp,clouds,rlat,tml)
      use mod_icestate
      implicit none
      type(t_ice), intent(in) :: icem
      real,        intent(in) :: relhum, tair, wind, qsw, slp, clouds, rlat, tml,t_old

      ! a(n) is coefficient of T^n. Equation to be solved is  sum(a(n)T^n) = 0
      real, dimension(0:4) :: a
      real :: rhoair,vpair,vpsrf,qair,qsrf
      real :: srftemp_ice

      ! Density, vapour pressure and humidities.
      rhoair  = slp / (gasconst * tair)
      vpair   = vapp(aice  ,bice  ,tair)
      vpsrf   = vapp(aice,bice,icem%tsrf)
      qair    = relhum*humid(slp,vpair)
      qsrf    = humid(slp,vpsrf)
 
      ! Components of the surface budget
      a=0
      a(0) = qsw
      a = a + turbulent(icem%tsrf,qsrf,tair,qair,rhoair,wind)
      a = a + longwave(tair,clouds,vpair,rlat)
      a = a + heatconduction(icem,tml)
      a = a + heatinertia(icem,t_old)

      if ( all(abs(a(3:4))<epsil1) ) then
         srftemp_ice= - a(0) / a(1)
      else
         srftemp_ice=-999
      end if
   end function









! ====================================================================
! Routine for calculating Net outgoing atmospheric heat flux.
! The routine uses function calls for different flux components
! and sets the heat balance up as a function of surface temperature. 
! The final flux is calculated as:
! 
! atm_heatflux=a(0) + a(1)*T^1 + a(2)*T^2 + a(3)*T^3 + a(4)*T^4 
! ====================================================================
   elemental function atm_heatflux(tsrf,relhum,tair,wind,qsw,slp,clouds,rlat,ice)
      use mod_icestate
      implicit none
       real,        intent(in) :: tsrf, relhum, tair, wind, qsw, slp, clouds, rlat
       logical,     intent(in) :: ice
       real :: atm_heatflux
       real :: a(0:4)

       real :: aa,bb,vpair,qair,vpsrf,qsrf,rhoair

  
      if (ice) then
         aa = aice
         bb = bice
      else
         aa = awater
         bb = bwater
      end if

      ! Density of air
      rhoair   = slp / (gasconst * tair) 

      ! Vapour pressure and humidity of air over the surface
      vpair   = vapp(aa  ,bb  ,tair)  
      qair    = relhum*humid(slp,vpair)

      ! Vapour pressure and humidity of the air over the surface
      vpsrf     = vapp(aa,bb,tsrf)
      qsrf      = humid(slp,vpsrf)
 
      !       

      a=0     
      ! Shortwave component
      a(0) = qsw

      ! Sensible and latent heat transfer
      a = a + turbulent(tsrf,qsrf,tair,qair,rhoair,wind)

      ! Longwave parametrization
      a = a + longwave(tair,clouds,vpair,rlat)

      ! Heat flux
      atm_heatflux = a(0) + tsrf*(a(1) + tsrf*(a(2) + tsrf*(a(3) + tsrf*a(4))))

   end function atm_heatflux


! ====================================================================
! Routine for calculating conductive heat flux.
! The routine uses function calls for different flux components
! and sets the heat balance up as a function of surface temperature. 
! The final flux is calculated as:
! 
! cond_heatflux=a(0) + a(1)*T^1 + a(2)*T^2 + a(3)*T^3 + a(4)*T^4 
! ====================================================================
   elemental function cond_heatflux(icem,tml)
      use mod_icestate
      implicit none
      type(t_ice), intent(in) :: icem
      real,        intent(in) :: tml
      real :: cond_heatflux
      real :: a(0:4)

      !       
      a=0     
      a = a + heatconduction(icem,tml)

      ! Heat flux
      cond_heatflux = a(0) + icem%tsrf*(a(1) + icem%tsrf*(a(2) + icem%tsrf*(a(3) + icem%tsrf*a(4))))
   end function cond_heatflux



! ====================================================================
! Routine for calculating thermal inertia heat flux.
! The routine uses function calls for different flux components
! and sets the heat balance up as a function of surface temperature. 
! The final flux is calculated as:
! 
! inertia_heatflux=a(0) + a(1)*T^1 + a(2)*T^2 + a(3)*T^3 + a(4)*T^4 
! ====================================================================
   elemental function inertia_heatflux(icem,tsrf_old)
      use mod_icestate
      implicit none
      type(t_ice), intent(in) :: icem
      real,        intent(in) :: tsrf_old
      real :: inertia_heatflux
      real :: a(0:4)

      !       
      a=0     
      a = a + heatinertia(icem,tsrf_old)

      ! Heat flux
      inertia_heatflux = a(0) + icem%tsrf*(a(1) + icem%tsrf*(a(2) + icem%tsrf*(a(3) + icem%tsrf*a(4))))
   end function inertia_heatflux





! ====================================================================
!   Turbulent heat flux routine . Calculate the coefficients 
!   for a polynomial in T (surface temperature) which describe 
!   conductive heat flux.
! ====================================================================
   pure function turbulent(tsrf,qsrf,tair,qair,rhoair,wind)
      use mod_icestate
      implicit none
      real, intent(in) :: wind,tair,tsrf,qair,rhoair,qsrf
      real, dimension(0:4) :: turbulent
      real :: clatent,csens,fqlat,fqsens,dtemp
      integer :: iturb,jturb


      ! Sensible and latent heat coeffs in solving for T
      iturb    = min(16,int(wind * .5) + 1)
      dtemp    = tair - tsrf
      jturb    = max(1,min(29,int((dtemp + 7.) *2.) +1))
      clatent = clat(iturb,jturb)
      csens   = .94 * clatent
      fqlat   = rhoair * hosubl * clatent * wind
      fqsens  = rhoair * cpair  * csens   * wind

      ! Coefficients 
      turbulent(0:4)=0.
      turbulent(0) = fqsens * tair - fqlat * max(0.,qsrf-qair)
      turbulent(1) = -fqsens
   end function turbulent



   ! Wrapper for flux calc -- diag
   elemental function turb_flux(tair,wind,clouds,tsrf,slp,relhum,ice)
      use mod_icestate
      implicit none
      real, intent(in) :: tair, wind, tsrf,clouds,slp,relhum
      logical, intent(in) :: ice
      real :: turb_flux,aa,bb,a(0:4),qair,vpair,qsrf,vpsrf,rhoair

      if (ice) then
         aa = aice
         bb = bice
      else
         aa = awater
         bb = bwater
      end if

      ! Density of air -- etc ...
      rhoair   = slp / (gasconst * tair) 
      vpair   = vapp(aa  ,bb  ,tair)  
      vpsrf     = vapp(aa,bb,tsrf)
      qair    = relhum*humid(slp,vpair)
      qsrf      = humid(slp,vpsrf)
 
      a = turbulent(tsrf,qsrf,tair,qair,rhoair,wind)
      turb_flux=a(0)+tsrf*(a(1)+tsrf*(a(2)+tsrf*(a(3)+tsrf*(a(4)))))
   end function turb_flux



      
      
! ====================================================================
!   Longwave heat flux routine . Calculate the coefficients 
!   for a polynomial in T (surface temperature) which describe 
!   conductive heat flux.
! ====================================================================
   pure function longwave(tair,clouds,vpair,rlat)
      use mod_icestate
      implicit none
      real,        intent(in) :: tair, clouds, vpair, rlat
      real :: fqlw, fqlwcc,fqlw1,fqlw2
      real :: longwave(0:4)

      ! --- upward net long wave radiation flux factor (w / (m^2 k))
      !  .    emiss       -                   emissivity
      !  .    stefanb     w / (m^2 k^4)       stefan-boltzman constant
      !  .    rlat        rad                 latitude
      fqlw  =emiss*stefanb*tair**3
      fqlwcc=1.-(.5+.246*abs(rlat))*clouds**1.2
      fqlw1=fqlw*tair*((.254-4.95e-5*vpair)*fqlwcc-4.)
      fqlw2=fqlw*4.

      ! This paramtetrization is linearly dep. upon temperature
      longwave=0.
      longwave(0) = -fqlw1
      longwave(1) = -fqlw2
   end function longwave


   ! Wrapper for longwave -- diag
   elemental function longwave_flux(tair,clouds,rlat,tsrf,ice)
      use mod_icestate
      implicit none
      real,        intent(in) :: tair, clouds, rlat, tsrf
      logical,     intent(in) :: ice
      real :: fqlw, fqlwcc,fqlw1,fqlw2,a(0:4),aa,bb,vpair
      real :: longwave_flux


      if (ice) then
         aa = aice
         bb = bice
      else
         aa = awater
         bb = bwater
      end if
      vpair   = vapp(aa  ,bb  ,tair)  
      a=longwave(tair,clouds,vpair,rlat)

      longwave_flux=a(0)+tsrf*(a(1)+tsrf*(a(2)+tsrf*(a(3)+tsrf*(a(4)))))
      end function longwave_flux





! ====================================================================
!   Conductive heat flux routine . Calculate the coefficients 
!   for a polynomial in T (surface temperature) which describe 
!   conductive heat flux.
! ====================================================================
   pure function heatconduction(icem,tml)
      use mod_icestate
#if defined(SSNOWD) && defined (HEAT_CONDUC) 
      use mod_icestate_tools
#endif      
      implicit none
      type(t_ice), intent(in) :: icem
      real,        intent(in) :: tml
      real :: rksnw,rkeq,tbot
      real :: heatconduction(0:4)
         
      ! Snow conductivity
      rksnw   = rkice*(icem%rhosnw / rhow)**1.885

      ! Effective heat conductivity factor
      rkeq=0.
      if (icem%hice>epsil1) then
      if (icem%nlay>1) then
#if defined(SSNOWD) && defined (HEAT_CONDUC)
         !Effects of snow cover heterogneities are included via the conduction
         !  correction factor
         icem%hsnw = average_depth(icem)
         rkeq    = corr_factor(icem)*rkice*rksnw / (rkice*icem%hsnw + rksnw*icem%hice/(2*icem%nlay))
#else         
         rkeq    = rkice*rksnw / (rkice*icem%hsnw + rksnw*icem%hice/(2*icem%nlay))
#endif         
         tbot    = icem%vtp(icem%nlay)
      else
#if defined(SSNOWD) && defined (HEAT_CONDUC) 
         icem%hsnw = average_depth(icem)
         rkeq    = corr_factor(icem) *rkice*rksnw / (rkice*icem%hsnw + rksnw*icem%hice)
#else         
         rkeq    = rkice*rksnw / (rkice*icem%hsnw + rksnw*icem%hice)
#endif         
         tbot    = tml
      end if
      end if

      ! Linear dependence...
      heatconduction    = 0.
      heatconduction(0) = rkeq*tbot
      heatconduction(1) = -rkeq
   end function heatconduction



! ====================================================================
!   Thermal inertia heat flux routine . Calculate the coefficients 
!   for a polynomial in T (surface temperature) which describe 
!   conductive heat flux.
! ====================================================================
   pure function heatinertia(icem,t_old)
      use mod_icestate
      implicit none
      type(t_ice), intent(in) :: icem
      real       , intent(in) :: t_old
      real :: fqinert,rkeq,rksnw
      real :: heatinertia(0:4)


      ! Snow conductivity
      rksnw   = rkice*(icem%rhosnw / rhow)**1.885

      ! Calculate heat conductivity  and heat inertia factors (cond. factor 
      ! is used because we treat snow with no heat capacity, we assume linear
      ! effective vtp in the snow and ice)
      fqinert=0.
      if (icem%fice>epsil1) then
      if (icem%nlay>1) then
         rkeq    = rkice*rksnw / (rkice*icem%hsnw + rksnw*icem%hice/(2*icem%nlay))
         fqinert = .5*(icem%hice/(2*icem%nlay))*rhoice*cpice*rkeq / (dtt*rkice/(icem%hice/(2*icem%nlay)))
      else
         rkeq    = rkice*rksnw / (rkice*icem%hsnw + rksnw*icem%hice)
         fqinert = .5*icem%hice*rhoice*cpice*rkeq / (dtt*rkice/icem%hice)
      end if
      end if

      ! Linear dependence..
      heatinertia=0. 
      heatinertia(0) = fqinert*t_old   ! Old surf temp for inertia ...
      heatinertia(1) = -fqinert
   end function heatinertia



! ====================================================================
!   Routine to get evaporation over water. 
! ====================================================================
   pure function evaporation(tair,tsrf,wind,slp,relhum,rhosw)
      use mod_icestate
      implicit none
      real, intent(in) :: wind,tair,tsrf,slp,rhosw,relhum
      real :: evaporation
      real :: clatent,fqlat,dtemp,rhoair,vpair,vpsrf,qair,qsrf
      integer :: iturb,jturb


      rhoair   = slp / (gasconst * tair) 

      ! Vapour pressure and humidity of air over the surface
      vpair   = vapp(awater  ,bwater  ,tair)  
      qair    = relhum*humid(slp,vpair)

      ! Vapour pressure and humidity of the air over the surface
      vpsrf     = vapp(awater,bwater,tsrf)
      qsrf      = humid(slp,vpsrf)
 
      dtemp    = tair-tsrf
      iturb    = min(16,int(wind * .5) + 1)
      jturb    = max(1,min(29,int((dtemp + 7.) *2.) +1))
      clatent  = clat(iturb,jturb)
      evaporation = rhoair * hosubl * clatent * wind * max(0.,qsrf-qair) * hocondi / rhosw

   end function evaporation





! =======================================================================
! Defines vapour pressure function (vapp has unit pa; tx is temp in k). 
!
! =======================================================================
   elemental FUNCTION vapp(ax,bx,tx)
     implicit none
     REAL ::  vapp
     REAL, INTENT(in) ::   ax,bx,tx
     vapp = 611. * 10.**(ax * (tx - 273.16)/(tx - bx))
   END FUNCTION vapp




! =======================================================================
! Defines specific humidity (dimensionless).
!
! =======================================================================
   elemental FUNCTION humid(px,ex)
     implicit none
     REAL ::              humid
     REAL, intent(in) ::  px,ex
     humid = .622 * ex / (px - .378 * ex)
   END FUNCTION humid


! =======================================================================
! =========================== clat_turbm ================================
! =======================================================================
!
!  Read numerical values of the latent transfer coefficient based
!  on the tabelled values in isemer et al, j clim., p. 1180, 1989
!  note:
!  i-index gives the wind at 10m height
!             from 0 to 30 m/s in intervals of 2 m/s
!  j-index gives the virtual air-sea temperature difference at 10m 
!  height
!             from -7 to +7 deg c in intervals of .5 deg c
!  for all but the equatorial and sub-equatorial waters, virtual 
!  temp is close to real temp (see gill, 1982, p. 41)
!
!  PURPOSE  
!  -------
!  Computes turbulent heat exchange coefficient as a function of wind
!  and the temperature difference between marine surface and 
!  atmospheric mixed layer. 
!
!  EXTERNALS 
!  ---------
!  Called by procedures : th/thermf 
!  Calls procedures     : none
!
!  AUTHOR
!  ------
!  Paul Budgell, 1994
!
!  MODIFIED
!  --------
!  David Salas y Melia, nov. 1996
!

   subroutine clat_turbm(clatx)
      use mod_icestate
      implicit none
      REAL, DIMENSION(16,29), INTENT(out) :: clatx
      CHARACTER(80)                       :: filename
      REAL                                :: d_wind,d_temp
      REAL, DIMENSION(16,29)              :: clat1
      integer :: i,j, iossum, ios
      logical :: ex

      ! READ data file
      filename = 'iwh_tabulated.dat'
      inquire(exist=ex,file=trim(filename))
      if (.not. ex) then
         if (mnproc==1) print *,'Can not find '//trim(filename)
         call xcstop('(mod_icestate_srfbdget:clat_turbm)')
      end if


      iossum=0
      OPEN(10,file = TRIM(filename))
      DO j = 1,29       
         DO i = 1,16
            READ(10,*,iostat=ios) d_wind,d_temp,clat1(i,j)
            clat1(i,j) = clat1(i,j) * 1.E-3
            iossum=iossum+ios
         END DO
      END DO
      CLOSE(10)
      if (iossum/=0) then
         if (mnproc==1) print *,'Error reading '//trim(filename)
         call xcstop('(mod_icestate_srfbdget:clat_turbm)')
      end if

      ! let  i = 1  represent wind speeds in the interval [0, 2) m/s
      !      i = 2  represent wind speeds in the interval [2, 4) m/s ...
      DO j = 1,29
         DO i = 1,15
            clatx(i,j) = .5 * (clat1(i,j) + clat1(i+1,j))
         END DO
         clatx(16,j) = clatx(15,j)
      END DO

      ! let  j = 1  represent temp differences in the interval [-7, -6.5) m/s
      !      j = 2  represent temp differences in the interval [-6.5, -6) m/s ...
      DO i = 1,16
         DO j = 1,28
            clatx(i,j) = .5 * (clat1(i,j) + clat1(i,j+1))
         END DO
         clatx(i,29) = clatx(i,28)
      END DO

  end subroutine clat_turbm



end module mod_icestate_srfbudget

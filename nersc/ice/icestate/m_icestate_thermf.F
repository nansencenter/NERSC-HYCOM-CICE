      module m_icestate_thermf

!NB!! vertical ice growth treatment and lead treatment should be separated for
! efficiency of the code

! =======================================================================================
! This version of the icestate thermodynamics drops all lead parametrizations from the 
! original code. In stead the parametrization from Drange and Simonsen (1996, MICOM thermf
! at the Nansen Center) is used.
! 
! In addition another implementation of the ice thermodynamics is used.
! Instead of solving the heat equation using the method of M&C (1973), it is here solved 
! with an explicit scheme. This means shorter timesteps, but with a vertical discretization
! of 30 cm. in the ice timesteps of ~ 1400. seconds will be "tolerated".
! =======================================================================================
      contains
      subroutine icestate_thermf(rt,thermo,restart,dtime,rforce)
      use mod_icestate
      use mod_icestate_diag
      use mod_icestate_tools
      use mod_icestate_fluxes
      use mod_year_info, only : year_info
      use mod_icestate_srfbudget
      use mod_icestate_hpar
      use m_icestate_heatsolve
      use m_icestate_exchange
      use m_icestate_solar
      use m_icestate_prec
      use m_thermf_nersc, only : qsw0
      use mod_common_ice, only : surf_albedo_sum

      IMPLICIT none 



      type(year_info),intent(in) :: rt
      logical,        intent(in) :: thermo
      logical,        intent(in) :: restart
      real*8,         intent(in) :: dtime
      character*5,    intent(in) :: rforce


      ! Local variables
      !****************

      integer :: i,j,l,hk, hh
      real :: rlat, wind,radfl, rainfall, fice_, hice_, qml_ice, 
     &   qtotw1, qsww, rhosw, cpml, qrest, ficetot, omficetot, ustr_i,
     &   fcor, qice_ml,  evap, emnp, flxfac, fX, fL, dmice, dmsnw, 
     &   smlXXX, smlLLL, tmlXXX, tmlLLL, rhoswXXX, rhoswLLL, hmlXXX, 
     &   hmlLLL, cpmlLLL, cpmlXXX, tice_fLLL, tice_fXXX, tice_f, hml, 
     &   sml, tml, qtotlead, mml, surpl_water, tauy, taux, tauyice, 
     &   tauxice, rh, prcp, tair, clouds, fac, cawdir,coalbw_eff, 
     &   tmp_ftot, slp, Xfrac, dsgds, dsgdt,excoeff, tau_trl, 
     &   tau_srl,ssss,ssts,exchng,flg,omflg,sal_corr,hmlEMNPXXX,
     &   smlEMNPXXX,tmlEMNPXXX, enrgy_diff,atmfl,fl_tmp(0:4),vptmpaw,
     &   vptmpw,qtmpaw,qtmpw,qtmpai,vptmpai,rhoair,dhice,tmp_ahice
#if defined(ICEAGE)
     &   ,vice
#endif

      real*8 :: timeref, timeqsw0

      real, dimension(nthick) ::         
     &   qsw,  condt, condb, qatm, qtotice, qbot, qfac, hofusn, hofusni,
     &   inod,qmax, trans, qswi, qsw_surf, qsw_store, qsw_trans, qsup, 
     &   coalbs,qsws, t_old, meltfactor,albi, coalbi,vptmps,qtmps,
     &   atmfl_tmp,snwflg
#if defined (SSNOWD) 
     &    ,hsnwfm,fscov
#endif

      type(t_ice ),dimension(nthick) :: icem, iceXXX, iceLLL, iceddd,  
     &   iceEMNP, ice_init
      type(t_istate_cell) :: icest,icest_tmp,icest_tmp2
      logical      :: rain, snow,errsnw,errflx, lsolar,lexchcoeff
!     &                hmelt_hak, hmelt_solar
      type(t_ice ) :: iceLLL_

      ! Diag !!
      real :: massdiff,saldiff,snwdiff,tmp2,fictot(idm,jdm)
      character(len=2) :: cnn

      logical, parameter :: lsynssr=.false.

      ! Needed for pg - compilers ...
      real, dimension(nthick) :: ones
      logical, dimension(nthick) :: trueones
      ones=1.
      trueones=.true.


!$OMP  PARALLEL DO PRIVATE(i,j) 
!$OMP&SCHEDULE(STATIC,1)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         Iustar  (i,j) = 0.
         Isalflx (i,j) = 0.
         Isalflx2(i,j) = 0.
         Isurflx (i,j) = 0.
         Isswflx (i,j) = 0.
      end do
      end do
!$OMP END PARALLEL DO

      if (.not.thermo) then
         return
      end if 
      errsnw=.false.
      errflx=.false.

      ! Activate Hakkinen & Mellor Horizontal melt - moved to infile
      !hmelt_hak=.true.
      !hmelt_hak=.false.

      ! Activate solar flux -> horizontal melt - moved to infile
      !hmelt_solar=.true.
      !hmelt_solar=.false.

      ! Activate use of exchange coeff in place on "infinite diffusivity"
      lexchcoeff=.true.
      !lexchcoeff=.false.

      if (hmelt_hak.and.mod(rt%ihh,6)==0 .and.rt%iss<dtt .and.mnproc==1)
     &    print *,'Hakkinen horizontal melt. param'

      if (hmelt_hak.and.mod(rt%ihh,6)==0 .and.rt%iss<dtt .and.mnproc==1)
     &    print *,'Solar horizontal melt. param'

      if (.not.lexchcoeff.and.mod(rt%ihh,7)==0 .and.rt%iss<dtt .and.
     &    mnproc==1 ) 
     &    print *,'"Infinite diffusivity" in icestate_thermf'


      if (mnproc==1) 
     &print *,'NB: remember qml_ice change',rt%ihh

      ! Set time from start of year [s]
      !timeref=rt%idd+rt%ihh/24.+rt%iss/3600.

      ! Set salinity and temperature relaxation [s]
      !tau_trl=trelax*86400.
      !tau_srl=srelax*86400.

      ! Set compute-solar-forcing-flag
      lsolar = restart .or. (rt%ihh*3600.+rt%iss<dtt)

      if (lsolar) then
         ! NB - this should be modified to account for average days in
         ! year over a 400-year period (when running with yearflag=0)
         timeqsw0=rt%idd + rt%iss/86400. + rt%ihh/24. + dtt/86400.

Cold     call qsw0(radfl_day,cawdir_day,Iclouds,rlatm,timeqsw0,
Cold &             rt%daysinyear)
         !TODO: clouds used in qsw0 calc from mod_forcing_nersc !
#if defined (DIURNAL)
c --- Coarse reproduction of the diurnal cycle
         call qsw0(radfl_day,cawdir_day,timeqsw0, rt%daysinyear,rt%ihh)
#else
c --- Daily average incoming shortwave radiation    
         call qsw0(radfl_day,cawdir_day,timeqsw0, rt%daysinyear)
#endif
!         These two routines (sun_radfl,qsw0) do the same stuff. Lets stick
!         to qsw0
!
!         if  (mnproc==1) print *,'KAL: Solar flux estimate'
!C$OMP PARALLEL DO PRIVATE(j,l,i,clouds)
!C$OMP&SCHEDULE(STATIC,jblk)
!         do j=1-margin,jj+margin
!         do l=1,isp(j)
!         do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
!
!            ! This should be ok when clouds are climatology!!!
!            clouds  = Iclouds (i,j)
!
!            call sun_radfl(cawdir_day(i,j),radfl_day(i,j),clouds,
!     &         rlatm(i,j),dtime,rforce)
!
!         end do
!         end do
!         end do
!C$OMP END PARALLEL DO
      end if




! --- -------------------------------------------------------
! --- 0 Start loop in space. 
! --- The following is added to variable names :
! ---
! --- XXX   -   Denotes variable calculated for ice-covered
! ---           section of the grid cell.
! ---
! --- LLL   -   Denotes  variable calculations for open-water
! ---           fraction of grid cell.
! --- -------------------------------------------------------

!$OMP  PARALLEL DO PRIVATE(i,j,l,hk, hh, slp, rlat,
!$OMP& wind,radfl, fice_, hice_, qml_ice, qtotw1, qswi, qsww,
!$OMP& qsw, Xfrac, rhosw, cpml, qrest, ficetot, omficetot, 
!$OMP& qice_ml,  evap, emnp, flxfac, fX, fL, dmice, dmsnw, 
!$OMP& smlXXX, smlLLL, tmlXXX, tmlLLL, rhoswXXX, rhoswLLL, 
!$OMP& hmlXXX, hmlLLL, cpmlLLL, cpmlXXX,tice_fLLL,tice_fXXX, 
!$OMP& tice_f, hml, sml, tml, qtotlead, mml, surpl_water, tauy, 
!$OMP& taux, tauyice, tauxice, rh, prcp, tair, clouds, fac, 
!$OMP& cawdir,coalbw_eff, tmp_ftot, condt, condb, qatm, 
!$OMP& qtotice, qbot, qfac, hofusn, hofusni, inod,qmax, trans, 
!$OMP& qsw_surf, qsw_store, qsw_trans, qsup, rain, snow, icem, 
!$OMP& iceXXX, iceLLL, iceddd, iceLLL_, iceEMNP, icest, dsgds,
!$OMP& dsgdt,t_old,albi,coalbi,meltfactor, coalbs, qsws, 
!$OMP& excoeff,ssss,ssts,exchng,ustr_i,fcor, rainfall,ice_init,  
!$OMP& massdiff,saldiff,snwdiff,icest_tmp2,icest_tmp,flg,omflg,
!$OMP& sal_corr, smlEMNPXXX,tmlEMNPXXX,hmlEMNPXXX,enrgy_diff,
!$OMP& atmfl,fl_tmp,vptmps,vptmpai, vptmpaw,vptmpw,qtmpai,qtmpw,
!$OMP& qtmps,qtmpaw,snwflg,tmp2,rhoair,dhice,tmp_ahice) 
!$OMP&SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do l=1,isp(j)
      do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))

! --- =============================================================
! --- Initialize program variables 
! --- =============================================================

      icest    = icestate(i,j)
      icem     = icest%ice
      ice_init = icem
      iceXXX   = icem
      iceLLL   = icem
      iceEMNP  = icem
#if defined(SSNOWD)
      iceddd   = clear_ice(iceddd,nlay,cvsnw)
#else
      iceddd   = clear_ice(iceddd,nlay)
#endif      

      coalbs = 1. - icem%albs

      ! Derived and tmp. variables
      ficetot   = sum(icem%fice)
      omficetot = 1. - ficetot

      ! Seawater properties
      hml       = icest%hml
      sml       = icest%sml
      tml       = icest%tml
      tice_f    = freeze_temp(sml)
      rhosw     = 1000. + sig0(tml-t0deg,sml)
      cpml      = cpsw * rhosw * hml
      !print *,'rhosw',i,j,rhosw,hml,tml,sml

      ! Seawater properties under leads
      smlLLL    = sml     
      tmlLLL    = tml
      hmlLLL    = hml     
      tice_fLLL = tice_f 
      rhoswLLL  = rhosw  
      cpmlLLL   = cpml   

      ! Seawater properties under ice
      smlXXX    = sml
      tmlXXX    = tml
      hmlXXX    = hml
      tice_fXXX = tice_f
      rhoswXXX  = rhosw
      cpmlXXX   = cpml

      ! Latitude  & Coriolis
      rlat   = rlatm(i,j)
      fcor   = 1.45e-4 * sin(rlat)

      ! Atmospheric forcing
      wind    = Iwndspd (i,j)
      rh      = Irelhum (i,j)
      tair    = Iairtmp (i,j)
      prcp    = Iprecip (i,j)
      taux    = Itaux   (i,j) 
      tauy    = Itauy   (i,j)
      clouds  = Iclouds (i,j)
      slp     = Islp    (i,j) * 100 ! mbar -> Pa

 
      ! Changed, now passed directly from EVP to ICESTATE
      tauxice = Itauxice(i,j)
      tauyice = Itauyice(i,j)

      ! Heats of fusion
      hofusn = hofusn0
      hofusni = 1. / hofusn

#if defined (ICEAGE)
      !Ice age increases for each category
      call increment_age(icem,rt)
      vice=sum(icem%fice*icem%hice)
      !ave_age refers to an volume averaged ice age. 
      if(vice>epsil1) then
         ave_age(i,j) = sum(icem%age*icem%fice*icem%hice)/vice
      else
         ave_age(i,j) = 0.
      end if
#endif      

! --- =============================================================
! --- Calculate albedos etc
! --- =============================================================

      ! Compute 24 hrs mean solar irradiance and average solar direction.
      !call sun_radfl(cawdir,radfl,clouds,timeref,rlat)
      ! Solar forcing, estimated above
      if (lsynssr) then
         ! Net downward shortwave flux  - already weighted by clouds
         radfl   = Iswflx   (i,j) 
      else
         radfl   = radfl_day(i,j)  !Also weighted by clouds
      end if

      cawdir  = cawdir_day(i,j)

      ! Effective co-albedo over water with contr from diffuse and direct light
      coalbw_eff =  coalbw * clouds + cawdir * (1.-clouds)

      ! Ice albedo. Albedo for thin ice from Maykut & Perovich (1987). Ice 
      ! albedo changed in case of surface melt.
      meltfactor  = step(icem%tsrf,tice_m)
      albi        = min(albi_max,.08+.44*icem%hice**.28)
      albi        = meltfactor*albi_mlt + (1.-meltfactor)*albi
      coalbi      = 1. - albi

      ! net short wave downward irradiance [w / m^2] (i=ice, s=snow, w=water)
      qswi     = radfl * coalbi
      qsws     = radfl * coalbs
      qsww     = radfl * coalbw_eff

#if defined (SSNOWD)
      ! snow cover fraction depends on snow distribution
      fscov = snow_frac(icem)
      qsw    = (1.-fscov) * qswi + fscov * qsws
      ! Compute grid cell average albedo
      ! Open water is included
       surf_albedo_sum (i,j) = (1.-ficetot)*(1.-coalbw_eff)+            
     &                 sum(fscov*icem%fice*icem%albs)+                 
     &                  sum((1-fscov)*icem%fice*albi)
      ! Compute snow cover fraction over sea ice
       if(ficetot>epsil1) then
         snow_cov_frac(i,j) = sum(fscov*icem%fice)/(ficetot)
      else
         snow_cov_frac(i,j) = 0.
      end if

#else       
      ! Short wave downward irradiance [W/m^2] on surface
      snwflg = step(icem%hsnw,snwlim_trans) 
      qsw    = (1.-snwflg) * qswi + snwflg * qsws
      surf_albedo_sum (i,j) = (1.-ficetot)*(1.-coalbw_eff)+            
     &                 sum(snwflg*icem%fice*icem%albs)+                 
     &                  sum((1-snwflg)*icem%fice*albi)
      
#endif



! --- =============================================================
! --- Handle incoming shortwave radiation  for ice
! --- =============================================================

      ! Shortwave radiation going through ice etc.
      CALL calc_trans(icem,inod,qmax,trans)
      qsw_surf   = (1. - inod) * qsw               ! Heat retained at surface
      qsw_store  = inod * (1. - trans) * qsw       ! Heat stored in ice
      qsw_trans  = inod * trans * qsw              ! Heat transmitted to ml
      icem%qstore= icem%qstore + dtt * qsw_store   ! Update brine reservoir

      ! Modified heat of fusion due to brine pockets. Sets in when 
      ! heat storage is full.
      where (icem%qstore>qmax-epsil1.and.icem%hice>epsil1
     &       .and.icem%fice>epsil1)
         qsw_surf = qsw_surf + (icem%qstore-qmax)/dtt
         hofusn = hofusn0 - icem%qstore/icem%hice
         hofusni = 1. / hofusn
      endwhere

      ! Compute the solar flux that crossed ice  to warm up the mixed layer  ...
      qice_ml    = sum(icem%fice * qsw_trans) / (ficetot + 1e-6)

      ! ...heat the ml under ice accordingly...
      !print *,i,j,tmlXXX,qice_ml,dtt,cpmlXXX,cpsw,hml
      tmlXXX     = tmlXXX + qice_ml  * dtt / cpmlXXX
      
      ! ...and compute mixed layer->ice heat flux, deduce new  ML temperature.
      ! Heat transfer coeff from Holland and Jenkins (1999). 
      ustr_i  = sqrt( sqrt(tauxice*tauxice + tauyice*tauyice)/rhoswXXX)
      ustr_i  = max(windi,ustr_i)  ! Minimum friction velocity is windi
      excoeff = heatexch_coeff(ustr_i,fcor) 
      if (lexchcoeff) then
!        qml_ice = max(0.,min(cpmlXXX/dtt,rhoswXXX*cpsw*excoeff)*
!    &            (tmlXXX-tice_fXXX) )

         ! Note that the heat flux can only act one way - melt..
         ! Numerical effects in the ocean model can lead to a net melt
         ! since the actual temperature can wob slightly around the melting
         ! point. The following correction does not melt ice unless the
         ! temperature is 0.01 degrees above the freezing point
         qml_ice = max(0.,min(cpmlXXX/dtt,rhoswXXX*cpsw*excoeff)*
     &            (tmlXXX-tice_fXXX -.05) )
      else
         qml_ice = cpmlXXX*(tmlXXX-tice_fXXX)/dtt
      end if

      qml_ice = qml_ice*step(ficetot,1e-10)
      

      ! Update ice ml properties
      mml     = hmlXXX * rhoswXXX
      tmlXXX  = tmlXXX - qml_ice*dtt /cpmlXXX
      rhoswXXX= 1000. + sig0(tmlXXX-t0deg,smlXXX)
      hmlXXX  = mml / rhoswXXX
      cpmlXXX = cpsw * rhoswXXX * hmlXXX



#if defined (ENERGY_BALANCE)
      icest_tmp%ice=icem
      icest_tmp%hml=hmlXXX
      icest_tmp%sml=smlXXX
      icest_tmp%tml=tmlXXX
      enrgy_diff=energy_difference(icest,icest_tmp,ficetot)/dtt
     &          -qml_ice*ficetot
      if (enrgy_diff>1e-9.and.ficetot>1e-8) then
         print *,'Qsolar energy_diff: ',enrgy_diff
         print *,'Qsolar Q tranmsmit  ',sum(qsw_trans*icem%fice)
         print *,'Qsolar Q qsw_surf   ',sum(qsw_surf *icem%fice)
         print *,'Qsolar Q qsw_store  ',sum(qsw_store*icem%fice)
         print *,'Qsolar Q qsw -> ml  ',qice_ml*ficetot
         print *,'Qsolar Q qml -> ice ',qml_ice*ficetot
      end if
#endif


! --- =============================================================
! --- Calculate Precipitation -- Update ML properties
! --- =============================================================

      ! Calculate precipitation [m / (m^2)]
      rainfall   = prcp *  dtt ! In New hycom, prec is already in m/s

      ! Rainfall over snow/ice
      call icestate_prec(icem,rainfall,tair,rain,snow,rhoswXXX,
     &                   surpl_water,sal_corr)
      if (rain.and.snow) then
         print *,'Na-na-nah ...'
         errsnw=.true.
      end if

      ! used in mass/salinity balance calcs
      snwdiff = sum(icem%fice*icem%hsnw*icem%rhosnw) - 
     &          sum(icest%ice%fice*icest%ice%hsnw*icest%ice%rhosnw) 
     &          + rainfall*1000.*ficetot-surpl_water

      ! Evaporation over leads [m /(s m^2)]
      evap      =  evaporation(tair,tml,wind,slp,rh,rhoswLLL)

      ! Evaporation minus precipitation over leads [ kg / (m^2 * m^2)]
      emnp      =  ( evap * dtt - rainfall ) * 1000.

      iceEMNP=icem ! Keep for later

      ! Modification of mixed layer under lead due to evaporation and precipitation. 
      mml       = hmlLLL * rhoswLLL
      smlLLL    = smlLLL * mml /(mml - emnp )
      rhoswLLL  = 1000. + sig0(tmlLLL-t0deg,smlLLL)
      hmlLLL    = (mml - emnp ) / rhoswLLL
      tice_fLLL = freeze_temp(smlLLL)
      cpmlLLL   = cpsw * rhoswLLL * hmlLLL


#if defined (MASS_SALINITY_BALANCE)
      if (omficetot>1e-6) then
        icest_tmp%ice=iceLLL
        icest_tmp%hml=hmlLLL
        icest_tmp%sml=smlLLL
        icest_tmp%tml=tmlLLL
        saldiff = salinity_difference(icest,icest_tmp,omficetot)
        massdiff = mass_difference(icest,icest_tmp,omficetot)
     &           +emnp*omficetot
        if (abs(saldiff)>1e-9  .or.abs(massdiff)>1e-10) then 
           print *,'EVAP1:omfice Salinity difference ', saldiff
           print *,'EVAP1:omfice MASS     difference ', massdiff
           print *,'EVAP1:omfice emnp                ', emnp*omficetot
        end if
      end if
#endif

#if defined (ENERGY_BALANCE)
      icest_tmp%ice=iceLLL
      icest_tmp%hml=hmlLLL
      icest_tmp%sml=smlLLL
      icest_tmp%tml=tmlLLL
      enrgy_diff=energy_difference(icest,icest_tmp,omficetot)/dtt 
     &           + emnp*omficetot*tmlLLL*cpsw/dtt
      if (enrgy_diff>1e-9.and.omficetot>1e-8) then
         print *,'QEVAP (ow) energy_diff: ',enrgy_diff
         print *,'QEVAP (ow) new ms enrg: ',
     &           emnp*omficetot*tmlLLL*cpsw/dtt
      end if
#endif


      ! Surplus water after precipitation on ice sent under ice
      mml       = hmlXXX * rhoswXXX
      smlXXX    = (smlXXX * mml - sal_corr) /(mml + surpl_water )
      rhoswXXX  = 1000. + sig0(tmlXXX-t0deg,smlXXX)
      hmlXXX    = ( mml + surpl_water ) / rhoswXXX
      tice_fXXX = freeze_temp(smlXXX)
      cpmlXXX   = cpsw * rhoswXXX * hmlXXX




#if defined (MASS_SALINITY_BALANCE)
      icest_tmp2=icest
      icest_tmp%ice=iceEMNP 
      icest_tmp%hml=hmlXXX
      icest_tmp%sml=smlXXX
      icest_tmp%tml=tmlXXX
      saldiff = salinity_difference(icest_tmp2,icest_tmp,ficetot)
      massdiff = mass_difference(icest_tmp2,icest_tmp,ficetot) 
     &         - rainfall*1000. * ficetot
      if (abs(saldiff)>1e-9 .or.abs(massdiff)>1e-10 .and. ficetot>1e-6) then 
        print *,'EVAP :fice Salinity difference ', saldiff
        print *,'EVAP :fice MASS     difference ', massdiff
        print *,'EVAP :fice rainfall            ',rainfall*1000.*ficetot
        print *,'EVAP :fice icediff snwdiff     ',
     &  sum(icem%hice-icest%ice%hice)*rhoice,     
     &  sum(icem%hsnw*icem%rhosnw-icest%ice%hsnw*icest%ice%rhosnw)
        print *,'EVAP :fice mldiff              ',mml-rhosw*hml
        print *,'EVAP : surplus_water           ',surpl_water
        print *,'EVAP : ice th diff             ', 
     &          sum(icest_tmp%ice%hice-icest_tmp2%ice%hice)
        print *,'EVAP :fice ficetot             ',ficetot
        print *,'EVAP :rhosnw                   ',icem%rhosnw
        print *,'EVAP :fice                     ',icem%fice
        print *,'EVAP :hsnw                     ',icem%hsnw
        print *,'EVAP :tair                     ',tair
        print *, '------------------------------'
      end if
#endif


      iceEMNP=icem ! Keep for later
      hmlEMNPXXX=hmlXXX
      smlEMNPXXX=smlXXX
      tmlEMNPXXX=tmlXXX

#if defined (ENERGY_BALANCE)
      icest_tmp%ice=icem
      icest_tmp%hml=hmlXXX
      icest_tmp%sml=smlXXX
      icest_tmp%tml=tmlXXX
      enrgy_diff=energy_difference(icest,icest_tmp,ficetot)/dtt       
     &         - qml_ice*ficetot                                      
     &         - surpl_water*tmlXXX*cpsw*ficetot/dtt                 
     &         + hofusn0                                             
     &          *sum(icest%ice%fice*icest%ice%hsnw*icest%ice%rhosnw -
     &               iceEMNP%fice*iceEMNP%hsnw*iceEMNP%rhosnw ) /   
     &               (rhosw*dtt)
      if (enrgy_diff>1e-9.and.ficetot>1e-8) then
        print *,'QEVAP (fice) energy_diff: ',enrgy_diff
        print *,'QEVAP (fice) energy_diff: ',
     &     energy_difference(icest,icest_tmp,ficetot)/dtt
        print *,'QEVAP (fice) new ms enrg (snow): ',          
     &       hofusn0                                             
     &      *sum(icest%ice%fice*icest%ice%hsnw*icest%ice%rhosnw -
     &           iceEMNP%fice*iceEMNP%hsnw*iceEMNP%rhosnw ) /    
     &           (rhosw*dtt*ficetot)
        print *,'QEVAP (fice) new ms enrg (surpl): ',            
     &         - surpl_water*tmlXXX*cpsw*ficetot/dtt             
        print *,'QEVAP (fice) new ms enrg (qsw  ): ',           
     &         - qml_ice*ficetot                               
        print *,'QEVAP (fice) new ms enrg (ice ): ',          
     &      sum(hofusn*(icest%ice%fice*icest%ice%hice
     &          -iceEMNP%fice*iceEMNP%hice))/dtt
        print *
      end if
#endif

! === ===============================================================================
! --- Calculate surface energy budget and vertical heat conduction
! === ===============================================================================



      condt   = 0.
      condb   = 0.
      qbot    = 0.
      qatm    = 0.
      qtotice = 0.
      qfac    = 0.
      qsup    = 0.
      t_old   = 0.

      !print *,'prior vcond',i,j

      where (icem%fice>epsil1) 
         
         ! Calculate new surface temperature based on previous timestep
         t_old = icem%tsrf
C        icem%tsrf = srftemp_ice(icem,t_old,rh,tair,wind,qsw_surf,slp,
C    &                           clouds,rlat,tml)
         icem%tsrf = srftemp_ice(icem,t_old,rh*ones,tair*ones,wind*ones,
     &                           qsw_surf,slp*ones,clouds*ones,
     &                           rlat*ones,tml*ones)

         ! Correct surface temp and adjust heat reservoir if necessary.
         where (icem%tsrf>tice_m)
            icem%tsrf=tice_m
         elsewhere ( icem%qstore>epsil1 .and. icem%tsrf < tice_m) 

            ! Let us use the heat store to raise surface temperature
            qsup      = icem%qstore/dtt
C           icem%tsrf = srftemp_ice(icem,t_old,rh,tair,
C    &                              wind,qsw_surf+qsup,
C    &                              slp,clouds,rlat,
C    &                              tml)
            icem%tsrf = srftemp_ice(icem,t_old,rh*ones,tair*ones,
     &                              wind*ones,qsw_surf+qsup,
     &                              slp*ones,clouds*ones,rlat*ones,
     &                              tml*ones)
            where (icem%tsrf>tice_m-epsil1)

               ! Whoops! Too much ... Keep at melting temperature, and 
               ! use only part of qstore
               icem%tsrf = tice_m
C              qsup      = -atm_heatflux(icem%tsrf,rh,tair,wind,
C    &                             qsw_surf,slp,clouds,rlat,ice=.true.)
               qsup      = -atm_heatflux(icem%tsrf,rh*ones,tair*ones,
     &                                   wind*ones,qsw_surf,slp*ones,
     &                                   clouds*ones,rlat*ones,
     &                                   ice=trueones)
C              qsup      = qsup - cond_heatflux(icem,tmlXXX)
C              qsup      = qsup - inertia_heatflux(icem,t_old)
               qsup      = qsup - cond_heatflux(icem,tmlXXX*ones)
               qsup      = qsup - inertia_heatflux(icem,t_old*ones)
            end where
         end where
         
         ! Qstore is modified by flux qsup
         icem%qstore  = icem%qstore - qsup*dtt

         ! Atmospheric fluxes
C        qatm    = atm_heatflux(icem%tsrf,rh,tair,wind,
C                               qsw_surf,slp,clouds,
C                               rlat,ice=.true.)
         qatm    = atm_heatflux(icem%tsrf,rh*ones,tair*ones,wind*ones,
     &                          qsw_surf,slp*ones,clouds*ones,
     &                          rlat*ones,ice=trueones)

         ! Surface and bottom boundary conditions
         ! NB ! cond_heatflux is for surfa
         where (icem%fice>epsil1)
            where (thin)
C              condb = -cond_heatflux(icem,tmlXXX)
C              condt = -condb  + inertia_heatflux(icem,t_old)
               condb = -cond_heatflux(icem,tmlXXX*ones)
               condt = -condb  + inertia_heatflux(icem,t_old*ones)
            elsewhere 
C              condt = inertia_heatflux(icem,t_old) 
C    &                +cond_heatflux(icem,tmlXXX)
               condt = inertia_heatflux(icem,t_old*ones) 
     &                +cond_heatflux(icem,tmlXXX*ones)
               condb = rkice*(icem%vtp(1)-tice_fXXX)*2*
     &                 icem%nlay/icem%hice
            endwhere
         endwhere
         
         ! Total heat flux towards top of ice
         qtotice = condt + qsup + qatm 

      end where ! Surface heat budget

            
      ! Set up linear temperature profile in the thin ice.
      forall (hk=1:hklim , icem(hk)%fice>epsil1)
         icem(hk) = vtplin(icem(hk),icem(hk)%tsrf,tice_fXXX)
      end forall 



      ! Set up temperature profile in thick ice
      do hk=hklim+1,nthick
      if (icem(hk)%fice>epsil1) then

         condt(hk)      = -condt(hk) ! heatsolve uses different positive directions
         icem(hk)%vtp(:) = 
     &      icestate_heatsolve(icem(hk),dtt,condt(hk), condb(hk))
         condt(hk)       = -condt(hk)  ! heatsolve uses different positive directions
      endif
      enddo


      ! Primarily used in daily average calculations
      icestate_swfl(i,j,:)=qsw_surf 
      icestate_swtr(i,j,:)=qsw_trans
      icestate_trb (i,j,:)=
     &   turb_flux(tair*ones,wind*ones,clouds*ones,icem%tsrf,
     &             slp*ones,rh*ones,ice=trueones)
      icestate_lw  (i,j,:)=
     &   longwave_flux(tair*ones,clouds*ones,rlat*ones,
     &                 icem%tsrf,ice=trueones)
      icestate_brfl(i,j,:)=-qsw_store+qsup
      icestate_mlfl(i,j,:)=qml_ice 
      icestate_ctop(i,j,:)=condt
      icestate_cbot(i,j,:)=condb
      icestate_ntop(i,j,:)=qtotice
      icestate_nbot(i,j,:)=condb + qml_ice
      

! === ===============================================================================
! --- Calculate New ice/snow thickness based on the heat fluxes
! === ===============================================================================


      ! Update snow albedo
      call snow_albedo_update(rain,snow,rainfall,qtotice,icem%albs)


      ! First we use  qtotice on snow surface to melt it
      where (qtotice > epsil1) 
#if defined (SSNOWD)
         prov=qtotice
         ! Snow is melting
         hsnwfm     = dtt * qtotice * hofusni0 *  rhoice / icem%rhosnw

         !Snow on the ice
           icem%hmelt = icem%hmelt+hsnwfm

           ! Compute the average snow depth in the grid cell
           icem%hsnw  = average_depth(icem)
           ! When snow cover becomes too thin, reset melt and acccumulation variable to 0
           where( icem%hsnw < 0.005)
              !Removed snow is added to the ice (energy conservation)
              icem%hice=icem%hice+icem%hsnw*icem%rhosnw/rhoice
              icem%hsnw  = 0.
              icem%hprcp = 0.
              icem%hmelt = 0.
           elsewhere
           ! qtotice has been used to melt snow
            qtotice = 0.
           end where
#else

         ! qtotice > 0 means that there is surface melting -- apply it ...
         icem%hsnw    = icem%hsnw - dtt * qtotice * hofusni0 *  rhoice 
     &                              / icem%rhosnw

         ! Negative ice thickness means only part of qtotice was used for snow melt
         qtotice = max(0.,- hofusn0 * icem%rhosnw * icem%hsnw 
     &             /(rhoice * dtt))
         icem%hsnw    = max(0.,icem%hsnw)
#endif         
         where (icem%hsnw<epsil1) icem%rhosnw=rhosnwmin
      end where



      ! Then we start melting ice ... 
      where (icem%fice>epsil1)

         ! Bottom flux
         qbot  = condb

         ! Global available heat flux
         qfac  = qtotice + qbot + qml_ice

         ! Apply melting or freezing (qfac)
         icem%hice  = icem%hice - dtt * qfac * hofusni

         ! We melted too much ice. Try melting all the snow first
         where (icem%hice < epsil1)

            ! Ice melting is cancelled :
            icem%hice = icem%hice + dtt * qfac * hofusni

            ! Snow melting instead :
#if defined (SSNOWD)
            hsnwfm     = dtt * qfac * hofusni0 *  rhoice / icem%rhosnw
            where(icem%hprcp>epsil1)
             ! Snow is melting
             icem%hmelt = icem%hmelt+hsnwfm
             ! Compute the average snow depth in the grid cell
             icem%hsnw  = average_depth(icem)
             ! When snow cover becomes too thin, reset melt and acccumulation variable to 0
             where( icem%hsnw < 0.005)
                ! Removed snow is added to the ice (energy conservation)
                icem%hice=icem%hice+icem%hsnw*icem%rhosnw/rhoice
                icem%hsnw  = 0.
                icem%hprcp = 0.
                icem%hmelt = 0.
             end where
             qfac   = max(0.,- hofusn0 * icem%rhosnw * icem%hsnw 
     &               / (rhoice * dtt) )
          elsewhere
              qfac   = hofusn0 * icem%rhosnw * hsnwfm
     &               / (rhoice * dtt)
          endwhere
#else             
            icem%hsnw  = icem%hsnw - dtt * qfac * hofusni0 * rhoice 
     &                               / icem%rhosnw   
            qfac   = max(0.,- hofusn0 * icem%rhosnw * icem%hsnw 
     &               / (rhoice * dtt) )
            icem%hsnw   = max(0.,icem%hsnw)
#endif            
            where (icem%hsnw<epsil1) icem%rhosnw = rhosnwmin 

            ! Proceed melting ice if any energy is left
            icem%hice = icem%hice - dtt * qfac * hofusni

         endwhere
      endwhere

      ! The ice can now be negative, use remaining energy towards the ML
      qtotw1 = - sum( min(0.,icem%hice) * hofusn * icem%fice ) / dtt
      icem%hice   = max(0.,icem%hice)
      tmlXXX  = tmlXXX + qtotw1*dtt /cpmlXXX


      ! Keep ice state properties after vertical freeze/melt
      iceXXX = icem

      snwdiff = sum(iceEMNP%fice*iceEMNP%hsnw*iceEMNP%rhosnw) - 
     &          sum(icem%fice*icem%hsnw*icem%rhosnw)

! === ===============================================================================
! --- Calculate New ML properties  --  Combine ice classes if necessary
! === ===============================================================================

      
      flg = step(ficetot,epsil1)
      omflg = step(omficetot,epsil1)

      ! Update the ml quantities under ice. Ice changes relative to last operation 
      ! (which was precipitation, iceEMNP) 
      dmice     = rhoice * ( sum(iceXXX%hice * iceXXX%fice) - 
     &                       sum(iceEMNP%hice   * iceEMNP%fice) )
      dmsnw     = min(0.,
     &             sum(iceXXX%fice  * iceXXX%hsnw  * iceXXX%rhosnw )
     &            -sum(iceEMNP%fice * iceEMNP%hsnw * iceEMNP%rhosnw))
!Nasty ... ficetot should be accounted for -- was not done before ...
      dmice     = flg*dmice/(ficetot+1.-flg)  ! For ice unit area
      dmsnw     = flg*dmsnw/(ficetot+1.-flg)  ! For ice unit area
      mml       = hmlXXX*rhoswXXX
      smlXXX    = (mml*smlXXX - dmice*sice) /  (mml - dmice - dmsnw)
      rhoswXXX  = 1000. + sig0(tmlXXX-t0deg,smlXXX)
      hmlXXX    = (mml - dmice - dmsnw) / rhoswXXX
      tice_fXXX = freeze_temp(smlXXX)
      cpmlXXX   = cpsw * rhoswXXX * hmlXXX


      ! Empty and join classes if necessary
#if defined(SSNOWD)
      iceddd   = clear_ice(iceddd,nlay,cvsnw)
#else
      iceddd   = clear_ice(iceddd,nlay)
#endif      
      call combine_ice(iceddd,iceXXX,tmlXXX)
      iceXXX=iceddd

      ! Finally adjust qstore -- some of it could be used when 
      ! we melted ice -- through hofusn
      CALL calc_trans(iceXXX,inod,qmax,trans)
      where (iceXXX%qstore>qmax-epsil1)
         iceXXX%qstore=qmax
      end where


#if defined (MASS_SALINITY_BALANCE)
      if (ficetot>1e-6) then
        icest_tmp2=icest
        icest_tmp%ice=iceXXX
        icest_tmp%hml=hmlXXX
        icest_tmp%sml=smlXXX
        icest_tmp%tml=tmlXXX
        saldiff = salinity_difference(icest_tmp2,icest_tmp,ficetot)
        massdiff = mass_difference(icest_tmp2,icest_tmp,ficetot)
     &     -rainfall*1000.*ficetot
        if (abs(saldiff)>1e-8  .or.abs(massdiff)>1e-10) then 
           print *,'VICE :Salinity difference ', saldiff
           print *,'VICE :MASS     difference ', massdiff
           print *,'VICE :rainfall            ',rainfall*1000.*ficetot
           print *,'VICE :fice icediff snwdiff     ',
     &       sum(iceXXX%hice*iceXXX%fice-icest%ice%hice*icest%ice%fice)
     &       *rhoice,    
     &       sum(iceXXX%hsnw*iceXXX%rhosnw-
     &           icest%ice%hsnw*icest%ice%rhosnw)
           print *,'VICE :fice mldiff ',
     &             ficetot*(rhoswXXX*hmlXXX-hml*rhosw)
        end if
      end if
#endif

! === ===============================================================================
! --- Start calculating parametrization of "horizontal" heat flux                          
! === ===============================================================================

      ! iceLLL = icesheet after vertical melt of ice (iceXXX)
      icem   = iceXXX

      ! Surface energy balance in lead. 
      qtotlead = atm_heatflux(tml,rh,tair,wind,qsww,slp,clouds,
     &   rlat,.false.)


      ! used in average calculations
      icestate_lead_tot(i,j)= qtotlead  
      icestate_lead_sw (i,j)= qsww
      icestate_lead_lw (i,j)= 
     &   longwave_flux(tair,clouds,rlat,tmlLLL,.false.)
      icestate_lead_trb(i,j)= 
     &   turb_flux(tair,wind,clouds,tmlLLL,slp,rh,.false.)
      

      ! New lead ML temperature
      tmlLLL = tmlLLL + qtotlead * dtt / cpmlLLL

      ! (Possible) new ice is initialized
      hice_=0.
      fice_=0.

      ! If new temperature below freezing, adjust and freeze some ice
      if (tmlLLL < tice_fLLL) then
         
         ! Energy available for freezing
         !qrest = (tmlLLL - tice_fLLL) * cpml / dtt
         qrest = (tmlLLL - tice_fLLL) * cpmlLLL / dtt

         ! Reset ML temperature
         tmlLLL = tice_fLLL
        
         ! Call the "horizontal" freezing routine
         call icestate_hfreeze(icem,qrest,omficetot,fice_,hice_)



      ! Case of lateral melting of the ice and snow
      else if ( hmelt_hak .or. hmelt_solar) then 

         if (hmelt_hak) then

            ! Mass reduction from vertical melt
            dhice = 
     &         sum(ice_init%fice*ice_init%hice-iceXXX%fice*iceXXX%hice)
            tmp_ftot = sum(icem%fice)
            tmp_ahice= sum(icem%fice*icem%hice)

            if (tmp_ftot>1e-2.and.tmp_ahice>1e-2) then 
               ! Call the "horizontal" melting routine -- Hakkinen & Mellor ! param (1990)
               call icestate_hmelt_hak(icem,smlLLL,tmlLLL,cpmlLLL,
     &                                 ficetot,hofusn,dhice)
            end if
         end if

         if (hmelt_solar) then

            if (qtotlead>epsil1) then 
               ! Call the "horizontal" melting routine
               call icestate_hmelt(icem,tmlLLL,cpmlLLL,ficetot,
     &                             qtotlead,hofusn)
            end if

            end if

            ! Finally, if melting calc. resulted in ice conc below fice_min, adjust fice
            tmp_ftot = sum(icem%fice)
            if (tmp_ftot>epsil1 .and. tmp_ftot<fice_min .and. 
     &          sum(icem%fice*icem%hice) > epsil1 ) then
               icem%hice = icem%hice * icem%fice * tmp_ftot / fice_min
               icem%hsnw = icem%hsnw * icem%fice * tmp_ftot / fice_min
               icem%fice =             icem%fice * fice_min / tmp_ftot
            end if

         endif 

! === ===============================================================================
! --- Combine ice classes that have fallen into the same category
! === ===============================================================================

      ! Clear iceclasses if necessary
      do hk = 1,nthick
         if (icem(hk)%fice < epsil1)
#if defined(SSNOWD)
     &       icem(hk) = clear_ice(icem(hk),nlay(hk),cvsnw(hk))
#else
     &       icem(hk) = clear_ice(icem(hk),nlay(hk))
#endif          
      enddo

      ! Merge new ice with an existing iceclass
      hh = nthick
      DO while (hh>=1)

        IF (hice_>thickl(hh)) then
           !print *,icem(hh)%fice*icem(hh)%hsnw*icem(hh)%rhosnw
#if defined(SSNOWD)
           iceLLL_      = clear_ice(iceLLL_,nlay(hh),cvsnw(hh))
#else
           iceLLL_      = clear_ice(iceLLL_,nlay(hh))
#endif           
           iceLLL_%fice = fice_
           iceLLL_%hice = hice_
           iceLLL_%tsrf = tair
#if defined(ICEAGE)
           ! Age of new ice
           iceLLL_%age  = dtt/(86400.*real(rt%daysinyear))
#endif            

           ! Let vtp of new ice be linear.
           iceLLL_ = VtpLin(iceLLL_,tair,tice_fLLL)

           ! Add iceclasses together, (+) is an operator defined for type t_ice
           icem(hh)  = icem(hh)+iceLLL_
           hh = 0 ! exit inner loop

        end if

        hh=hh-1
      END DO

      ! The ice frozen in leads can change the distribution of ice in a gridcell,
      ! the iceclasses are updated here.
#if defined(SSNOWD)
      iceddd   = clear_ice(iceddd,nlay,cvsnw)
#else
      iceddd   = clear_ice(iceddd,nlay)
#endif       
      call combine_ice(iceddd,icem,tmlLLL)
      icem   = iceddd
      iceLLL = iceddd

      ! Finally adjust qstore -- some of it could be used when 
      ! we melted ice -- through hofusn
      CALL calc_trans(iceLLL,inod,qmax,trans)
      where (iceLLL%qstore>qmax-epsil1)
         iceLLL%qstore=qmax
      end where


! === ===============================================================================
! --- Update the ML properties                                        
! === ===============================================================================

      flg = step(ficetot,epsil1)
      omflg = step(omficetot,epsil1)

      ! Update the ml quantities under the lead if lat freezing occurs
      mml       = hmlLLL*rhoswLLL
      dmice     = rhoice *  max(0., sum( iceLLL%fice * iceLLL%hice ) -  
     &                              sum( iceXXX%fice * iceXXX%hice ) )
      dmice     = omflg*dmice/(omficetot+1.-omflg)  ! For ice unit area
      !dmice     = omflg*rhoice*fice_*hice_/(1.-omflg+omficetot) !For OW unit area
      smlLLL    = ( mml*smlLLL - dmice*sice ) / (mml - dmice)
      rhoswLLL  = 1000. + sig0(tmlLLL-t0deg,smlLLL)
      hmlLLL    = (mml - dmice ) / rhoswLLL
      tice_fLLL = freeze_temp(smlLLL)
      cpmlLLL   = cpsw * rhoswLLL * hmlLLL

#if defined (MASS_SALINITY_BALANCE)
      if (omficetot>1e-6.and.dmice>epsil1) then
        icest_tmp2=icest
        icest_tmp%ice=iceLLL
        icest_tmp2%ice%fice=0.
        icest_tmp%ice%fice=0.
        icest_tmp%hml=hmlLLL
        icest_tmp%sml=smlLLL
        icest_tmp%tml=tmlLLL
        saldiff = salinity_difference(icest_tmp2,icest_tmp,omficetot)
     &           +dmice*sice*omficetot
        massdiff= mass_difference(icest_tmp2,icest_tmp,omficetot)
     &           +emnp*omficetot+dmice*omficetot
        if (abs(saldiff)>1e-9  .or.abs(massdiff)>1e-10) then 
           print *,'HICE1:Salinity difference ', saldiff
           print *,'HICE1:MASS     difference ', massdiff
           print *,'HICE1:rainfall            ',emnp*omficetot
           print *,'HICE1:dmice,dmsnw         ',dmice*omficetot
           print *
        end if
      end if
#endif


      ! In addition the lateral melt influences conditions initially under ice
      mml       = hmlXXX*rhoswXXX
      dmice     = rhoice *  min(0., sum( icem%fice * icem%hice ) - 
     &                              sum( iceXXX%fice * iceXXX%hice ))
      dmsnw     = min(0., sum( icem%fice   *icem%hsnw   *icem%rhosnw )
     &                   -sum( iceXXX%fice *iceXXX%hsnw *iceXXX%rhosnw))
      dmice     = flg*dmice/(ficetot+1.-flg)  ! For ice unit area
      dmsnw     = flg*dmsnw/(ficetot+1.-flg)  ! For ice unit area
      smlXXX    = ( mml*smlXXX - dmice*sice) /  ( mml - dmice - dmsnw)
      rhoswXXX  = 1000. + sig0(tmlXXX-t0deg,smlXXX)
      hmlXXX    = ( mml - dmice - dmsnw ) / rhoswXXX
      tice_fXXX = freeze_temp(smlXXX)
      cpmlXXX   = cpsw * rhoswXXX * hmlXXX



#if defined (MASS_SALINITY_BALANCE)
      if (ficetot>1e-6.and.dmice<-epsil1) then
        icest_tmp2=icest
        icest_tmp%ice=iceLLL
        icest_tmp%hml=hmlXXX
        icest_tmp%sml=smlXXX
        icest_tmp%tml=tmlXXX
        saldiff = salinity_difference(icest_tmp2,icest_tmp,ficetot)
        massdiff = mass_difference(icest_tmp2,icest_tmp,ficetot)
     &            -rainfall*1000.*ficetot
        if (abs(saldiff)>1e-9  .or.abs(massdiff)>1e-10) then 
           print *,'HICE2:Salinity difference ', saldiff
           print *,'HICE2:MASS     difference ', massdiff
           print *,'HICE2:rainfall            ',rainfall*1000.*ficetot
           print *,'HICE2:dmice,dmsnw         ',dmice
        end if
      end if
#endif






! === ===================================================================
! === Step 4:
! === Final computations of grid cell thermodynamic state from ice-covered
! === (XXX) and open water (LLL) components.
! === ===================================================================
! --- 4.2 Final computation of ml properties
! --- ------------------------------------------------------------------------

      ! Area weighted salinities, temperatures and densities
      fX    = ficetot   * rhoswXXX * hmlXXX
      fL    = omficetot * rhoswLLL * hmlLLL
      sml   = (fX * smlXXX + fL * smlLLL) / (fX+fL)
      tml   = (fX * tmlXXX + fL * tmlLLL) / (fX+fL)
      rhosw = 1000. + sig0(tml-t0deg,sml)
      hml   = (fX + fL) / rhosw
      cpml  = cpsw * rhosw * hml



! --- ------------------------------------------------------------------------
! --- 4.3 Derive salinity and buoyancy fluxes used in HYCOM and MICOM routines
! --- 
! --- surflx     [W/(m^2)]
! --- salflx     [gr-salt/(m^2 sec)]
! --- buoyflx    [m^2/sec^3]
! ---
! --- MICOM: Positive flux -> ml is made cooler/fresher
! --- HYCOM: Positive flux -> ml is made warmer/saltier
! ---
! --- Buoyancy flux direction is the same for both models. If buoyfl>0 then the 
! --- coloumn is destabilized. Note that dimensions and direction of the 
! --- fluxes are taken care of in interface routines later.
! --- ------------------------------------------------------------------------

      ! Coeff used in calculations below. 
      flxfac= rhoref*hml/dtt

      ! Sigma differentiated by salinity  and temperature [q/(m^3* (K|ppm)]
      dsgds = dsig0ds(icest%tml-t0deg,icest%sml)  ! dsigds from HYCOM
      dsgdt = dsig0dt(icest%tml-t0deg,icest%sml)  ! dsigdt from HYCOM

      ! Total heat and salinity flux [W/m^2 | gr-salt/(m^2 * sec)]
      Isurflx(i,j)=(tml-icest%tml)*flxfac*cpsw
      Isalflx(i,j)=(sml-icest%sml)*flxfac

      Isalflx2(i,j)=(sml*hml*rhosw 
     &             -icest%sml*icest%hml*
     &              (1000.+sig0(icest%tml-t0deg,icest%sml)))/dtt

      ! The above is salinity LOST to the ice (due to mass loss), so for
      ! salinity to be conserved we must reverse the direction
      Isalflx2(i,j) = - Isalflx2(i,j)



      icestate_grw (i,j,:) = (-icest%ice%hice*icest%ice%fice 
     &                     + iceLLL%hice*iceLLL%fice) *86400/dtt
      icestate_lgrw(i,j,:)  = (-icest%ice%fice + iceLLL%fice)
     &                     *86400/dtt


      ! Flux contribution from shortwave fluxes 
      Isswflx(i,j) = omficetot*qsww

      ! Turbulent Kinetic Energy due to friction from ice and atmosphere [m/s]
      Iustar(i,j) = ustr_i*ficetot 
     &             +sqrt( sqrt(taux*taux+tauy*tauy) /  rhosw )*omficetot

         
! === ===================================================================
! === Step 5:
! === Surface relaxation     -- Handled by hycom thermf
! === ===================================================================
! ---    KAL --- All buoyancy flux AND relaxation calcs performed in HYCOM


! === ===================================================================
! === Step 6:
! === Finally update thermodynamic model values
! === ===================================================================



      ! Final state of ice in gridcell
      icestate(i,j)%ice(:)   = iceLLL(:)

      ! Seawater properties
      icestate(i,j)%hml = hml
      icestate(i,j)%sml = sml
      icestate(i,j)%tml = tml

     
#if defined (MASS_SALINITY_BALANCE)
      ! Final Diagnostic tests for salinity and mass
      saldiff = salinity_difference(icest,icestate(i,j),1.)
      massdiff= mass_difference(icest,icestate(i,j),1.)+
     &          evap*dtt*1000.*omficetot - rainfall*1000.

      if (abs(saldiff)>1e-8 .or. abs(massdiff)>1e-10) then
         print '(a,2i4,a,e12.5)', 'FINAL -- Differences at :', i,j, 
     &         '   salinity  ',saldiff
         print '(a,2i4,a,e12.5)','FINAL -- Differences at :',i,j, 
     &         '   mass-emnp ',massdiff

         print *,'Icediff:',rhoice*(sum(icestate(i,j)%ice%fice*
     &                                  icestate(i,j)%ice%hice) - 
     &                              sum(icest%ice%fice*icest%ice%hice) )
         print *,'snwdiff:',
     &      sum(icestate(i,j)%ice%fice*icestate(i,j)%ice%hsnw*
     &          icestate(i,j)%ice%rhosnw) 
     &     -sum(icest%ice%fice*icest%ice%hsnw*icest%ice%rhosnw)
         print *,'mldiff :',
     &      hml*(1000.+sig0(tml-t0deg,sml))
     &      -icest%hml*(1000.+sig0(icest%tml-t0deg,icest%sml))
         print *,'evap   :',evap*dtt*1000.*omficetot
         print *,'rain   :',rainfall*1000.
         print *,'surpl_w:',surpl_water*ficetot
         print *,'sal_crr:',sal_corr*ficetot
         print *
      end if
#endif



#if defined (ENERGY_BALANCE)
      enrgy_diff=energy_difference(icest,icestate(i,j),1.)/dtt
      atmfl_tmp=atm_heatflux(icest%ice%tsrf,rh,tair,wind,qsw,slp,
     &                       clouds,rlat,.true.)
      atmfl=sum(atmfl_tmp*icest%ice%fice)
      atmfl=atmfl+omficetot*
     & atm_heatflux(icest%tml,rh,tair,wind,qsww,slp,clouds,rlat,.false.)

      do hk=1,nthick
         ! Sensible and latent heat transfer
         ! Vapour pressure and humidity of the air over the surface
         vptmps(hk)= vapp(aice,bice,icest%ice(hk)%tsrf)
      end do
      qtmps = humid(slp,vptmps)
      vptmpai= vapp(aice,bice,tair)
      qtmpai = humid(slp,vptmpai)
      vptmpaw= vapp(awater,bwater,tair)
      qtmpaw=humid(slp,vptmpaw)
      vptmpw= vapp(awater,bwater,tml)
      qtmpw=humid(slp,vptmpw)
      rhoair  = slp / (gasconst * tair)

      atmfl=0.
      do hk=1,nthick
         fl_tmp = turbulent(icest%ice(hk)%tsrf,qtmps(hk),tair,qtmpai,
     &                      rhoair,wind)
         atmfl=atmfl+(fl_tmp(0)+fl_tmp(1)*icest%ice(hk)%tsrf)
     &               *icest%ice(hk)%fice
      end do
      fl_tmp = turbulent(icest%tml,qtmpw,tair,qtmpaw,rhoair,wind)*
     &         omficetot
      atmfl=atmfl+(fl_tmp(0)+fl_tmp(1)*icest%tml)*omficetot


      atmfl=0.
      do hk=1,nthick
         fl_tmp=longwave(tair,clouds,vptmpai,rlat)
         atmfl=atmfl+(fl_tmp(0)+fl_tmp(1)*icest%ice(hk)%tsrf)
     &               *icest%ice(hk)%fice
      end do
      fl_tmp=longwave(tair,clouds,vptmpaw,rlat)
      atmfl=atmfl+(fl_tmp(0)+fl_tmp(1)*icest%tml)*omficetot
#endif


      ! Error check -- qstore
      CALL calc_trans(icestate(i,j)%ice,inod,qmax,trans)
      if (any(icestate(i,j)%ice%qstore>qmax-epsil1.and.
     &     icestate(i,j)%ice%fice>epsil1)) then
         print *,'Error -- qstore:'
         print '(a,5g12.5)','qmax = ',qmax
         print '(a,5g12.5)','qst  = ',icestate(i,j)%ice%qstore
         print '(a,5g12.5)','diff = ',(qmax-icestate(i,j)%ice%qstore)
      end if


!KAL  if (i0+i==39 .and. j0+j==69) then
!KAL     !KAL
!KAL     open(33,file='icetst.dat',form='formatted',position='append')
!KAL     write(33,'(f14.3,15e16.7)') dtime,
!KAL &      iceXXX(3)%fice,
!KAL &      iceXXX(3)%hice, 
!KAL &      iceXXX(3)%hice-ice_init(3)%hice, 
!KAL &      qtotice(3),
!KAL &      condb(3),
!KAL &      qml_ice,
!KAL &      freeze_temp(sml) - freeze_temp(icest%sml),
!KAL &      tmlXXX-icest%tml,
!KAL &      tmlLLL-icest%tml,
!KAL &      tml   -icest%tml,
!KAL &      icest%tml   -freeze_temp(icest%sml),
!KAL &      tml   -freeze_temp(sml),
!KAL &      qfac(3),
!KAL &      Isurflx(i,j),
!KAL &      Isswflx(i,j)
!KAL     close(33)
!KAL     open(33,file='icetst2.dat',form='formatted',position='append')
!KAL     write(33,'(f14.3,5e15.5)') dtime,
!KAL &     iceXXX%fice
!KAL     close(33)
!KAL  end if


      end do
      end do
      end do


      if (errsnw.or.errflx) then
         call xcstop('(errsnw .or.errflx)')
         stop '(errsnw .or.errflx)'
      end if

      end subroutine icestate_thermf
      end module m_icestate_thermf

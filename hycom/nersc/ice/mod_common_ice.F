      module mod_common_ice
      use mod_xc
      implicit none

! --- common blocks for the ice and snow part of the model
! --- KAL: Taken from v 2.2.37. Allow for "RELO" (dynamic allocation)
! --- KAL: Removed common blocks
! --- TODO: Remove CPP options, replace with fortran ifs

      real, save ::
     &     albi_m,albi,albs,albs_m,albw,albw_d, albi_d,                 &
     &     rhoice,rhosnw,rhowat,emiss,gasconst,rkice,hocond,            &
     &     fusi,fuss,fice_max,tice_m,tsnw_m,hice_min,epsmol,sice,       &
     &     rksnw,cpair,cpsw,stefanb,aice,bice,awater,bwater,t0deg       &
     &    ,albsa(12)

      real, dimension(:,:), save, allocatable ::                        &
     &  ficem,                                                          &
     &  hicem,                                                          &
     &  tsrfm,                                                          &!LB Temperature of the surface (snow/ice/water)
     &  hsnwm,                                                          &
     &  ticem,                                                          &!LB Temperature of the ice top surface 
     &  qbrine,                                                         &
     &  cawdir,                                                         &
     &  pemnp,                                                          &
     &  clat,                                                           &
     &  radfl0,                                                         &
     &  tauxice,                                                        &
     &  tauyice,                                                        &
     &  qfrz,                                                           &
     &  iceu,                                                           &
     &  icev,                                                           &
     &  delta_icevol,                                                   &
     &  delta_snwvol,                                                   &
     &  surf_albedo_sum,                                                &
     &  surf_qsw_sum                                                   

      integer, dimension(:,:), save, allocatable ::                     &
     &  MIZ_MASK

#if defined(ICE_DYN_DIAG)
      real, dimension(:,:), save, allocatable ::                        &
     &  stressp,                                                        &
     &  stressm,                                                        &
     &  stress12,                                                       &
     &  pice,                                                           &
     &  strainI,                                                        &
     &  strainII

#endif

#if defined (ICE_NEST) && defined (ICE)
      real, dimension(:,:,:), save, allocatable ::                      &
     &  uicenest,                                                       &
     &  vicenest,                                                       &
     &  hicenest,                                                       &
     &  ficenest,
     &  hsnwnest
#endif
#if defined (TEST_ICE_AGE)
      real, dimension(:,:), save, allocatable ::                        &
     &  fy_frac ,                                                       &
     &  fy_age,                                                         &
     &  rdg_frac
#endif
#if defined (PARAM_EST)
      real, dimension(:,:), save, allocatable ::                        &
     &  msshb ,                                                         &
     &  sstb
#endif
#if defined (ALBSNW_EVOL)
      real albs_min, albs_max, hsnw_lim
      real, dimension(:,:), save, allocatable ::                        &
     & albsnwm,                                                         &
     & surf_fscov  
      real taua,                                                        &
     & tauf,                                                            &
     & swe_newalb                                                       
#endif      
#if defined (SSNOWD_ICE)
      real, dimension(:,:), save, allocatable ::                        &
     &   hmelt                                                           ! Cumulated melt depth (m)
      real, dimension(:,:), save, allocatable ::                        &
     &   hprcp                                                           ! Cumulated precipitations depth (m)
      real cv_snw                                                        ! coefficient of variation : depends on the surface.
                                                                         ! Based on SHEBA dataset for sea ice.  
#endif    

      !KAL 20160715 New variables in stand alone
      real, dimension(:,:), save, allocatable ::                        &
     &  corio, pang, pasp, depths,                                      &
     &  util1, util2, util3,                                            &
     &  plon,plat,                                                      &
     &  ulon,ulat,                                                      &
     &  vlon,vlat,                                                      &
     &  qlon,qlat,                                                      &
     &  scpx,scpy,                                                      &
     &  scux,scuy,                                                      &
     &  scvx,scvy,                                                      &
     &  scu2,scv2,                                                      &
     &  scp2,scq2,                                                      &
     &  scp2i,scq2i,                                                    &
     &  scuxi,scvyi,                                                    &
     &  aspux, aspuy,                                                   &
     &  aspvx, aspvy,                                                   &
     &  scqx,scqy
      integer, save :: mapflg
      integer, save :: itest=-1
      integer, save :: jtest=-1
      real, parameter :: epsil  = 1.0d-11
      real, parameter :: hugel  = 2.0**100



      character(len=*), parameter :: flnmgrd="regional.grid"
      character(len=*), parameter :: flnmdep="regional.depth"


      contains

         
      subroutine icedat
      implicit none
      albi_m   =.60      ! --                 max albedo over melting ice
      albi_d   =.73      ! --                 max albedo over dry ice
      albs     =.80      ! --                 albedo over snow
      albs_m   =.71      ! --                 albedo over melting snow
      albw_d   =.065     ! --                 albedo over water; diff light
      rhoice   =900.     ! kg / m^3           density of ice
      rhosnw   =330.     ! kg / m^3           density of snow
      rhowat   =1000.    ! kg / m^3           density of pure water 
      emiss    =.97      ! --                 emissivity of water
      gasconst =.287e3   ! pa m^3 / (k kg)    gas constant
      rkice    =2.04     ! w / (m k)          ice conductivity
      hocond   =2.5e6    ! j / kg             heat of condensation/vap
      fusi     =3.02e8   ! j / m^3            heat of fusion of ice
      fuss     =1.10e8   ! j / m^3            heat of fusion of snow
      fice_max =.995     ! --                 maximum fractional ice cover
      tice_m   =273.05   ! k                  melting point of ice
      tsnw_m   =273.15   ! k                  melting point of snow
      hice_min =.1       ! m                  minimum ice thickness
      epsmol   =.622     ! --                 molecular weight of h2o/dry air
      !albsa    =12*.75 
      albsa    =albs     ! Seasonal albedo -- set to fixed value albs above 

      sice     = 6.      ! per mil            salinity of seaice
      rksnw    =.31      ! w / (m k)          snow conductivity
      cpair    =1004.    ! j / (k kg)         specific heat of dry air
      cpsw     =3987.    ! j / (k kg)         specific heat of seawater
      stefanb  =5.67e-8  ! w / (m^2 k^4)      stefan-boltzman constant
      aice     =9.5      ! --                 vapor pressure parameters
      bice     =7.66     ! k                    ..
      awater   =7.5      ! --                   ..
      bwater   =35.86    ! k                    ..
      t0deg    =273.15   ! k                  zero deg celcius in k
#if defined (ALBSNW_EVOL)
      albs_min =.71      ! --                 minimum snow albedo
      albs_max =.85      ! --                 maximum snow albedo 
      taua     =.008     ! --                 constant for linear decrease of albedo for dry snow
      tauf     =.24      ! --                 constant for exponential decrease of wet snow due to wet metamorphism
      swe_newalb =0.002  ! m                  snowfall water equivalent depth necessary to refresh albedo back to maximum value
      hsnw_lim =.02      ! m                  limit snow depth to compute snow cover fraction (CICE formulation)
#endif
#if defined (SSNOWD_ICE)
      cv_snw   =.68      ! --                 coefficient of variation of the lognormal snow depth distribution 
                         !                    It depends on the ice type and age. Based on SHEBA dataset for sea ice. 
#endif       

      end subroutine icedat
C
C
C
      subroutine iceinit(sss,sst)
      implicit none
      real, intent(in) :: sst(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) !KAL 20160715
      real, intent(in) :: sss(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) !KAL 20160715
      integer :: i,j
      real    :: flagi1,flagi2,tfrz, flagi, flaghemi, margin
      !KAL - include 'common_blocks.h'
      allocate(ficem          (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(hicem          (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(tsrfm          (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(hsnwm          (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(ticem          (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(qbrine         (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(cawdir         (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(pemnp          (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(clat           (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(radfl0         (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(tauxice        (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(tauyice        (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(qfrz           (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(iceu           (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(icev           (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(delta_icevol   (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(delta_snwvol   (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(surf_albedo_sum(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(surf_qsw_sum   (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))

      allocate(MIZ_MASK       (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))

#if defined(ICE_DYN_DIAG)
      allocate(stressp  (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(stressm  (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(stress12 (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(pice     (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(strainI  (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(strainII (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(MIZ_MASK (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
#endif

#if defined (ICE_NEST) && defined (ICE)
      allocate(uicenest (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2))
      allocate(vicenest (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2))
      allocate(hicenest (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2))
      allocate(ficenest (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2))
      allocate(hsnwnest (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2))
#endif

#if defined (TEST_ICE_AGE)
      allocate(fy_frac (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(fy_age  (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(rdg_frac(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
#endif

#if defined (PARAM_EST)
      allocate(msshb (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(sstb (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
#endif

#if defined (ALBSNW_EVOL)
      allocate(albsnwm   (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(surf_fscov(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
#endif      

#if defined (SSNOWD_ICE)
      allocate(hmelt(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(hprcp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
#endif


      !KAL 20160715 New variables for stand-alone setup
      allocate(scpx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(scpy(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(scux(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(scuy(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(scvx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(scvy(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(scqx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(scqy(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(plon(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(plat(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(ulon(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(ulat(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(vlon(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(vlat(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(qlon(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(qlat(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
c
      allocate(scp2 (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(scq2 (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(scu2 (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(scv2 (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(scp2i(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(scq2i(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(scuxi(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(scvyi(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(aspux(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(aspvy(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(aspvx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(aspvy(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(scqx (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(scqy (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(util1 (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(util2 (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(util3 (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(depths(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(corio (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(pang  (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(pasp  (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))

      margin=nbdy


      if (mnproc==1) then
         write(lp,*) 'Calling iceinit '
      end if
!$OMP PARALLEL DO PRIVATE(j,i,flagi,flaghemi,tfrz) 
!$OMP&      SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         !tfrz=273.216-.057*saln(i,j,1,1)
         tfrz=273.216-.057*sss(i,j)
         tfrz=tfrz-t0deg
c
         ! 1 if temperature < freezing point + .3
         ! 0 otherwise
         !flagi=(1.+sign(1.,tfrz+.3 -temp(i,j,1,1)))/2.
         flagi=(1.+sign(1.,tfrz+.3 -sst(i,j)))/2.
c
         ! 0  in -50 < lat < 50
         flagi=flagi*(1.+sign(1.,abs(plat(i,j))-50.))/2.
c
         ! Hemisphere flag, 1 if plat > 0
         flaghemi=(1.+sign(1.,plat(i,j)))/2.
c
         ! Initial ice thickness - different for N/S hemispheres
         hicem(i,j)=flagi*(flaghemi*2.5 + (1.-flaghemi)*1.0)
c
         ficem(i,j)=flagi*.95
         hsnwm(i,j)=0.
         !KAL 20160715 ticem(i,j)=temp(i,j,1,1)-flagi*4.+273.15
         ticem(i,j)=sst(i,j)-flagi*4.+273.15
         tsrfm(i,j)=ticem(i,j)
         qfrz(i,j)=0.0
         pemnp(i,j)=0.
#if defined (TEST_ICE_AGE)
         fy_age  (i,j)=0.
         fy_frac (i,j)=0.
         rdg_frac(i,j)=0.
#endif
#if defined (ALBSNW_EVOL) 
         albsnwm(i,j)=(albs_min+albs_max)*.5
#endif
#if defined(SSNOWD_ICE)
        hmelt(i,j)=0. 
        hprcp(i,j)=hsnwm(i,j)
#endif         
      enddo
      enddo
!$OMP END PARALLEL DO
c

      tauxice=0.
      tauyice=0.

      MIZ_MASK = 0!!test variable

      end subroutine iceinit
c

      subroutine clat_turb
c
c --- -------------------------------------------------------------------
c --- read numerical values of the latent transfer coefficient based
c --- on the tabelled values in isemer et al, j clim., p. 1180, 1989
c ---  n o t e :
c --- i-index gives the wind at 10m height
c ---            from 0 to 30 m/s in intervals of 2 m/s
c --- j-index gives the virtual air-sea temperature difference at 10m height
c ---            from -7 to +7 deg c in intervals of .5 deg c
c --- for all but the equatorial and sub-equatorial waters, virtual temp is
c --- close to real temp (see gill, 1982, p. 41)
c --- -------------------------------------------------------------------
c


      use mod_xc
      implicit none

      REAL d_wind,d_temp,clat1(16,29)
      integer skip
      integer i,j
      logical :: ex
c
c      dimension clat1(16,29)
c
c --- read data file
      inquire(exist=ex,file='iwh_tabulated.dat')
      if (.not.ex) then
         if (mnproc==1) 
     &   write(lp,*) 'clat_turb: iwh_tabulated.dat does not exist'
         call xcstop ('(clat_turb)')
         stop '(clat_turb)'
      end if
      open(19,file='iwh_tabulated.dat',STATUS= 'OLD')
      if (mnproc==1) print *,'reading wh_tabulated.dat'

      do 100 j=1,29       
      do 100 i=1,16
         read(19,*) d_wind,d_temp,clat1(i,j)
         clat1(i,j)=clat1(i,j)*1.e-3
  100 continue
      close(19)
c
c --- let  i=1  represents wind speeds in the interval [0, 2) m/s
c ---      i=2  represents wind speeds in the interval [2, 4) m/s ...
c
      do 111 j=1,29
      do 110 i=1,15
  110 clat(i,j)=(clat1(i,j)+clat1(i+1,j))*.5
  111 clat(16,j)=clat(15,j)
c
c
c --- let  j=1  represents temp differences in the interval [-7, -6.5) deg c
c ---      j=2  represents temp differences in the interval [-6.5, -6) deg c...
c
      do 121 i=1,16
      do 120 j=1,28
  120 clat(i,j)=(clat1(i,j)+clat1(i,j+1))*.5
  121 clat(i,29)=clat(i,28)
c
      end subroutine clat_turb
c
c
c
      real function icevolume()
      implicit none
      !KAL - include 'common_blocks.h'

      integer i,j,l
      real,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: vol
      real*8 :: icevol8

      vol=0.0
!$OMP PARALLEL DO PRIVATE(j,l,i) 
!$OMP&SCHEDULE(STATIC,jblk)
      do j=1,jj
      do l=1,isp(j)
      do i=ifp(j,l),ilp(j,l)
         vol(i,j)=ficem(i,j)*hicem(i,j)*scpx(i,j)*scpy(i,j)
      end do
      end do
      end do
!$OMP END PARALLEL DO

      !Gather tot ice volume from tiles
      call xcsum(icevol8,vol,ip)
      icevolume=icevol8
      end function icevolume
c
c
      real function iceextent()
      implicit none
      !KAL - include 'common_blocks.h'
      integer i,j,l
      real extent(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      real*8 :: iceext8
      extent=0.0
!$OMP PARALLEL DO PRIVATE(j,l,i) 
!$OMP&SCHEDULE(STATIC,jblk)
      do 100 j=1,jj
      do 100 l=1,isp(j)
      do 100 i=ifp(j,l),ilp(j,l)
         extent(i,j)=ficem(i,j)*scpx(i,j)*scpy(i,j)
100   continue
!$OMP END PARALLEL DO
c --- Gather ice extent from tiles to get tot extent
      call xcsum(iceext8,extent,ip)
      iceextent=iceext8
      end function iceextent



CKAL c --- Ice volume/ice area diagnostics
CKAL      subroutine diag_icevol(dtime)
CKAL      !KAL use mod_year_info, only: daysinyear
CKAL      !KAL use mod_hycom_nersc
CKAL#if defined(ICESTATE)
CKAL      use mod_icestate , only : icestate
CKAL      use mod_icestate_diag
CKAL#endif
CKAL      implicit none
CKAL      real*8,           intent(in) :: dtime
CKAL      integer :: iyear, iday, ihour, diy 
CKAL      integer*8 :: ihour2,isec
CKAL      real    :: icevol,iceext
CKAL      real, dimension(itdm,jtdm) :: modlat
CKAL      !KAL include 'common_blocks.h'
CKALc
CKAL      call forday(dtime, yrflag, iyear,iday,ihour)
CKAL      call xcaget(modlat,plat,0)
CKAL      ihour2=nint(dtime*86400.d0,kind=8)/3600
CKAL      isec =nint(dtime*86400.d0,kind=8)-ihour2*3600
CKAL      diy=daysinyear(iyear,yrflag)
CKAL      if (mod(ihour,6)==0 .and. isec < baclin) then
CKAL#if defined (ICESTATE) 
CKAL         icevol=icestate_volume()
CKAL         iceext=icestate_extent(.15)
CKAL#else
CKAL         icevol=icevolume()
CKAL         iceext=iceextent()
CKAL#endif
CKAL         if (mnproc==1) then
CKAL            open(10,file=rungen//'icevolume.dat',position='append')
CKAL            write(10,'(i5,g15.9,3g15.9)') 
CKAL     &         imem, iyear+(iday+(ihour+isec/3600.d0)/24.d0)/diy,
CKAL     &         dtime,icevol*1.0d-9,iceext*1.0d-12
CKAL            close(10)
CKAL         end if
CKAL         if (icevol < 1.0d12 .and. maxval(modlat) > 85) then 
CKAL            call xcstop('(diag_icevol:Arctic Ice disappeared)')
CKAL         end if
CKAL      end if
CKAL      end subroutine diag_icevol


#if defined (ICESTATE) 
! ===================================================================================
! =========================== icestate2ice ==========================================
! ===================================================================================
! Routine transfers values from icestate to 'common_ice.h' variables. 
! For diagnostic purposes (HYCOM prtsol)
      subroutine icestate2ice
      use mod_icestate  , only:icestate
      implicit none
      integer i,j,l, margin
      margin=nbdy
C$OMP PARALLEL DO PRIVATE(j,l,i) 
C$OMP&SCHEDULE(STATIC,1)
      do j=1-margin,jj+margin
      do l=1,isp(j)
      do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
         ficem(i,j) = sum(icestate(i,j)%ice%fice)
         hicem(i,j) = sum(icestate(i,j)%ice%fice*icestate(i,j)%ice%hice)
         hsnwm(i,j) = sum(icestate(i,j)%ice%fice*icestate(i,j)%ice%hsnw)
         ticem(i,j) = sum(icestate(i,j)%ice%fice*icestate(i,j)%ice%tsrf)
         ticem(i,j) = ticem(i,j) / (ficem(i,j) + 1e-6)
         hicem(i,j) = hicem(i,j) / (ficem(i,j) + 1e-6)
         hsnwm(i,j) = hsnwm(i,j) / (ficem(i,j) + 1e-6)
         tsrfm(i,j) = ticem(i,j)*ficem(i,j) + (1.-ficem(i,j))*
     &      icestate(i,j)%tml
      end do
      end do
      end do
C$OMP END PARALLEL DO
      end subroutine
#endif

c   Ice advection routine
c      
c
      
c     Ice advection routine
c      
c
      subroutine iceadv(h,u,v,scuy,scvx,scp2i,scp2,dt)
!
! --- ------------------------------------------------------------------
! --- Advection is done with flux limited 3rd order WENO in space and
! --- 2nd order Runge-Kutta in time
! --- ------------------------------------------------------------------
!
      use mod_xc
!
      implicit none
!
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::               &
     &  h,u,v,scuy,scvx,scp2i,scp2
      real dt
!
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: sao,hp
      real dtm
      integer i,j,l,margin

       sao(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)=0
!
! --- Use a modified time step since velocities are in m/s while scale
! --- factors are in cm
!      dtm=dt*1.e2
       dtm=dt
! --- Prediction step
      call weno3pd(h,sao,u,v,scuy,scvx,scp2i,scp2,dtm)
      margin=nbdy 
!$OMP PARALLEL DO
      do j=1-margin,jj+margin
        do l=1,isp(j)
        do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
          hp(i,j)=h(i,j)+dtm*sao(i,j)
        enddo
        enddo
      enddo
!$OMP END PARALLEL DO
!
! --- Correction step
      call weno3pd(hp,sao,u,v,scuy,scvx,scp2i,scp2,dtm)
!$OMP PARALLEL DO
      do j=1-margin,jj+margin
        do l=1,isp(j)
        do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
          h(i,j)=.5*(h(i,j)+hp(i,j)+dtm*sao(i,j))
        enddo
        enddo
      enddo
!$OMP END PARALLEL DO
      return
      end subroutine iceadv





       subroutine weno3pd(g,sao,u,v,scuy,scvx,scp2i,scp2,dt)
!
! --- ------------------------------------------------------------------
! --- By a weighted essentially non-oscillatory scheme with up to 3th
! --- order accuracy, obtain the spatial advective operator of a
! --- 2-dimensional field defined at the scalar points of a C-grid. The
! --- fluxes are limited to make the scheme positive definite.
! --- Advective velocities in the i- and j-direction are defined at u-
! --- and v-points, respectively.
! --- ------------------------------------------------------------------
!
      use mod_xc
!
      implicit none
!
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::               &
     &  g,sao,u,v,scuy,scvx,scp2i,scp2
      real dt
!
      real cq00,cq01,cq10,cq11,ca0,ca1,eps
      parameter (cq00=-1./2.,cq01= 3./2.,                               &
     &           cq10= 1./2.,cq11= 1./2.,                               &
     &           ca0=1./3.,ca1=2./3.,                                   &
     &           eps=1.e-12)
!
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::               &
     &  ful,fuh,fvl,fvh,gt
      real q0,q1,a0,a1,q
      integer i,j,l,im1,im2,ip1,jm1,jm2,jp1
!
! --- Compute grid cell boundary fluxes. Split in a low order flux
! --- (donor cell) and a high order correction flux.
!
!$OMP PARALLEL DO
      do j=0,jj+2
        do i=0,ii+2
          ful(i,j)=0.
          fuh(i,j)=0.
          fvl(i,j)=0.
          fvh(i,j)=0.
        enddo
      enddo
!$OMP END PARALLEL DO
!
!      call xctilr(g, 1,1, 3,3, halo_ps)
!
!$OMP PARALLEL DO PRIVATE(im1,im2,q0,q1,a0,a1,ip1)
      do j=0,jj+1
        do l=1,isu(j)
        do i=max(0,ifu(j,l)),min(ii+2,ilu(j,l))
c       do i=1,ii
          im1=i-1
!
          if (u(i,j).gt.0.) then
            im2=im1-iu(im1,j)
!
            q0=cq00*g(im2,j)+cq01*g(im1,j)
            q1=cq10*g(im1,j)+cq11*g(i  ,j)
!
            a0=ca0
            a1=ca1*(abs(g(im2,j)-g(im1,j))+eps)                         &
     &            /(abs(g(im1,j)-g(i  ,j))+eps)
!
            ful(i,j)=u(i,j)*g(im1,j)*scuy(i,j)
!

          else
            ip1=i+iu(i+1,j)
!
            q0=cq11*g(im1,j)+cq10*g(i  ,j)
            q1=cq01*g(i  ,j)+cq00*g(ip1,j)
!
            a0=ca1
            a1=ca0*(abs(g(im1,j)-g(i  ,j))+eps)                         &
     &            /(abs(g(i  ,j)-g(ip1,j))+eps)
!
            ful(i,j)=u(i,j)*g(i  ,j)*scuy(i,j)
!
          endif
!
          fuh(i,j)=u(i,j)*(a0*q0+a1*q1)/(a0+a1)*scuy(i,j)-ful(i,j)
!
        enddo
        enddo
      enddo
!$OMP END PARALLEL DO
!
!$OMP PARALLEL DO PRIVATE(jm1,q0,q1,a0,a1,jm2,jp1)
      do j=0,jj+2
        jm1=j-1
        do l=1,isv(j)
        do i=max(0,ifv(j,l)),min(ii+1,ilv(j,l))
c       do i=1,ii
!
          if (v(i,j).gt.0.) then
            jm2=jm1-iv(i,jm1)
!
            q0=cq00*g(i,jm2)+cq01*g(i,jm1)
            q1=cq10*g(i,jm1)+cq11*g(i,j  )
!
            a0=ca0
            a1=ca1*(abs(g(i,jm2)-g(i,jm1))+eps)                         &
     &            /(abs(g(i,jm1)-g(i,j  ))+eps)
!
            fvl(i,j)=v(i,j)*g(i,jm1)*scvx(i,j)
!
          else
            jp1=j+iv(i,j+1)
!
            q0=cq11*g(i,jm1)+cq10*g(i,j  )
            q1=cq01*g(i,j  )+cq00*g(i,jp1)
!
            a0=ca1
            a1=ca0*(abs(g(i,jm1)-g(i,j  ))+eps)                         &
     &            /(abs(g(i,j  )-g(i,jp1))+eps)
!
            fvl(i,j)=v(i,j)*g(i,j  )*scvx(i,j)
!
          endif
!
          fvh(i,j)=v(i,j)*(a0*q0+a1*q1)/(a0+a1)*scvx(i,j)-fvl(i,j)
!        
        enddo
        enddo
      enddo
!$OMP END PARALLEL DO
!
! --- Update field with low order fluxes.
!$OMP PARALLEL DO
      do j=0,jj+1
        do l=1,isp(j)
        do i=max(0,ifp(j,l)),min(ii+1,ilp(j,l))
!        do i=1,ii
          gt(i,j)=g(i,j)-dt*(ful(i+1,j)-ful(i,j)                        &
     &                      +fvl(i,j+1)-fvl(i,j))*scp2i(i,j)
        enddo
        enddo
      enddo
!$OMP END PARALLEL DO
!
! --- Obtain fluxes with limited high order correction fluxes.
      q=.25/dt
!$OMP PARALLEL DO
      do j=1,jj
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii+1,ilu(j,l))
c       do i=1,ii
          fuh(i,j)=ful(i,j)+max(-q*gt(i  ,j)*scp2(i  ,j),               &
     &                      min( q*gt(i-1,j)*scp2(i-1,j),fuh(i,j)))
        enddo
        enddo
      enddo
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
      do j=1,jj+1
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
c       do i=1,ii
          fvh(i,j)=fvl(i,j)+max(-q*gt(i,j  )*scp2(i,j  ),               &
     &                      min( q*gt(i,j-1)*scp2(i,j-1),fvh(i,j)))
        enddo
        enddo
      enddo
!$OMP END PARALLEL DO
!
! --- Compute the spatial advective operator.
!$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c       do i=1,ii
          sao(i,j)=-(fuh(i+1,j)-fuh(i,j)+fvh(i,j+1)-fvh(i,j))*scp2i(i,j)
        enddo
        enddo
      enddo
!$OMP END PARALLEL DO
!
      return
      end subroutine weno3pd

#if defined(SSNOWD_ICE) 

      subroutine fpsolver(cv,D_a,D_ave,D_m)

      use m_get_erfc

c This uses a "fixed-point" iteration algorithm to solve for the
c   depth of melt when the cell-average snow depth is known.

      implicit none

      integer maxiter,i
      real cv,D_a,D_m,sca_Dm,zeta,D_ave,xlambda,z_Dm,tol,D_m_old

c Define the melt-depth tolerance to be 0.1 mm.
      tol = 1.0e-4
      maxiter = 10

c Set the initial guess to a small number (this does not seem to
c   affect the iterations required to gain convergence, and
c   starting with a large number can lead to divergence of the
c   solution).
      D_m_old = 1.0e-9

      zeta = sqrt(log(1.0 + max(cv,0.001)**2))
      xlambda = log(D_a) - 0.5 * zeta**2

      do i=1,maxiter
        z_Dm = (log(D_m_old) - xlambda) / zeta
        sca_Dm = 0.5 * get_erfc(z_Dm/sqrt(2.0))
        D_m = (0.5 * exp(xlambda + 0.5*zeta**2) *
     &    get_erfc((z_Dm - zeta)/sqrt(2.0)) - D_ave) / sca_Dm

c       print *, i,cv,D_m,sca_Dm,D_a,D_ave

        if (abs(D_m - D_m_old).lt.tol) return
        D_m_old = D_m

      return
      end subroutine fpsolver
#endif

      end module mod_common_ice

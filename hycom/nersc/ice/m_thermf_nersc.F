      module m_thermf_nersc
      contains
      subroutine thermf_nersc(m,n,rt)
c
c === ==================================================================
c --- derive heat- and salt-fluxes for input to mixing routines
c --- net       heat flux surflx  in  w / m^2
c --- shortwave heat flux sswflx  in  w / m^2 (component of net flux)
c --- net salt flux salflx  in  10-3 kg / (m^2 sec)
c ---  ***flx > 0  means that the ml is varmed/made saltier 
c ---
c --- This routine partly replaces thermf in standard hycom - the
c --- exception is surface and wall relaxation which is handled by
c --- thermf
c === ==================================================================
c === ==================================================================
c --- computation of the thermodynamic forcing over open water, sea ice 
c --- and snow-covered sea ice
c ---
c --- this routine is based on si-units, and all of the temperatures are
c --- expressed in kelvin, with one exception;
c ---             temp(i,j,1)  is in deg celcius    
c ---
c --- authors:  helge drange and knud simonsen, nersc
c --- version:  1.0
c --- date:     sep 17, 1997
c ---
c --- documentation:
c --- drange, h, and k. simonsen, 1997: formulation of air-sea fluxes in 
c ---          the esop2 version of micom, tech. rep. no. 125, 
c ---          nansen environmental and remote sensing center, 
c ---          edv. griegsv. 3a, n-5037 solheimsviken, norway
c ---
c --- Modyfied in several places by Knut Liseter knut.liseter@nersc.no
c === ==================================================================
c
      use mod_xc
      use mod_common_ice
      use mod_forcing_nersc
      use mod_year_info, only: year_info
#if defined (SSNOWD_ICE)
      use m_get_erfc
#endif 

      implicit none
      include 'common_blocks.h'
      include 'stmt_fns.h'
      integer,         intent(in) :: m,n
      type(year_info), intent(in) :: rt
c
      real :: qturb,qlw,lqfrac,smlinit,tmlinit,hmlinit,dpinit
      real, parameter :: seadns=1.025E-3
      real, parameter :: geps=0.1
      real, parameter :: snwth=0.40 ! Snow thickness threshold
c --  due to new implicit none  -- safer 
      integer :: i,j,l
      real :: rhosw, rhoswi, cpfac, tsrf, cpfaci, srfmlt, tice, 
     &        fice, hice, coalbs, coalbi, rhoair, fice_i, vpair_i, 
     &        vpair_w, fice_m, qair_i, qair_w, rice, fice_fac,
     &        qml, vpml, fice_ , vpsrf, qsrf, vpimlt, qimlt, 
     &        ficexx, fice_f, qswi0, dtemp, qsww, qsws, omfice,
     &        tice_fdegc, fachice, fqlat1,fqlat,
     &        hice_, clatent, fqlw, fqlwcc, csensible, fqsens,
     &        prcp, cc, rh, hsnw, fqlww1, fqlwi2, fqlwi1,
     &        fqlww2, fcondi, facice, winda, tair, hml_, 
     &        hmlxx, wind, hicexx, omfachice, evap, prec, snwfac, 
     &        emnp, f1, f2, hml, tmlxx, tml_, fconds, fac1, 
     &        omfice_fac, tsrf1, hsnwsfm, tsrf_s, tice1, qfac, 
     &        qrest_s, dvsnw, omfacsnw, facsnw, qtot1, omfacice, 
     &        omfacmlt, omsrfm, facmlt, tsrf_i, qice,
     &        tml, sml, sml_, tice_f, smlxx, qrest_i, f3,
     &        f1mf2mf3, qtotw, qcorr, dvice, omhifac,
     &        hifac, tsrfxx, f_snw_ice, qbot, qtotn, qtotp,
     &        albsnw, rhos_rhoi, h2o_snw,
     &                    dti, dt, rzero,         emp_fac, 
     &        spr, windi, vol_ice, tml_f, freeze, omfacf, 
     &        facf, fice1, hml_f, sml_f, tsrf_, f1mf2, fush,
     &        dfice, ficef, hice1, fac, fussi, hice2, vol_f,
     &        omqfac, qrest, smn, dsgdt, ficei, tmn, tmps,
     &        tmp, dsgds, ffac, omffac, tml_m, cawdif, 
     &        fusii,  fheat, f4, omfac, tmpt, hice_m,
     &        sml_m, hml_m, swfl,  gfac, h_eff,
     &        fyfac, delta_fice, finalfy
#if defined (ALBSNW_EVOL) || defined(SSNOWD_ICE)
     &        ,fsnwcov  
#endif
#if defined (ALBSNW_EVOL)
     &        ,israin  
#endif    
#if defined (SSNOWD_ICE)  
     &        ,hmlt,hpcp,xi,lam,zz,snet,hsnowfall
#endif      
      real*8 :: timeqsw0
c --- 'cormn4' = 4 times minimum coriolis magnitude
      real, parameter ::   cormn4=4.0e-5  ! corio(4N) is about 1.e-5
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::  
     &   ficem_old, hicem_old, hsnwm_old
     




c --- ------------------------------------------------------------------
c --- definition of inline functions
c --- ------------------------------------------------------------------

      real :: step, xx1, xx2
      real :: vapp, ax,bx,tx
      real :: humid,px,ex
      real :: qtot,qswx, fqsensx,fqlatx,qsrfx,qairx,qlw1x,qlw2x,
     &        fcondx,tfreezx,tsrfx,tairx
      real :: tnew

      integer :: iturb, jturb
      real    :: xtemp,xwind

c --- define step function;  step=1  if  xx1.ge.xx2 ; =0  if xx1.lt.xx2
      step(xx1,xx2)=(1.+sign(1.,xx1-xx2))*.5
c 
c --- define vapour pressure ( vapp  has unit pa;  tx  is temp in k)
      vapp(ax,bx,tx)=611.*10.**(ax*(tx-273.16)/(tx-bx))

c --- define specific humidity (dimensionless)
      humid(px,ex)=.622*ex/(px-.378*ex)

c --- define total surface flux at the marine boundary layer (w/m^2)
      qtot(qswx,fqsensx,fqlatx,qsrfx,qairx,qlw1x,qlw2x,fcondx,
     .     tfreezx,tsrfx,tairx)=
     .     qswx
     .    -fqsensx*(tsrfx-tairx)
     .    -fqlatx *max(rzero,qsrfx-qairx)
     .    -qlw1x
     .    -qlw2x*tsrfx
     .    +fcondx*(tfreezx-tsrfx)

c --- define new surface temperature at the marine boundary layer (k)
      tnew(qswx,fqsensx,fqlatx,qsrfx,qairx,qlw1x,qlw2x,fcondx,
     .     tfreezx,tairx)=
     .    (qswx
     .    +fqsensx*tairx
     .    -fqlatx *max(rzero,qsrfx-qairx)
     .    -qlw1x
     .    +fcondx*tfreezx)/
     .    (fqsensx
     .    +qlw2x
     .    +fcondx)

c --- define appropriate i-index of turbulent transfer coefficients
      iturb(xwind)=min0(16,int(xwind*.5)+1)
c
c --- define appropriate j-index of turbulent transfer coefficients
      jturb(xtemp)=max0(1,min0(29,int((xtemp+7.)*2.)+1))
c
c --- Follow "thermo flag" set by blkdat
      if (.not.thermo) then
         surflx(:,:)=0.0
         salflx(:,:)=0.0
         sswflx(:,:)=0.0
         hekman(:,:)=0.0
         print *,'Warning: Surface thermal forcing is off!'
         return
      endif
c
c --- TODO make sure thermf_nersc+fluxes plays nice (at least partially) when 
c --- flxflg.ne. 99
c
c --- ------------------------------------------------------------------
c --- set some quantities
c --- ------------------------------------------------------------------

      dt      =baclin
      dti     =1./dt
      spr     =1.e5      ! surface pressure (pa)   1 atm.
      rzero   =0.        ! real zero

      emp_fac =1./(seadns*spcifh*onecm*thref)

      cawdif=1.-albw_d

      fusii =1./fusi
      fussi =1./fuss

      h2o_snw=1000./rhosnw
      rhos_rhoi=rhosnw/rhoice

      windi=0.006               !m/s  friction velocity below ice; corresponds
                                !      to a wind vel. of ~4 m/s over open water
c --- ------------------------------------------------------------------
c --- compute the 24-hrs mean downward irradiance field (w/m^2) once a day
c --- ------------------------------------------------------------------
c
c --- TODO - adhere to hycom time
      ! NB: hycom uses a reference date from 1901 (!) for "real" years.
      ! We use data from the "rt" variable as input to qsw now, in order
      ! to avoid a shift in the seasons...
      timeqsw0=rt%idd + rt%iss/86400. + rt%ihh/24. + dt/86400.
      !print *,'Thermf_nersc:time :',time, timeqsw0
c --- Radian is now pi/180 , before it was 180/pi
#if defined (DIURNAL)
c --- Incoming shortwave radiation is updated every 3 hours
c --- Coarse reproduction of the diurnal cycle
      if (mod(rt%ihh,3)==0.and.rt%iss<baclin) then
        call qsw0(radfl0,cawdir,timeqsw0,rt%daysinyear,rt%ihh)
#if defined (CALANUS)
        call calc_radmax(plat/radian,timeqsw0,rt%daysinyear)
#endif 
      end if
#else
c --- Daily average incoming shortwave radiation      
      if (mod(nstep-1,int(86400/dt))==0) then
        call qsw0(radfl0,cawdir,timeqsw0,rt%daysinyear)
#if defined (CALANUS)
        call calc_radmax(plat/radian,timeqsw0,rt%daysinyear)
#endif 
      end if
#endif      
      if (nstep-nstep1.eq.1) then
        call clat_turb
#if defined (DIURNAL)
        call qsw0(radfl0,cawdir,timeqsw0,rt%daysinyear,rt%ihh)
#else
        call qsw0(radfl0,cawdir,timeqsw0,rt%daysinyear)
#endif
#if defined (CALANUS)
        call calc_radmax(plat/radian,timeqsw0,rt%daysinyear)
#endif  
      endif
c
c --- ------------------------------------------------------------------
c --- Seasonally dependent snow albedo
c --- ------------------------------------------------------------------
c
#if defined (ALBSNW_EVOL)
c ---  Snow albedo is a prognostic variable
#else
      albsnw=albsa(l0)*w0+albsa(l1)*w1+albsa(l2)*w2+albsa(l3)*w3
#endif
      margin=nbdy
c
c
c --- ------------------------------------------------------------------
c --- start of main loop
c --- ------------------------------------------------------------------
C$OMP   PARALLEL DO PRIVATE(
C$OMP&  j, l, i, sml, sml_, smlXX, tice_f, tml, tml_, tmlXX,
C$OMP&  hml, hmlXX, hml_, tair, wind, cc, rh, prcp, rice,
C$OMP&  fice_fac, omfice_fac, tice_fdegc, fachice, omfachice, hice,
C$OMP&  hsnw, hiceXX, hice_, fice, omfice, ficeXX, fice_, fice_f,
C$OMP&  fice_m, tice, tsrf, rhosw, rhoswi, cpfac, cpfaci, srfmlt,
C$OMP&  omsrfm, albi, coalbi, coalbs, rhoair, vpair_i, vpair_w, qair_i,
C$OMP&  qair_w, vpml, qml, vpsrf, qsrf, vpimlt, qimlt, qsww, qsws,
C$OMP&  qswi0, dtemp, clatent, fqlat1, fqlat, csensible, fqsens, fqlw,
C$OMP&  fqlwcc, fqlwi1, fqlwi2, fqlww1, fqlww2, winda, facice, fcondi,
C$OMP&  fac1, fconds, snwfac, prec, evap, f1, f2, emnp, tsrf1,qtot1,
C$OMP&  facsnw, omfacsnw, facmlt, omfacmlt, omfacice, tice1, tsrf_s,
C$OMP&  hsnwsfm, dvsnw, qrest_s, qfac, qbot, qice, tsrf_i, qrest_i,
C$OMP&  hifac, omhifac, f_snw_ice, tsrfXX, dvice, f3, f1mf2mf3, qcorr,
C$OMP&  qtotw, qtotp, qtotn, omqfac, qrest, freeze, tml_f, vol_ice,
C$OMP&  fice1, facf, omfacf, hice2, vol_f, fac, tsrf_, f1mf2, sml_f,
C$OMP&  hml_f, ficef, hice1, fush, dfice, ffac, omffac, tml_m, fheat,
C$OMP&  hice_m, sml_m, hml_m, f4, omfac, ficei,         tmn, smn,
C$OMP&  dsgdt, dsgds, qlw, qturb, lqfrac, tmpt,tmps, tmp, hmlinit,
C$OMP&  dpinit,smlinit,tmlinit,swfl,gfac,h_eff,fyfac,
C$OMP&  delta_fice,finalfy)
C$OMP&  SCHEDULE(STATIC,jblk)
      do 100 j=1-margin,jj+margin
      do 100 l=1,isp(j)
      do 100 i=max(ifp(j,l),1-margin),min(ilp(j,l),ii+margin)

CKAL - Fix this when new restart setup is ready (for re-producability)
      !if (nint(delt1/baclin)==2) then  ! Leapfrog
      !   !if (i==thermi.and.j==thermj) print *,'Leapfrog'
      !   tmlinit   =(temp(i,j,1,n)+temp(i,j,1,m))/2. 
      !   smlinit   =(saln(i,j,1,n) + saln(i,j,1,m))/2.
      !   hmlinit   =.5*(dp(i,j,1,n)+dp(i,j,1,m))/onem          
      !else  ! Euler
      !   !if (i==thermi.and.j==thermj) print *,'Euler'
      !   smlinit   =saln(i,j,1,n)
      !   hmlinit   =dp(i,j,1,n)/onem          
      !   tmlinit   =temp(i,j,1,n)
      !endif
      tmlinit   =(temp(i,j,1,n)+temp(i,j,1,m))/2. 
      smlinit   =(saln(i,j,1,n) + saln(i,j,1,m))/2.
      hmlinit   =.5*(dp(i,j,1,n)+dp(i,j,1,m))/onem          

      sml = smlinit
      sml_  =sml
      smlXX  =sml

c --- freezing point of sea water (in k)

      tice_f=273.216-.057*sml
      
      tml  = tmlinit+t0deg
      tml_  =tml
      tmlXX  =tml

      hml = hmlinit
      dpinit  = hmlinit*onem
      hmlXX  =hml
      hml_  =hml

c --- ------------------------------------------------------------------
c --- interpolate the climatological forcing fields ( radfl  from  -qsw0-)
c  .    tair             k       air temperature     
c  .    wind             m/s     wind speed
c  .    cc               -       fractional cloud cover
c  .    rh               -       relativ humidity
c  .    prcp   (> 0)     m/s     precipitation
c     ( radfl  (> 0)     w/m^2   short wave downward irradiance )
c --- ------------------------------------------------------------------
c
c
c ---   airt = air temperature (C)
        tair=airtmp(i,j,l0)*w0+airtmp(i,j,l1)*w1
     &      +airtmp(i,j,l2)*w2+airtmp(i,j,l3)*w3
c ---   prcp = precipitation (m/sec; positive into ocean)
        prcp=precip(i,j,l0)*w0+precip(i,j,l1)*w1
     &      +precip(i,j,l2)*w2+precip(i,j,l3)*w3
c ---   relhum = relhum ()
        rh  =relhum(i,j,l0)*w0+relhum(i,j,l1)*w1
     &      +relhum(i,j,l2)*w2+relhum(i,j,l3)*w3
c ---   wind = windspeed ()
        wind=wndspd(i,j,l0)*w0+wndspd(i,j,l1)*w1
     &      +wndspd(i,j,l2)*w2+wndspd(i,j,l3)*w3
c ---   clouds 
        cc  =clouds(i,j,l0)*w0+clouds(i,j,l1)*w1
     &      +clouds(i,j,l2)*w2+clouds(i,j,l3)*w3
c --- Air temp to Kelvin
      tair=tair+t0deg
c --- No ice relaxation
      rice=0.0
c --- ------------------------------------------------------------------
c --- adjust  ssts  such that  ssts  equlas freezing temperature minus
c --- 4 degrees (arbitrary choice) if ice is present according to observations 
c --- with concentration higher than 30%
c ---
c --- ml surface temperatures below freezing temperature is set to
c --- freezing temp in _mxlayr_ and _tsadvc_
c --- ------------------------------------------------------------------

C     if (trelax .or.srelax) then
C
C        ssts=max(ssts,tice_fdegc+1.e-8)*fice_fac     !no ice
C    .    +   (tice_fdegc-4.)      *omfice_fac     !ice
C
C --- turn off ice-relaxation if  hice>3m
C        fice_fac=(1.+sign(1.,3.-hicem(i,j)))*5.
C        ssts=ssts*fice_fac+tmlinit*(1.-fice_fac)
C        ssts = ssts2
C     end if

      fice_fac=(1.+sign(1.,.3-rice))*.5
      omfice_fac=1.-fice_fac
      tice_fdegc=tice_f-t0deg

c --- ------------------------------------------------------------------
c --- continue initialization
c --- ------------------------------------------------------------------
c
c --- KAL: Keep to compute growth
      ficem_old(i,j)=ficem(i,j)
      hicem_old(i,j)=hicem(i,j)
      hsnwm_old(i,j)=hsnwm(i,j)
c
      fachice  =step(hicem(i,j),1.e-4)
      omfachice=1.-fachice
c
      hice  =hicem(i,j)*fachice
      hsnw  =hsnwm(i,j)*fachice
c
      hiceXX =hice
      hice_ =0.
c
      fice  =ficem(i,j)*fachice
      omfice=1.-fice
c
      ficeXX =fice
      fice_ =omfice
      fice_f=0.
      fice_m=0.
c
      tice  =min(tice_m,ticem(i,j))*fachice+tml*omfachice
C   tsrf  =min(tice,tsrfm(i,j))*fachice+tml*omfachice
C    Vincent V allows surface temp to be equal to the value for melting snow
C    273.15 K
      tsrf  =tsrfm(i,j)*fachice+tml*omfachice
#if defined (ALBSNW_EVOL) 
      albsnw=albsnwm(i,j)*fachice
#endif      
#if defined (SSNOWD_ICE)
      !if(.not.(hmelt(i,j).ge.0)) then
      !   write(*,*) 'hmelt NaN av th',i0+i,j0+j
      !   write(*,*) 'hmelt NaN av th',hmelt(i,j)
      !endif
      hmlt=hmelt(i,j)*fachice
      hpcp=hprcp(i,j)*fachice
#endif 
c
c --- ------------------------------------------------------------------
c --- definition of various quantities
c --- ------------------------------------------------------------------

      rhosw =(1000.+sig(tml-t0deg,sml))             ! kg-sw / m^3
      rhoswi=1./rhosw

      cpfac =cpsw*rhosw*hml                          ! w / (m^2 k)
      cpfaci=1./cpfac                                ! m^2 k / w

cdiag sal0w=sml*hml*rhosw
cdiag sal0i=sice*hice*rhoice*fice
cdiag sal0 =sal0w+sal0i
cdiag wat0= hml*rhosw+(hice*rhoice+hsnw*rhosnw)*fice

c --- ------------------------------------------------------------------
c --- compute ice and snow albedo based on the previous time step.
c --- melting is defined to occur if  (tsrf > tice_m)
c --- ------------------------------------------------------------------
      srfmlt=step(tsrf,tice_m)
      omsrfm=1.-srfmlt

      albi  =min(albi_d,.08+.44*hice**.28)
      coalbi=1.-(albi*omsrfm+min(albi_m,albi)*srfmlt)


#if defined (SSNOWD_ICE)
c --- Snow cover fraction is computed       
      xi=(log(1+cv_snw**2))**(0.5)
c --- No snow on the ground at the previous time step
      if(hpcp .lt. 1e-5) then
         fsnwcov=0.
c --- Snow on the ground at the previous time step
      else
c --- Melt has occured
      if(hmlt.gt.0.) then
             lam=log(hpcp)-0.5*xi**2
             zz=(log(hmlt)-lam)/xi
c --- Compute the snow cover fraction under the assumption of the erf distribution curve
             fsnwcov=0.5*get_erfc(zz/2**0.5)
c --- the ground is covered by snow
          else 
             fsnwcov=1.
          endif
      endif
      surf_fscov(i,j) = fsnwcov
#endif 

#if defined (ALBSNW_EVOL)
#if defined (SSNOWD_ICE)
c --- we use SSNOWD formulation for the snow cover fraction
#else      
c --- snow cover fraction (CICE model formulation)
      fsnwcov=hsnw/(hsnw+hsnw_lim)
      surf_fscov(i,j) = fsnwcov
#endif       
      coalbs=1.-fsnwcov*albsnw-(1.-fsnwcov)*(1.-coalbi)
#else
      coalbs=1.-(albsnw*omsrfm+albs_m*srfmlt)
#endif

c --- ------------------------------------------------------------------
c  -  compute various quantities that occur in the heat flux expressions
c --- ------------------------------------------------------------------

c --- density of air (kg / m^3)
c  .    spr         pa                  surface pressure (1 atm ~ 1.e5 pa)
c  .    gasconst    pa m^3 / (k kg)     gas constant
      rhoair=spr/(gasconst*tair)

c --- vapour pressure (pa) for  temp=temp_air  over (ice,snow) and water,
c --- and the corresponding specific humidities
      vpair_i=vapp(aice  ,bice  ,tair)
      vpair_w=vapp(awater,bwater,tair)
c 
      qair_i =rh*humid(spr,vpair_i)
      qair_w =rh*humid(spr,vpair_w)

c --- vapour pressure (pa) for  temp=tml  (over water),
c --- and the corresponding specific humidity
      vpml=vapp(awater,bwater,tml)
      qml =humid(spr,vpml)

c --- vapour pressure (pa) for  temp=tsrf over (ice,snow),
c --- and the corresponding specific humidity
      vpsrf=vapp(aice,bice,tsrf)
      qsrf =humid(spr,vpsrf)

c --- vapour pressure (pa) for  temp=tice_m over (ice,snow),
c --- and the corresponding specific humidity
      vpimlt=vapp(aice,bice,tice_m)
      qimlt =humid(spr,vpimlt)

c --- vapour pressure (pa) for  temp=tice_f  over (ice,snow),
c --- and the corresponding specific humidities
c...  vptf_i=vapp(aice,bice,tice_f)
c...  qtf_i =humid(spr,vptf_i)

c --- net short wave downward irradiance (w / m^2) (i=ice, s=snow, w=water)
c  .    coalbx      -                   co-albedox = (1 - albedox)

      qsww =radfl0(i,j)*(cawdir(i,j)*(1.-cc)+cawdif*cc)
      qsws =radfl0(i,j)*coalbs
      qswi0=radfl0(i,j)*coalbi*.932    !nb: 932

      facsnw  =step(hsnw,1.e-4)

      ! KAL - added variable to keep surface albedo and sw flux
      surf_qsw_sum   (i,j) = fice*(facsnw*qsws + (1.-facsnw)*qswi0) +
     &   (1.-fice)*qsww
      surf_albedo_sum(i,j) = 
     &   fice*(facsnw*(1.-coalbs) + (1.-facsnw)*(1.-coalbi)) +
     &   (1.-fice)*(1.-(cawdir(i,j)*(1.-cc)+cawdif*cc))

c932  qswi =qswi0*.83
c932  qswi0=qswi0*.17


c --- upward latent heat flux factor (w / m^2)
c  .    clatent     -                   latent heat transfer coefficient 
c                                       (isemer, 1989)
      dtemp=tair-tsrf
      !print *,i,j,wind,wndspd(i,j,:)
      clatent=clat(iturb(wind),jturb(dtemp))
      fqlat1=rhoair*clatent*wind
      fqlat =fqlat1*hocond

c --- upward sensible heat flux factor (w / ( m^2 k))
c  .    cpair       J / (k kg)          specific heat of dry air
c  .    csensible   -                   sensible heat transfer coefficient 
c                                       (isemer, 1989)
      csensible=0.94*clatent
      fqsens=rhoair*cpair*csensible*wind

c --- upward net long wave radiation flux factor (w / (m^2 k))
c  .    emiss       -                   emissivity
c  .    stefanb     w / (m^2 k^4)       stefan-boltzman constant
c  .    rlat        rad                 latitude
      fqlw  =emiss*stefanb*tair**3
cH<2.1fqlwcc=1.-(.5+.246*abs(plat(i,j)/radian))*cc**1.2
      fqlwcc=1.-(.5+.246*abs(plat(i,j)*radian))*cc**1.2
c
      fqlwi1=fqlw*tair*((.254-4.95e-5*vpair_i)*fqlwcc-4.)
      fqlwi2=fqlw*4.
c
      fqlww1=fqlw*tair*((.254-4.95e-5*vpair_w)*fqlwcc-4.)
      fqlww2=fqlwi2
c
c --- ------------------------------------------------------------------
c --- wind generated tke that enters the ml is modified in the presence
c --- of ice ( winda  and  windi  have dim  m/s)
c --- ------------------------------------------------------------------
      winda=sqrt(cd*rhoair*rhoswi)*wind
      ustar(i,j)=winda*omfice+windi*fice

c --- ------------------------------------------------------------------
c --- conductive factors for ice and snow
c --- ------------------------------------------------------------------
      facice=step(hice,1.e-4)
      fcondi=rkice*facice/(hice+(1.-facice))

      fac1=rksnw*hice+rkice*hsnw
      fconds=rkice*rksnw*facice/(fac1+(1.-facice))

c === ==================================================================
c --- KAL - Correction for unresolved ice thickness distribution.
c --- From Fichefet and Morales Maqueda (1997)
      h_eff = (hsnw*rkice + hice*rksnw) / (rkice+rksnw)
      h_eff = max(h_eff,.5*exp(1.)*geps)
      gfac  = 2*h_eff/(geps*exp(1.))
      gfac  = 1 + 0.5*log(gfac)
      fcondi=fcondi*gfac
      fconds=fconds*gfac
      
c === ==================================================================
c --- compute evaporation minus precipitation fluxes over open water and ice
c ---  o  evaporation    only over open water
c ---  o  precipitation  over open water and ice/snow
c ---
c ---  prec                         is  in m/month
c ---  evap  and  emnp              are in m/s
c ---  pemnp                        is  in pa/s   (used in -barotp-)
c ---  emnp > 0  means that there is a net transport of water
c ---            from the ocean to the atmosphere
c === ==================================================================

c --- prec in the form of snow (snwfac=1) if surface temp is below melting,
c --- ice exist, and tair is below melting
      snwfac=step(tice_m,tsrf)*step(fice,1.e-4)*step(tsnw_m,tair)         

c --- precipitation and evaporation in  m/s
      prec =prcp ! Hycom reads in m/s now (from 2.1)
      evap =fqlat1*max(0.,qml-qair_w)*rhoswi

c --- update lead ml salinity and thickness
      f1=hml_*rhosw
      f2=(evap-prec)*dt*1000.
      sml_=f1*sml_/(f1-f2)
      hml_=(f1-f2)/rhosw

c --- update ice fraction ml salinity and thickness
      f1=hmlXX*rhosw
      f2=-prec*dt*1000.*(1.-snwfac)
      smlXX=f1*smlXX/(f1-f2)
      hmlXX=(f1-f2)/rhosw

c --- Update snow thickness
#if defined(SSNOWD_ICE)
c --- Snow fall is included later
      hsnowfall=dt*prec*h2o_snw*snwfac
#else
      hsnw=hsnw+dt*prec*h2o_snw*snwfac
#endif

c --- compute pemnp
      emnp=(evap *omfice
     .     -prec*(omfice+fice*(1.-snwfac)))

      pemnp(i,j)=emnp*onem
      
#if defined (ALBSNW_EVOL)
c --- Albedo evolution based on Douville et al (95): includes temporal
c     decrease of the snow albedo due to snow aging. New snowfall
c     refresh the snow albedo
c
c --- israin=1 if precipitation are falling in the form of rain      
      israin=step(prec,1e-15)*(1.-snwfac)
c --- Albedo decreases with time
      albsnw=albsnw-taua*dt/86400.     
c --- If rain, albedo decreases exponentially :
      albsnw=(1.-israin)*albsnw+israin*((albsnw-albs_min)*              &
     & exp(-tauf*dt/86400.)+albs_min)
c --- If snow is falling : albedo is refresh :
      albsnw=albsnw+snwfac*dt*prec*(albs_max-albsnw)/swe_newalb
c --- If melting occurs and rain is not falling
      albsnw=(1.-israin)*(1.-srfmlt)*albsnw+israin*albsnw+(1.-israin)*  &
     &   srfmlt*((albsnw-albs_min)*exp(-tauf*dt/86400.)+albs_min)
c --- Albedo is restricted
      albsnw=min(max(albsnw,albs_min),albs_max)
#endif

cdiag wat0=wat0-(evap*omfice-prec)*dt*1000.

c === ==================================================================
c --- do the thermodynamics for the initial snow and ice covered fraction 
c --- of the grid cell. final ice-quantities labelled by 'XX' 
c === ==================================================================

c --- -------------------------------------------------------------------
c --- --->  s n o w
c --- -------------------------------------------------------------------

c --- compute the snow surface temperature (k),
      tsrf1=tnew(qsws,fqsens,fqlat,qsrf,qair_i,fqlwi1,fqlwi2,fconds,
     .           tice_f,tair)

c --- ... and the heat budget for melting snow (w/m^2)
      qtot1=qtot(qsws,fqsens,fqlat,qimlt,qair_i,fqlwi1,fqlwi2,fconds,
     .           tice,tsnw_m,tair)

c --- set some switches
      facsnw  =step(hsnw,1.e-4)
      omfacsnw=1.-facsnw

      facmlt  =step(tsrf1,tsnw_m)
      omfacmlt=1.-facmlt

      omfacice=1.-facice

c --- compute the ice-snow interface temperature (k)
      tice1=(rksnw*hice*tsrf+rkice*hsnw*tice_f)
     .      *facice/(fac1+(1.-facice))
      tice =facice*tice1 + omfacice*tsrf

c --- set the snow surface temperature
      tsrf_s=facsnw*(tsnw_m*facmlt+tsrf1*omfacmlt) + omfacsnw*tsrf

c --- update the snow thickness
      hsnwsfm=qtot1*facsnw*facmlt*dt*fussi     !surface melt
#if defined(SSNOWD_ICE)
c --- Compute net snow gain or loss for this time step
      snet=hsnowfall-hsnwsfm
c --- Avoid accumulate extremely small snow fall
      if(snet.gt.0.0 .and. snet.lt.1.0e-6) snet=0.0
c --- Define constant xi
      xi=(log(1+cv_snw**2))**(0.5)
c
c --- No Snow at the previous time step    
      if(hpcp .eq. 0.) then
c ---Accumulation
         if(snet.gt.0.) then
            hpcp=snet
            hmlt=0.
            hsnw=hpcp
         else
            hpcp=0.
            hmlt=0.
         endif
c --- Snow on ground at previous time step
      else 
c --- Accumulation 
         if(snet.gt.0.) then
c --- 100% snow cover at previous time step. No melting has occured
            if(hmlt .eq. 0.)then
               hpcp=hpcp+snet
               hsnw=hpcp
c --- Less than 100% snow cover at previous time step, with new accumulation.
c     Two options are possible (described in Liston , 2004)
c     We use option 2 : new accumulation decreases the melt depth, 
c     "pushing" the depletion curve back towards 100 % 
            else
                lam=log(hpcp)-0.5*xi**2
                zz=(log(hmlt)-lam)/xi
                hsnw=0.5*hpcp*get_erfc((zz-xi)/2**0.5)-                 &
     & hmlt*0.5*get_erfc(zz/2**0.5)
c --- Two different cases :
c       - there is enough new accumulation to push the melt depth to zero
c       - there is not enough (the melt depth is still positive) : hsnw is updated
c     and hmlt is reduced based on Eq 18 from Liston, 2004. (iterative solution) 
               if(hsnw+snet.ge.hpcp) then
                  hpcp=hsnw+snet
                  hmlt=0.
                  hsnw=hpcp
               else
                  hsnw=hsnw+snet
                  call fpsolver(cv_snw,hpcp,hsnw,hmlt)
               endif
           endif

c --- Melt 
         elseif(snet.lt.0.) then
            hmlt=hmlt-snet
            lam=log(hpcp)-0.5*xi**2
            zz=(log(hmlt)-lam)/xi
c --- Compute the average snow depth in the grid cell
            hsnw=0.5*hpcp*get_erfc((zz-xi)/2**0.5)-                     &
     & hmlt*0.5*get_erfc(zz/2**0.5) 
c --- When snow cover becomes too thin, reset melt and acccumulation variable to 0
            if(0.5*get_erfc(zz/2**0.5) .lt. 0.005) then
              hsnw=0.
              hpcp=0.
              hmlt=0.
            endif
c --- No change since the last time step
         else
         endif
      endif
      dvsnw  =max(0.,dvsnw-hsnw)               !change in snow volume
!      
c --- Check NaN apparition. Uncomment for debugging      
      !if(.not.(hmlt.ge.0)) then
      !   write(*,*) 'hmelt NaN ap up',i0+i,j0+j
      !   write(*,*) 'hmelt NaN ap up',hmelt(i,j),hmlt,'voix', voix
      !   write(*,*) 'hmelt NaN ap up',hprcp(i,j),hpcp
      !   write(*,*) 'hmelt NaN ap up',hsnw,snet
      !   write(*,*) 'hmelt NaN ap up',lam,xi,'fscov',fsnwcov
      !endif
#else
      dvsnw  =min(hsnw,hsnwsfm)                !change in snow volume
      hsnw   =hsnw-hsnwsfm                     !update snow thickness
#endif       

c --- if  hsnw<0  , let the excess heat (w/m^2) go to melt ice
      qrest_s=-min(rzero,hsnw)*fuss*dti
      hsnw=max(rzero,hsnw)

c --- set  tsrf_s  to  tice_m  if  (qrest_s>0)
      qfac=step(qrest_s,1.e-5)
      tsrf_s=qfac*tice_m + (1.-qfac)*tsrf_s

c --- -------------------------------------------------------------------
c --- --->  i c e  
c --- -------------------------------------------------------------------

c --- compute the ice surface temperature (k),
      tsrf1=tnew(qswi0,fqsens,fqlat,qsrf,qair_i,fqlwi1,fqlwi2,fcondi,
     .           tice_f,tair)

c --- and the corresponding switch for melting
      facmlt  =step(tsrf1,tice_m)
      omfacmlt=1.-facmlt

c --- compute the heat budget for the ice layer in case of melting ice,...
      qtot1=qtot(qswi0,fqsens,fqlat,qimlt,qair_i,fqlwi1,fqlwi2,fcondi,
     .           tice_f,tice_m,tair)

c --- ... the heat flux at the bottom of the ice (w m^2),
c ---  qbot>0  denote melting,...
c932  flight=exp(-1.5*hice)*omfacsnw
      qbot=(tml-tice_f)*cpfac*dti+fcondi*(tice-tice_f)
c932 .    +qswi0*flight
      qbot=min(qbot,(hice*fusi+hsnw*fuss)*dti)

c ---- ... and the total heat flux applied to the ice layer
      qice=qrest_s+qtot1*omfacsnw*facmlt+qbot

c --- update the brine heat reservoir
c932  qbrine(i,j)=qbrine(i,j)+qswi0*(omfacsnw-flight)

c --- if  tsrf1<tmlt_i , extract heat from  qbrine  so that 
c ---  tsrf1  is slightly above  tice_m

c --- set the new surface temperature
      tsrf_i=facsnw*tsrf_s+omfacsnw*(facmlt*tice_m+omfacmlt*tsrf1)

c --- update the ice thickness
      hiceXX=hice-qice*dt*fusii

c --- if  hiceXX<0  , let the excess heat (w/m^2) go to warm ml water
      qrest_i=-min(rzero,hiceXX)*fusi*dti

      hiceXX=max(rzero,hiceXX)

      hifac  =step(hiceXX,1.e-7)
      omhifac=1.-hifac

c --- -------------------------------------------------------------------
c --- convert snow to ice if snow load is larger than the updrift of ice
c --- -------------------------------------------------------------------

c -1- if the snow load is larger than the updrift of ice, then the snow is
c -1- converted to ice. optionally this may also happen when air temperature
c -1- is above melting temp, which implies that precep is in liquid form.
c -1- note that this procedure equals a salt source, since the salinity of the
c -1- ice is non-zero, while the salinity of the snow is zero
c 
c     f_snw_ice=step(hiceXX,hsnw*rhos_rhoi)            !=0 ==>  conversion
c    .         *step(tice_m,tair)    
c     hiceXX    =hiceXX + (1.-f_snw_ice)*rhos_rhoi*hsnw
c     hsnw     =hsnw  *f_snw_ice                      !new snow thickness

c -2- if the snow load is larger than the updrift of ice, then the snow is
c -2- sweeped into the ocean (hice is not changed, heat changes are negl)

      f_snw_ice=step(hiceXX,hsnw*rhos_rhoi)
c     .          *step(tsnw_m,tair)
      dvsnw=dvsnw+(1.-f_snw_ice)*hsnw                 !update change in snow vol
      hsnw     =hsnw*f_snw_ice                        !new snow thickness
#if defined(SSNOWD_ICE)
      hpcp     =hpcp*f_snw_ice
      hmlt     =hmlt*f_snw_ice
#endif      

c --- -------------------------------------------------------------------
c --- adjustment of surface albedo
c --- -------------------------------------------------------------------

c -1- The snow albedo is decreased gradually to the melting surface albedo,
c -1- during the summer in order to mimic the melt pond effect in case
c -1- the thermo-snow model is not able to obtain melting surfaces.
c -1- This is cheating, but required due to shortcomming of the snow model.
c    
c -1- snow is converted in one time step
c -1-
c      f_snw_ice= step(albsnw,albs_m+.1)   
c      hiceXX    =hiceXX + (1.-f_snw_ice)*rhos_rhoi*hsnw
c      hsnw     =hsnw  *f_snw_ice
c 
c -2- the difference between the snow albedo and the melting surface albedo
c -2- is used to calculate the fraction of snow, which is converted   into ice
c -2- in the late summer. since the snow albedo albedo is changing gradually,
c -2- this proceudre provides a smooth conversion, in contrast to an abrupt
c -2- conversion over one time step.
c -2-
c     f_snw_ice=10.*max(albs_m+.1-albsnw,rzero)       !weight E[0;1]
c                                                     !if weight is pos. and
c                                                     !snow thickness less than 
c                                                     !1 cm, set weight=1 
c     f_snw_ice=max(f_snw_ice,step(f_snw_ice,1.e-4)*step(0.01,hsnw))
c     hiceXX    =hiceXX +f_snw_ice*hsnw*rhos_rhoi
c     hsnw     =hsnw *(1.-f_snw_ice)
C
c -3- snow thicker than 20 cm is assumed to melt and enter the ocean through 
c -3- rifts and leads in the ice (season dependent albedo)
C
C     f_snw_ice=min(1.,5.*max(albs_m+.2-albsnw,rzero)) !weight E[0;1]
C     f_snw_ice=f_snw_ice*step(hsnw,0.20)
C     hsnw     =hsnw *(1.-f_snw_ice)
C     dvsnw=dvsnw+f_snw_ice*hsnw
C
CKAL  Bug discovered by VV.  The above approach (3) is troublesome at best. We set
C     the limit to  snwth (snow threshold, typically 0.40 m) here
      dvsnw=max(0.,hsnw-snwth) ! Snow thickness > snow threshold
      hsnw =min(hsnw,snwth)
#if defined(SSNOWD_ICE)
C --- When the threshold is reached : ideal snow distribution with a
C     snow cover fraction equal to 1  and an accumulated snow depth
C     equal to the threshold depth     
      if(dvsnw .gt. 0.) then
         hpcp = snwth
         hmlt = 0.
      end if      
      !if(.not.(hmlt.ge.0)) then
      !   write(*,*) 'hmelt NaN ap thi',i0+i,j0+j
      !   write(*,*) 'hmelt NaN ap thi',hmelt(i,j),hmlt,dvsnw
      !endif
#endif      
c --- -------------------------------------------------------------------
c --- --->  m l   w a t e r
c --- -------------------------------------------------------------------

c --- warm the ml water if  (qrest_i>0)
      tmlXX=tice_f+qrest_i*dt*cpfaci

c --- update various quantities
      tsrfXX=hifac*tsrf_i + omhifac*tmlXX

      tice =hifac*tice   + omhifac*tmlXX

      ficeXX=hifac*fice

      hsnw =hsnw*step(ficeXX,1.e-6)

c --- change in ice volume
      dvice=hiceXX-hice

c --- updated ml salinity and ml thickness
      f1= hml  *rhosw
      f2= dvice*rhoice
      f3=-dvsnw*rhosnw
      f1mf2mf3=f1-f2-f3   

      smlXX=(f1*sml-f2*sice)/f1mf2mf3

      hmlXX=f1mf2mf3/rhosw

c === ==================================================================
c --- do the thermodynamics for the initial open water fraction of the 
c --- grid cell. final quantities labelled by '_'
c === ==================================================================

c --- if the mixed layer temperature is below the freezing temperature of
c --- seawater (may occur because of changes in the ml salinity),
c --- set  tml  equal to  tice_f , and reduce the heat flux over open
c --- water according to the correction
      qcorr=min(rzero,tml-tice_f)*cpfac*dti
      tml_=max(tml,tice_f)

c --- compute the total heat flux at the surface (w / m^2).  
c ---  (qtotw > 0)  means warming of the lead water
      qtotw=qtot(qsww,fqsens,fqlat,qml,qair_w,fqlww1,fqlww2,rzero,
     .           tml_,tml_,tair)
     .     +qcorr


      qtotp=max(rzero,qtotw)
      qtotn=min(rzero,qtotw)

      qfac  =step(qtotw,rzero)
      omqfac=1.-qfac

cKAL  Fractional value of fluxes used to warm ml.
      qlw    = - fqlww1 - fqlww2*tsrf
      qturb  = -fqlat*max(qml-qair_w,rzero) - fqsens*(tml_-tair)
      lqfrac = 1.*omfice

c -1- if  (qtotw < 0) , cool the ml
      tml_=tml_+qtotn*dt*cpfaci

c --- if  (tml_ < tice_f) , derive the (negative) heat used to freeze ice
cKAL  qrest=min(rzero,tml_-tice_f)*cpfac*dti
      qrest=min(rzero,tml_-tice_f)*cpfac*dti*omfice

      freeze=step(-qrest,rzero)

      tml_ =max(tml_,tice_f)
      tml_f=tml_

c --- determine ice volume to be frozen
      vol_ice=-qrest*dt/fusi

cKAL  Modify fractional value of fluxes used to cool ml in case of freeze
cKAL  lqfrac = lqfrac * ( 1. - abs(qrest)/(abs(qtotn*omfice)+1e-7))

c --- determine starting ice thickness to be used during freezing of ice
      hice_=hice_min

c --- updated fractional ice cover
      fice1=-qrest*dt/(fusi*hice_+1.e-10)

c --- if  (fice1+ficeXX > fice_max) , set  (fice_ = fice_max-ficeXX), and
c --- increase ice thickness up to 1 m, or until the heat budget is closed
      facf =step(fice1+ficeXX,fice_max)*freeze
      omfacf=(1.-facf)                *freeze

      fice_f=fice1*omfacf+(fice_max-ficeXX)*facf

      hice2 =-qrest*dt/(fice_f*fusi+1.e-10)
      hice_ =min(1.,hice_*omfacf+hice2*facf)

c --- actual ice volume that is frozen
      vol_f=hice_*fice_f

c --- determine additional ice to be frozen in order to close heat budget
      vol_ice=max(rzero,vol_ice-vol_f)

c --- update the brine heat reservoir for vanishing  hice_
c932  ficef=ficeXX/(ficeXX+fice_f+1.e-10)
c932  qbrine(i,j)=qbrine(i,j)*ficef

c --- if ice, set  tsrf_=tice_f , otherwise  tsrf_=tml_
      fac=step(hice_,1.e-6)
      tsrf_=tice_f*fac+tml_*(1.-fac)

c --- updated ml salinity and ml thickness
      f1=hml_*rhosw
      f2=hice_*rhoice
      f1mf2=f1-f2
      sml_f=(f1*sml_-f2*sice)/f1mf2

      hml_f=f1mf2/rhosw

c --- increase  hiceXX  if needed, and if ice is present initially, and
c --- correct for the different areas of ice covered and open water fraction
c --- of the grid cell

      ficef=step(ficeXX,0.01)
      !hice1=vol_ice*omfice*ficef/(ficeXX*fice*ficef+(1.-ficef))
      hice1=vol_ice*ficef/(ficeXX*fice*ficef+(1.-ficef))

c --- new ice thickness
      hiceXX=hiceXX+hice1

c --- update ml salinity and ml thickness accordingly
      f1= hmlXX *rhosw
      f2= hice1*rhoice
      f1mf2=f1-f2

      smlXX=(f1*smlXX-f2*sice)/f1mf2

      hmlXX=f1mf2/rhosw

c -2- if  (qtotw > 0) , let the fraction  (1 - fice)  go to warm the lead,
c --- and the fraction  fice  to melt ice laterally

c -2a warming of the lead

c --- new lead water and surface temperature 
      tml_=(tml+omfice*qtotp*dt*cpfaci)*qfac+tml_*omqfac
      tsrf_=tml_*qfac+tsrf_*omqfac


c KAL Modify lqfrac if qtotp is used to melt ice
cKAL  lqfrac = lqfrac*(1. - step(qtotp,rzero)*fice )
c  
c -2b lateral melting of ice and snow

c --- new fractional ice cover
      fush=fusi*hiceXX+fuss*hsnw
      dfice=-omfice*qtotp*dt/(fush+1.e-8)

      fice1=ficeXX+dfice

      ffac  =step(fice1,1.e-2)
      omffac=1.-ffac

c --- if  (fice1 < 0) , determine the excess heat

      qrest=-min(rzero,fice1)*fush*dti
      fice1= max(rzero,fice1)

c --- ... and warm the ml accordingly
      tml_m=tmlXX+qrest*dt*cpfaci

cKAL  Modify lqfrac if ml is re-heated
cKAL  lqfrac= lqfrac*(1 + max(dfice-ficeXX,0.)/(dfice+1e-7))

      fheat=step(qrest,rzero)
      tmlXX =tml_m*fheat+tmlXX*(1.-fheat)

c --- change in ice and snow volume
      dvice=hiceXX
      dvsnw=hsnw

c --- store change in grid-fraction of ice due to lateral melting
      fice_m=ficeXX-fice1

c --- update various quantities
      hice_m=hiceXX

c --- updated ml salinity and ml thickness
      f1=hmlXX*rhosw
      f2=-dvice*rhoice
      f3=-dvsnw*rhosnw
      f1mf2mf3=f1-f2-f3
      sml_m=(f1*smlXX-f2*sice)/f1mf2mf3

      hml_m=f1mf2mf3/rhosw

c === ==================================================================
c --- final computation of the thermodynamic variables of the grid cell
c === ==================================================================

      f1=(fice -fice_m)*hmlXX *rhosw   
      f2= fice_f       *hml_f*rhosw  
      f3=(fice_-fice_f)*hml_ *rhosw
      f4= fice_m       *hml_m*rhosw
      f1mf2=f1+f2+f3+f4

      sml=(f1*smlXX+f2*sml_f+f3*sml_+f4*sml_m)/f1mf2

      hml=(fice -fice_m)*hmlXX+fice_f*hml_f
     .   +(fice_-fice_f)*hml_+fice_m*hml_m

      fice=min(fice_max,ficeXX-fice_m+fice_f)
      omfice=1.-fice
      fac=(1.+sign(1.,fice-1.e-10))*.5
      omfac=1.-fac
      ficei=fac/(fice+omfac)

      hice=((ficeXX-fice_m)*hiceXX+fice_f*hice_)*ficei
      hsnw=hsnw*(ficeXX-fice_m)*ficei
#if defined (SSNOWD_ICE)
c --- Ice concentration is changed. To conserve snow volume, snow depth is changed 
c     accordingly. Accumulated and melt depths are adjusted. Snow cover
c     fraction is conserved
      hmlt=hmlt*(ficeXX-fice_m)*ficei
      hpcp=hpcp*(ficeXX-fice_m)*ficei
#endif      

      ffac=step(fice_f,1.e-6)
      tml=(f1*tmlXX+f2*tml_f+f3*tml_+f4*tml_m)/f1mf2
      tml=tml*(1.-ffac)+tice_f*ffac

      tsrf=((ficeXX-fice_m)*tsrfXX+fice_f*tsrf_ )*ficei+tml*omfac
      tice=((ficeXX-fice_m)*tice +fice_f*tice_f)*ficei
      tice =tice*facsnw+tsrf*omfacsnw

c --- ------------------------------------------------------------------
c --- if  hice  is below a certain min thickness, set the ice/snow 
c --- thicknesses to zero
c --- ------------------------------------------------------------------
      facice=step(hice,1.e-4)
      omfacice=1.-facice

      ficem(i,j)=fice*facice
      hicem(i,j)=hice*facice
      hsnwm(i,j)=hsnw*facice
      ticem(i,j)=tice*facice+tml*omfacice
      tsrfm(i,j)=tsrf*facice+tml*omfacice
      
#if defined (ALBSNW_EVOL)
      albsnwm(i,j)=albsnw*facice
#endif
#if defined (SSNOWD_ICE)
      hmelt(i,j)=hmlt*facice
      hprcp(i,j)=hpcp*facice
     ! if(.not.(hmelt(i,j).ge.0)) then
     !    write(*,*) 'hmelt NaN ap th',i0+i,j0+j
     !    write(*,*) 'hmelt NaN ap th',hmelt(i,j)
     ! endif
#endif 

c --- ------------------------------------------------------------------
c --- if  hice  is below a certain min thickness, set the ice/snow 
c --- thicknesses to zero
c --- ------------------------------------------------------------------

      fice_fac=step(ficem(i,j),fice_max-.01)

      hice=max(hicem(i,j),hice_min)*fice_fac
     .    +hicem(i,j)              *(1.-fice_fac)
      fice=ficem(i,j)

      ficem(i,j)=hicem(i,j)*ficem(i,j)/(hice+1.e-10)
      hicem(i,j)=hice

      fac=step(ficem(i,j),1.e-2)
      hsnwm(i,j)=hsnwm(i,j)*fice*fac/(ficem(i,j)*fac+(1.-fac))
#if defined (ALBSNW_EVOL) 
      facsnw=step(hsnwm(i,j),1.e-4)  
      albsnwm(i,j)=albsnw*facsnw
#endif      
#if defined (SSNOWD_ICE)
      hmelt(i,j)=hmelt(i,j)*fice*fac/(ficem(i,j)*fac+(1.-fac))
      hprcp(i,j)=hprcp(i,j)*fice*fac/(ficem(i,j)*fac+(1.-fac))
      !if(.not.(hmelt(i,j).ge.0)) then
      !   write(*,*) 'hmelt NaN ap th2',i0+i,j0+j
      !   write(*,*) 'hmelt NaN ap th2',hmelt(i,j)
      !   stop'Problem with NaN'
      !endif
#endif 
c === ==================================================================
c --- derive heat- and salt-fluxes that are consistent with the new
c --- temperature and salinity fields
c ---  surflx  in  w / m^2
c ---  salflx  in  gr-salt / (m^2 sec)
c ---  ***flx > 0  means that the ml is varmed/made saltier 
c --- 
c --- NOTE: The direction of salinity and heat fluxes have been changed 
c --- from MICOM.
c === ==================================================================

      tml=tml-t0deg
c
      fac=dpinit/(g*dt)
c
      surflx(i,j)=(tml-tmlinit)*spcifh*fac 
      salflx(i,j)=(sml-smlinit)       *fac
      sswflx(i,j)=qsww*lqfrac

c --- Ekman mixing depth in pressure point - also set in mixing
      hekman(i,j)=4.*cekman*ustar(i,j)/
     .(max(cormn4,
     . (abs(corio(i,j  ))+abs(corio(i+1,j  ))+
     .  abs(corio(i,j+1))+abs(corio(i+1,j+1)))))
c
      if (i==itest.and.j==jtest) then
         print '(a,2f14.5)','ficem, ficem_old :',
     &         ficem(i,j),ficem_old(i,j)
         print '(a,2f14.5)','hicem, hicem_old :',
     &         hicem(i,j),hicem_old(i,j)
         print '(a,2f14.5)','tml, tmlinit     :',tml,tmlinit
         print '(a,2f14.5)','sml, smlinit     :',sml,smlinit
         print '(a,4f14.5)','tair, wind, rh,cc:',tair,wind,rh,cc
         print '(a,3f14.5)','prcp             :',prcp
         print '(a,f14.5)' ,'salflx           :',salflx(i,j)
         print '(a,f14.5)' ,'surflx           :',surflx(i,j)
         print '(a,f14.5)' ,'sswflx           :',sswflx(i,j)
      end if
c
c --- ------------------------------------------------------------------
c --- Eqv. freshwater flux from the ML during the time step emnp (m/s)
c --- Negl. the density difference bwteen old and new time step
c --- ------------------------------------------------------------------
      emnp= salflx(i,j)*thref/saln(i,j,1,n) 
c
c --- ------------------------------------------------------------------
c --- Relaxation is moved to a relaxation routine --- see original
c --- HYCOM thermf
c --- ------------------------------------------------------------------
c
c --- ------------------------------------------------------------------
c --- KAL - Compute ice growth 
c --- ------------------------------------------------------------------
      delta_icevol(i,j)=ficem    (i,j)*hicem    (i,j) -
     &                  ficem_old(i,j)*hicem_old(i,j)
      delta_snwvol(i,j)=ficem    (i,j)*hsnwm    (i,j) -
     &                  ficem_old(i,j)*hsnwm_old(i,j)

#if defined (TEST_ICE_AGE)
c --- ------------------------------------------------------------------
c --- KAL - Update fraction and age of young, level ice
c --- ------------------------------------------------------------------

      delta_fice=ficem(i,j)-ficem_old(i,j)
      fy_age (i,j)=fy_age (i,j) + dt/86400. 

!      if (i==122.and.j==174) then
!         print *,'fy_age(i,j) -- 1',fy_age(i,j),delta_fice,fy_frac(i,j)
!      end if

      ! Ice is melting laterally
      if (delta_fice<0.) then


         ! This assumes fy and my melt at similar rates. Very simplified
         fy_frac (i,j)=max(fy_frac (i,j)+delta_fice*fy_frac (i,j),0.)
         rdg_frac(i,j)=max(rdg_frac(i,j)+delta_fice*rdg_frac(i,j),0.)

      else

         ! Update age
         finalfy=fy_frac(i,j)+delta_fice
         fyfac=step(finalfy,1e-4)

!      if (i==122.and.j==174) then
!         print *,'fy_age(i,j) -- 1.5',fyfac,fy_frac(i,j),finalfy
!      end if
         
         fy_age (i,j)=
     &     fy_age(i,j) * (fyfac*fy_frac(i,j)/(finalfy + (1.-fyfac)))
     &   + fy_age(i,j) * (1.-fyfac) 

         fy_frac(i,j)=finalfy

      end if

!      if (i==122.and.j==174) then
!         print *,'fy_age(i,j) -- 2',fy_age(i,j),delta_fice,fy_frac(i,j)
!      end if

      ! Transfer from fy to my. This is a simple step function
      fyfac=step(fy_age(i,j),365.)
      fy_frac(i,j)=(1.-fyfac) * fy_frac(i,j)
      fy_age (i,j)=(1.-fyfac) * fy_age (i,j)

c --- ------------------------------------------------------------------
c --- if  hice  is below a certain min thickness, set the FYI age 
c --- and FYI and ridged fractions to zero
c --- ------------------------------------------------------------------
      facice=step(hice,1.e-4)  
      fy_frac(i,j) = fy_frac(i,j) *facice
      fy_age(i,j)  = fy_age(i,j)  *facice
      rdg_frac(i,j) = rdg_frac(i,j) *facice

      ! More advanced - ice has a age (a) pdf linearly dependent on a in the 
      ! interval 0 to 1 year. The Following is g(1) 
      !ga_1= max(0.,6 * fy_age(i,j)/365. - 2)

      ! modify in case of negative values in linear pdf
      !if (fy_frac(i,j)/365. >.66)
      !   ga_weight = (6*fy_frac(i,j)/365. - 4.) / (12*fy_frac(i,j)/365.-6)
      !   ga_weight= 1./(1.-ga_weight)
      !else
      !   ga_weight= 1.
      !end if

      ! In time dt, this much of the ice becomes older than 365 days
      ! (time moves part of the age pdf across day 365.)
      !myloss=ga_1*ga_weight*dt/(365.*86400.)

      ! Corresponding loss of fy ice
      ! fy_frac (i,j) = fy_frac(i,j) * (1. - myloss)

      ! Ice Age is not changed by this method
#endif

c --- ------------------------------------------------------------------
c --- end of main loop
c --- ------------------------------------------------------------------

  100 continue
      !write(lp,*) 'Fix relaxation!!'

!#if defined (TEST_ICE_AGE)
!      if (mnproc==1) then
!         print *,'Test fy age, cfy'
!         print *,maxval(fy_age),maxval(fy_frac)
!      end if
!#endif
   99 format(1x,6f16.5)
      end subroutine thermf_nersc
C
C
C
      subroutine qsw0(radfl0,cawdir,qswtime,daysinyear,nhour)
c
c --- -------------------------------------------------------------------
c --- compute 24 hrs mean solar irrradiance at the marine surface layer
c --- (unit: w/m^2)
c --- -------------------------------------------------------------------
c
      use mod_xc
      use mod_forcing_nersc
      implicit none

      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), intent(out) :: 
     &  cawdir,radfl0
      real*8, intent(in) :: qswtime
      integer, intent(in) :: daysinyear
      integer, intent(in), optional :: nhour

      real pi2,deg,rad,eepsil
      integer ifrac,npart
      real fraci,absh2o,s0,day,dangle,decli,sundv
      real cc,sin2,cos2,cosz,scosz,stot,bioday,biohr,hangle
      real srad,sdir,sdif,altdeg,cfac,ssurf
      integer i,j,l
      real*8 :: daysinyear8
#if defined(DIURNAL)
      integer nh
#endif      
c
c --- Average number of days in year over a 400-year cycle (Gregorian Calendar)
      real, parameter :: daysinyear400=365.2425
      include 'common_blocks.h' ! New
c
c
      daysinyear8=daysinyear
c --- set various quantities
      pi2=8.*atan(1.)          !        2 times pi
      deg=360./pi2             !        convert from radians to degrees
      rad=pi2/360.             !        convert from degrees to radians
      eepsil=1.e-9             !        small number
      ifrac=24                 !        split each 12 hrs day into ifrac parts
      fraci=1./ifrac           !        1 over ifrac
      absh2o=0.09              ! ---    absorption of water and ozone
      s0=1365.                 ! w/m^2  solar constant
#if defined(DIURNAL)
      nh = int(nhour/3)
#endif       
c
c --- -------------------------------------------------------------------
c --- compute 24 hrs mean solar radiation at the marine surface layer 
c --- -------------------------------------------------------------------
C --- KAL: TODO - adhere to hycom time setup
      day=mod(qswtime,daysinyear8)    !0 < day < 364
      day=floor(day)
c
      dangle=pi2*day/float(daysinyear)   !day-number-angle, in radians 
      if (day<0. .or. day>daysinyear+1) then
         if (mnproc==1) then
            print *,'qsw0: Error in day for day angle'
            print *,'Day angle is ',day,daysinyear,qswtime
         end if
         call xcstop('(qsw0)')
         stop
      end if
c --- compute astronomic quantities -- 
      decli=.006918+.070257*sin(dangle)   -.399912*cos(dangle)
     $             +.000907*sin(2.*dangle)-.006758*cos(2.*dangle)
     $             +.001480*sin(3.*dangle)-.002697*cos(3.*dangle)

      sundv=1.00011+.001280*sin(dangle)   +.034221*cos(dangle)
     $             +.000077*sin(2.*dangle)+.000719*cos(2.*dangle)

c --- scan through each of the grid cells
c
!$OMP PARALLEL DO PRIVATE(j,l,i,cc,sin2,cos2,scosz,stot,bioday,
!$OMP&     npart,biohr,hangle,srad,sdir,sdif,altdeg,cfac,ssurf,
!$OMP&     cosz)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do l=1,isp(j)
      do i=max(1-margin,ifp(j,l)),min(ilp(j,l),ii+margin)

c --- compute cloudiness fraction
      cc=clouds(i,j,l0)*w0+clouds(i,j,l1)*w1
     &  +clouds(i,j,l2)*w2+clouds(i,j,l3)*w3

c --- compute astronomic quantities

      sin2=sin(plat(i,j)*radian)*sin(decli)
      cos2=cos(plat(i,j)*radian)*cos(decli)
c
c --- split each day into ifrac parts, and compute the solar radiance for 
c --- each part. by assuming symmetry of the irradiance about noon, it
c --- is sufficient to compute the irradiance for the first 12 hrs of
c --- the (24 hrs) day (mean for the first 12 hrs equals then the mean
c --- for the last 12 hrs)
c
c --- TODO - This routine can also return daily varying solar heat flux
      scosz=0.
      stot=0.
      do npart=1,ifrac
#if defined(DIURNAL)
c --- Modified shorwave radiation : 3h-average instead of 24 hour
c --- Coase reproduction of the diurnal cycle 
      if(nh==0 .or. nh==7) then
            bioday=day+(npart-.5)*fraci*.125
      else if(nh==1 .or. nh==6) then
            bioday=day+(npart-.5)*fraci*.125+0.125
      else if(nh==2 .or. nh==5) then
            bioday=day+(npart-.5)*fraci*.125+0.25
      else    
            bioday=day+(npart-.5)*fraci*.125+.375
      endif
#else                                                    
         bioday=day+(npart-.5)*fraci*.5
#endif      
         biohr=bioday*86400.                !hour of day in seconds
         biohr=amod(biohr+43200.,86400.)    !hour of day;  biohr=0  at noon
         hangle=pi2*biohr/86400.            !hour angle, in radians
c
         cosz=amax1(0.,sin2+cos2*cos(hangle)) !cosine of the zenith angle
         scosz=scosz+cosz                     !  ..accumulated..
         srad =s0*sundv*cosz                  !extraterrestrial radiation
c
!         sdir=srad*0.7**(1./(cosz+eepsil))    !direct radiation component
!         sdir=srad * exp(-0.356674943938732447/(cosz+eepsil))         
c ---    KAL prevent underflow - .7^100 = 3x10^-16 
         sdir=srad*0.7**(min(100.,1./(cosz+eepsil)))    !direct radiation component
c
         sdif=((1.-absh2o)*srad-sdir)*.5      !diffusive radiation component
         altdeg=amax1(0.,asin(min(1.0,sin2+cos2)))*deg !solar noon altitude in degrees
         cfac=(1.-0.62*cc+0.0019*altdeg)      !cloudiness correction 
         ssurf=(sdir+sdif)*cfac
         stot=stot+ssurf
      enddo
      scosz=scosz*fraci                    !24-hrs mean of  cosz
      radfl0(i,j)=stot*fraci               !24-hrs mean shortw rad in w/m^2
c
c --  Original formula was wrong ...
      !cawdir(i,j)=1.-amax1(0.15,0.05/(scosz+0.15)) 
      !cawdir(i,j)=1.-amax1(0.03,0.05/(scosz+0.15))  !Correction   - Mats
      cawdir(i,j)=1.-amin1(0.15,0.05/(scosz+0.15))   !Correction 2 - KAL
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO
c
      end subroutine qsw0
c
#if defined (CALANUS)
      subroutine calc_radmax(rlat,qswtime,daysinyear)
      use mod_xc
      use mod_forcing_nersc

      use mod_CAL06_calnusparticles, ONLY : radmax

      implicit none

      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),intent(in)::rlat
      real*8, intent(in) :: qswtime
      integer, intent(in) :: daysinyear

      real pi2,deg,rad,eepsil
      integer ifrac,npart
      real fraci,absh2o,s0,day,dangle,decli,sundv
      real cc,sin2,cos2,cosz,scosz,stot,bioday,biohr,hangle
      real srad,sdir,sdif,altdeg,cfac,ssurf
      integer i,j,l
      real*8 :: daysinyear8
      include 'common_blocks.h' ! New



c --- -------------------------------------------------------------------
c --- compute 24 hrs max solar irrradiance at the marine surface layer
c --- (unit: w/m^2)
c --- -------------------------------------------------------------------

      daysinyear8=daysinyear
c --- set various quantities
      pi2=8.*atan(1.)          !        2 times pi
      deg=360./pi2             !        convert from radians to degrees
      rad=pi2/360.             !        convert from degrees to radians
      eepsil=1.e-9             !        small number

      ifrac=24                 !        split each 12 hrs day into ifrac parts
      fraci=1./real(ifrac)     !        1 over ifrac

      absh2o=0.09              ! ---    absorption of water and ozone
      s0=1365.                 ! w/m^2  solar constant

c --- -------------------------------------------------------------------
c --- compute 24 hrs mean solar radiation at the marine surface layer 
c --- -------------------------------------------------------------------

      day=amod(qswtime,daysinyear8)    !0 < day < 364
      day=floor(day)
      dangle=pi2*day/float(daysinyear)   !day-number-angle, in radians 
      if (day<0. .or. day>daysinyear+1) then
         print *,'qsw0: Error in day for day angle'
         print *,'Day angle is ',day,daysinyear
         call xcstop('(qsw0)')
         stop
      end if

      


c --- compute astronomic quantities -- 
      decli=.006918+.070257*sin(dangle)   -.399912*cos(dangle)
     $             +.000907*sin(2.*dangle)-.006758*cos(2.*dangle)
     $             +.001480*sin(3.*dangle)-.002697*cos(3.*dangle)

      sundv=1.00011+.001280*sin(dangle)   +.034221*cos(dangle)
     $             +.000077*sin(2.*dangle)+.000719*cos(2.*dangle)

c --- scan through each of the grid cells

!    #if defined (CALANUS)
      radmax=0
!#endif

!$OMP PARALLEL DO PRIVATE(j,l,i,cc,sin2,cos2,scosz,stot,bioday,
!$OMP&     npart,biohr,hangle,srad,sdir,sdif,altdeg,cfac,ssurf,
!$OMP&     cosz)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do l=1,isp(j)
      do i=max(1-margin,ifp(j,l)),min(ilp(j,l),ii+margin)

c --- compute cloudiness fraction
      cc=clouds(i,j,l0)*w0+clouds(i,j,l1)*w1
     &  +clouds(i,j,l2)*w2+clouds(i,j,l3)*w3

c --- compute astronomic quantities

      sin2=sin(rlat(i,j))*sin(decli)
      cos2=cos(rlat(i,j))*cos(decli)

c --- split each day into ifrac parts, and compute the solar radiance for 
c --- each part. by assuming symmetry of the irradiance about noon, it
c --- is sufficient to compute the irradiance for the first 12 hrs of
c --- the (24 hrs) day (mean for the first 12 hrs equals then the mean
c --- for the last 12 hrs)

      scosz=0.
      stot=0.

      do npart=1,ifrac
         bioday=day+(npart-.5)*fraci*.5

         biohr=bioday*86400.                !hour of day in seconds
         biohr=amod(biohr+43200.,86400.)    !hour of day;  biohr=0  at noon
         hangle=pi2*biohr/86400.            !hour angle, in radians

         cosz=amax1(0.,sin2+cos2*cos(hangle)) !cosine of the zenith angle
         scosz=scosz+cosz                     !  ..accumulated..

         srad =s0*sundv*cosz                  !extraterrestrial radiation

         !print *,i,j,npart,srad,cosz,eepsil
         ! obs: .7^100 = 3x10^-16 , an already ridicolously low number ...
         sdir=srad*0.7**(min(100.,1./(cosz+eepsil)))    !direct radiation component
!         sdir=srad*0.7**(1./(cosz+eepsil))    !direct radiation component
!         sdir=srad * exp(-0.356674943938732447/(cosz+eepsil))         

         sdif=((1.-absh2o)*srad-sdir)*.5      !diffusive radiation component

         altdeg=amax1(0.,asin(min(1.0,sin2+cos2)))*deg !solar noon altitude in degrees

         cfac=(1.-0.62*cc+0.0019*altdeg)      !cloudiness correction 

         ssurf=(sdir+sdif)*cfac
         radmax(i,j)=max(ssurf,radmax(i,j))       
         stot=stot+ssurf
      enddo

      scosz=scosz*fraci                    !24-hrs mean of  cosz
!      radfl0(i,j)=stot*fraci               !24-hrs mean shortw rad in w/m^2
      
!      cawdir(i,j)=1.-amax1(0.15,0.05/(scosz+0.15)) 
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO
! CAS090106 To convert from W/m2 to uE/m2s.  
! CAS090106 At the moment this variable is not used, but it may be used later.
! CAS090106 Still not quite sure about the unit 
      radmax = radmax*4.15
      !print *,qswtime,day,dangle,decli,sundv,radfl0(63,95)
                                           !co-albedo over water for dir light
      end subroutine calc_radmax
#endif
      end module m_thermf_nersc

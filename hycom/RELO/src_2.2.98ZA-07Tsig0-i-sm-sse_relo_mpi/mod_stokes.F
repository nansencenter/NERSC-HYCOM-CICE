#if defined(ROW_LAND)
#define SEA_P .true.
#define SEA_U .true.
#define SEA_V .true.
#elif defined(ROW_ALLSEA)
#define SEA_P allip(j).or.ip(i,j).ne.0
#define SEA_U alliu(j).or.iu(i,j).ne.0
#define SEA_V alliv(j).or.iv(i,j).ne.0
#else
#define SEA_P ip(i,j).ne.0
#define SEA_U iu(i,j).ne.0
#define SEA_V iv(i,j).ne.0
#endif
      module mod_stokes
      use mod_xc  ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      use m_get_erfc
      implicit none
c
c --- HYCOM Stokes Drift from external data files
c --- this module include array declarations and data forcing
      logical, parameter, private :: debug_stokes=.true.  !usually .false.
c
c     Logical flags for Stokes Drift Effects
c     The modification of the turbulent viscosity near the surface 
c     has been implemented  following:
c          - Craig & Banner + Ardhuin-Filipot-Perenne
c          - Uchiyama-McWilliams-Shchepetkin
c
c     REAL, PARAMETER , private:: PI = 3.1415927
c      
      logical, save, public ::
     & stdflg,   ! Stokes Drift Velocities:           TRUE/FALSE
     & stdsur,   ! Stokes Drift Surface Stresses:     TRUE/FALSE (dissipation due to breaking and rolling)
     & stdbot,    ! Stokes Drift Bottom Stress:        TRUE/FALSE (dissipation due to bottom friction)
     & stdarc,   ! Stokes Drift Velocities in Archive TRUE/FALSE
     & stdtau,   ! Removing from the total wind stress the part transfered to the wave field subjected to dissipation: TRUE/FALSE
     & stdwom    ! Adding the wave-to-ocean momentum flux due to wave dissipation/breaking: TRUE/FALSE
C%%%      & stpres,   ! Stokes Wave Induced Mean Pressure: TRUE/FALSE
C%%%      & stdtau,   ! Removing from the total wind stress the part transfered to the wave field subjected to dissipation: TRUE/FALSE
C%%%      & stustr    ! if stmxsr=1,2: velocity scale for the wave vertical mixing in KPP given by WW3 (TRUE) or by Hyc

      logical, save, public :: 
     & altrlangmr ! altlng'  = STOKES: Use alternative Lang def. from Kai et al (0=F,1=T)
c
c    Constant values
c     Number of fixed interface depths for Stokes input
c
      integer,   save, public   ::
     &  nsdzi,   ! Number of fixed interface depths for Stokes input 
     &  langmr    ! Langmuir turb enhancement (KPP) 0: None 1:McWilliams-Sulliva 2:Smyth 3:McWilliams-Harcourt 4:Takaya
C%%%      & stmxsr,   ! Surface wave mixing turb enhancement 0: None, 1,2:Craig Banner Perenne et al (KPP,Yamada), 3: Uchiyama (KPP only), 4: replace diffusion by an analytical profile
C%%%      & stalfa,   ! if stmxsr=1,2: wave energy factor 0: Constant, 1:alpha from phi_oc (WW3) and u*(WW3 or Hycom) 
C%%%       real,   save, public   ::
C%%%      & al_cst,   ! if stmxsr=1,2: minimum value (stalfa=/0) or value (stalfa=0) of the wave energy factor
C%%%      & z0ctop    ! if stmxsr=1,2: minimum value (stz0tp=/0) or value (stz0tp=0) of the surface roughnes 

c     Arrays holding Surface Stokes Velocities

      real,    allocatable, dimension(:,:),
     &         save, public  ::
     & usds,    ! Surface Stokes Drift U Velocity, p-grid
     & vsds     ! Surface Stokes Drift V Velocity, p-grid

C$$ c   Arrays holding paramters for Langmuir mixing in KPP   
C$$       real,    dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
C$$      &         save, public  ::
C$$      & usdshat,     ! Surface  U hat , p-grid
C$$      & vsdshat,     ! Surface  V hat, p-grid 
C$$      & langnum     ! Langmuir number =sqrt(ustar_speed/usds_speed)

c     -U- & -V- grid Arrays holding Vertically averaged Stokes Velocities

      real,    allocatable, dimension(:,:),
     &         save, public  ::
     & usdbavg,    ! Vertical average Stokes Drift U Velocity, u-grid
     & vsdbavg     ! Vertical average Stokes Drift V Velocity, v-grid

c     Arrays holding Stokes Drift Velocities on U and V grids

      real,    allocatable, dimension(:,:,:),
     &         save, public  ::
     & usd,    ! Stokes Drift U Velocity, u-grid
     & vsd     ! Stokes Drift V Velocity, v-grid

c     Arrays holding Stokes Drift Velocities on -p- grids

      real,    allocatable, dimension(:,:,:),
     &         save, public  ::
     & usdp,    ! Stokes Drift U Velocity, p-grid
     & vsdp     ! Stokes Drift V Velocity, p-grid


c     Arrays to hold Input Stokes Drift Velocity

      real,    allocatable, dimension(:),
     &         save, private ::
     & sdzi       ! Input Stokes Drift layer depths

      real,    allocatable, dimension(:,:,:,:),
     &         save, private ::
     & usdz,    ! Stokes Drift U Velocity for fixed layers, p-grid
     & vsdz     ! Stokes Drift V Velocity for fixed layers, p-grid

c     Array holding the wave-induced mean pressure

      real,    allocatable, dimension(:,:,:),
     &         save, public  ::
     & sj     ! wave-induced mean pressure on pressure grid
c     Arrays holding the horizontal derivatives of the vertical position of the middle of the layer
c
c      real,    allocatable, dimension(:,:,:),
c     &         save, public  ::
c     & dzdx,    ! derivative of z with respect to x
c     & dzdy     ! derivative of z with respect to y

c     Arrays used to calculate the dissipations

      real,    allocatable, dimension(:,:,:),
     &         save, public  ::
     & sdka,        ! frequency as a function ot time and space (=sdk)
     & wave_dir,    ! main wave propagation direction
     & h_sig,       ! significant wave height
     & eps_brk      ! wave breaking dissipation rate

c     Dissipations (bottom and breakin)

      real,    allocatable, dimension(:,:,:),
     &         save, public  ::
     & wave_bdx(:,:,:),     ! Bottom dissipation  
     & wave_bdy(:,:,:),     ! Bottom dissipation  
     & wave_brkx(:,:,:),     ! Dissipation  by breaking
     & wave_brky(:,:,:)      ! Dissipation  by breaking
c
c     Arrays used to calculate transport

      real,    allocatable, dimension(:,:,:),
C     real,  dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2),
     &     save, private  ::
     &  tusd,      ! Stokes transport u direction, p-grid
     &  tvsd       ! Stokes transport v direction, p-grid 
      
c     Arrays used to calculate the wind stresses
 
      real,    allocatable, dimension(:,:,:),
C      real,  dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2),
     &         save, private  ::
     & tauwx,    ! Wind stress transfered to waves and subject to dissipation
     & tauwy     ! Wind stress transfered to waves and subject to dissipation

      real,    allocatable, dimension(:,:,:),
C     real,  dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2),
     &         save, private  ::
     & twomx,    ! Wave to ocean  momentum flux 
     & twomy     ! Wave to ocean momentum flux
      
C$$       real,  dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2),
C$$      &         save, public  ::
C$$      & tm01    ! Wave mean period from m0/m1
      


c     Arrays used for the parameterization of the vertical mixing

      real,    allocatable, dimension(:,:,:),
     &         save, public  ::
     & phi_ocw,  ! Turbulent kinetic energy flux (notation from WW3)
     & z0topw,   ! wind sea significant height(if stz0tp=1) or surface roughness
     & ustarw3w, ! Velocity scale from WW3
     & difx_wave ! Vertical diffusion enhancement due to wave breaking (Uchiyama et al.)
      real,    allocatable, dimension(:,:),
     &         save, public  ::
     & phi_oc,   ! Turbulent kinetic energy flux (notation from WW3)
     & alpha,    ! Wave energy factor from WW3
     & z0top,    ! Surface roughness from WW3
     & ustarw3   ! Velocity scale from WW3
c
C$$$$ c Be careful: When using radiation stresses, u must contain the eulerian and the Stokes velocity!
C$$$$       real,    allocatable, dimension(:,:,:),
C$$$$      &         save, public  ::
C$$$$      & d_sx,     ! d/dx (sxx) + d/dy (sxy) - sx. components of the radiation stress
C$$$$      & d_sy      ! d/dy (syy) + d/dx (syx) - sy. components of the radiation stress
C$$$$        
C$$$$ #endif



c    Weights for holding values during a barotropic time step.

      real,    save, public ::
     &  ws0,
     &  ws1


c    Stokes indexes
      integer,    save, public ::
     &  ls0,
     &  ls1

C$$$$ c    Units assigned to I/O     
C$$$$ C{{{{{{{{{{{{{{{{{{{{{{{{{{c
C$$$$ #if defined(STOKES_3D)
C$$$$       integer, parameter, private :: 
C$$$$ C$$$$      .    iusdk  = 980,      !  Unit used to read usdk
C$$$$ C$$$$      .    ivsdk  = 981,      !  Unit used to read vsdk
C$$$$      .    isj    = 982,      !  Unit used to read sj
C$$$$ C$$$$      .    isdk   = 983,      !  Unit used to read sdka
C$$$$ C$$$$      .    iwadir = 984,      !  Unit used to read wave_dir
C$$$$ C$$$$      .    ihsig  = 985,      !  Unit used to read h_sig
C$$$$      .    iepsbk = 986,      !  Unit used to read eps_brk
C$$$$      .    itauwx = 987,      !  Unit used to read x wind stress lost by wave dissipation
C$$$$      .    itauwy = 988,      !  Unit used to read y wind stress lost by wave dissipation
C$$$$      .    iustar = 989,      !  Unit used to read the velocity scale of WW3
C$$$$      .    iphioc = 990,      !  Unit used to read the turbulent kinetic energy from WW3
C$$$$ C$$$$      .    iz0top = 991       !  Unit used to read the surface roughness or Hsw
C$$$$ #endif
C$$$$ #
C$$$$ #if defined(STOKES_2D)
C$$$$       integer, parameter, private :: 
C$$$$      .    istrad = 980       !  Unit used to read radiation stresses
C$$$$ #endif
C$$$$ C}}}}}}}}}}}}}}}}}}}}}}}}}}


      contains 
c*********************************************************************************************
     

      subroutine stokes_set(dtime)
      implicit none
c
      real*8 dtime
c
c --- Stokes Drift velocity setup

      if     (mnproc.eq.1) then
      write(6,*)'=================================================='
      write(6,*)'  In mod_stokes.F  allocating arrays!'
      write(6,*)'nbdy,idm,jdm,kdm = ',nbdy,idm,jdm,kdm
      write(6,*)'=================================================='
      endif !1st tile

*
*     Allocate Stokes Drift Velocity arrays used in rest of HYCOM
*
      allocate( usd(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm),
     &          vsd(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm))
      call mem_stat_add( 2*(idm+2*nbdy)*(jdm+2*nbdy)*kdm )
      usd(:,:,:) = 0.0
      vsd(:,:,:) = 0.0

      allocate( usdp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm),
     &          vsdp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm))
      call mem_stat_add( 2*(idm+2*nbdy)*(jdm+2*nbdy)*kdm )
      usdp(:,:,:) = 0.0
      vsdp(:,:,:) = 0.0

      allocate( usdbavg(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &          vsdbavg(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      call mem_stat_add( 2*(idm+2*nbdy)*(jdm+2*nbdy) )
      usdbavg(:,:) = 0.0
      vsdbavg(:,:) = 0.0

      allocate( usds(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &          vsds(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      call mem_stat_add( 2*(idm+2*nbdy)*(jdm+2*nbdy) )
      usds(:,:) = 0.0
      vsds(:,:) = 0.0

c     Arrar for transport

       allocate( tusd(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2),
     &           tvsd(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2))
       call mem_stat_add( 2*(idm+2*nbdy)*(jdm+2*nbdy)*2 )
       tusd(:,:,:) = 0.0
       tvsd(:,:,:) = 0.0
c
       allocate( tauwx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2),
     &           tauwy(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2))
       call mem_stat_add( 2*(idm+2*nbdy)*(jdm+2*nbdy)*2 )
       tauwx(:,:,:) = 0.0
       tauwy(:,:,:) = 0.0

       allocate( twomx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2),
     &           twomy(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2))
       call mem_stat_add( 2*(idm+2*nbdy)*(jdm+2*nbdy)*2 )
       twomx(:,:,:) = 0.0
       twomy(:,:,:) = 0.0
C$$ C 
C$$   allocate( tm01(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2),
C$$   call mem_stat_add( 1*(idm+2*nbdy)*(jdm+2*nbdy)*2 )
C$$    tm01(:,:,:)    = 0.0
 
c      allocate( dzdx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm),
c     &          dzdy(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm))
c      call mem_stat_add( 2*(idm+2*nbdy)*(jdm+2*nbdy)*kdm )
c      dzdx(:,:,:) = 0.0
c      dzdy(:,:,:) = 0.0

c
      ws0 = -99.0
      ws1 = -99.0

      if     (stdflg) then
        ls0 = 1
        ls1 = 2
*
*       Read in Stokes Drift profile arrays for current model start time
*
        allocate( sdzi(nsdzi) )
        call mem_stat_add(nsdzi )
       
        allocate( usdz(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2,nsdzi),
     &            vsdz(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2,nsdzi))
        call mem_stat_add( 2*(idm+2*nbdy)*(jdm+2*nbdy)*2*nsdzi )
        usdz(:,:,:,:) = 0.0
        vsdz(:,:,:,:) = 0.0
C$$ c          
C$$       allocate(usdshat(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
C$$                vsdshat(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
C$$      &         langnum(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
C$$       call mem_stat_add( 3*(idm+2*nbdy)*(jdm+2*nbdy) )
C$$       usdshat(:,:) = 0.0
C$$       vsdshat(:,:) = 0.0
C$$       langnum(:,:) = 0.0

        call stokes_forfun(dtime,1)

        if     (mnproc.eq.1) then
          write (lp,*) '...finished initializing Stokes Drift ',
     &                 'velocity Fields'
        endif !1st tile
        call xcsync(flush_lp)
      else
        nsdzi=0
        if     (mnproc.eq.1) then
          write (lp,*)'Stokes drift version of HYCOM called with ',
     &                'stdflg == 0'
          write (lp,*)'All Stokes Drift fields are set to zero!'
          write (lp,*)'No attempt is made to read WW3 Stokes Drift ',
     &                'Data'
        endif !1st tile
        call xcsync(flush_lp)
      endif     !stdflg:else
      return
      end subroutine stokes_set
c      
c----------------------------------------------------------------------
c      
      subroutine stokes_forfun(dtime,n)
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
c
      real*8    dtime
      integer   n
c
c     Stokes Drift Velocity Fields 
c
c --- input for fixed interface depths on the p-grid.
c
c --- I/O and array I/O units 927 and 928 are reserved for the entire run.
c
c --- all input fields much be defined at all grid points
c
C$$       real*8    stime(927:928)
      real*8    stime(927:935)
c
      real*8    dtime0,dtime1
      save      dtime0,dtime1
c
      real      scl1,scl2
      save      scl1,scl2
c
      character preambl(5)*79,clinex*80,cliney*80
      integer   i,ios,iunit,j,lgth,nrec,k,l,f,fmax
      logical   sdprnt, hycom_isnaninf
      real      dpthin,sum_u,sum_v,pzb,pzt
c
c --- ws0 negative on first call only.
      if     (ws0.lt.-1.0) then
c
c ---   initialize Stokes fields
c
c ---   open all stokes files.
c
        if     (mnproc.eq.1) then
        write (lp,*) ' now initializing Stokes Drift fields ...'
        endif !1st tile
        call xcsync(flush_lp)
c
        lgth = len_trim(flnmfor)
c
        call zaiopf(flnmfor(1:lgth)//'forcing.stokex.a', 'old', 927)
        if     (mnproc.eq.1) then  ! .b file from 1st tile only
        open (unit=uoff+927,file=flnmfor(1:lgth)//'forcing.stokex.b',
     &        status='old', action='read')
        read (uoff+927,'(a79)') preambl
        endif !1st tile
        call preambl_print(preambl)
c
        call zaiopf(flnmfor(1:lgth)//'forcing.stokey.a', 'old', 928)
        if     (mnproc.eq.1) then  ! .b file from 1st tile only
        open (unit=uoff+928,file=flnmfor(1:lgth)//'forcing.stokey.b',
     &        status='old', action='read')
        read (uoff+928,'(a79)') preambl
        endif !1st tile
        call preambl_print(preambl)
c
c    Stokes transport        
        call zaiopf(flnmfor(1:lgth)//'forcing.transx.a','old',929)
        if     (mnproc.eq.1) then ! .b file from 1st tile only
           open (unit=uoff+929,file=flnmfor(1:lgth)/
     &          /'forcing.transx.b',status='old', action='read')
           read (uoff+929,'(a79)') preambl
        endif                   !1st tile
        call preambl_print(preambl)

        call zaiopf(flnmfor(1:lgth)//'forcing.transy.a','old',930)
        if     (mnproc.eq.1) then ! .b file from 1st tile only
           open (unit=uoff+930,file=flnmfor(1:lgth)/
     &          /'forcing.transy.b',status='old', action='read')
           read (uoff+930,'(a79)') preambl
        endif                   !1st tile
        call preambl_print(preambl)
        
c    wind stress to be substrated (atmos-to-wave mom. flux)        
        if (stdtau) then
           call zaiopf(flnmfor(1:lgth)//'forcing.tauwx.a','old',931)
           if     (mnproc.eq.1) then ! .b file from 1st tile only
           open (unit=uoff+931,file=flnmfor(1:lgth)/
     .             /'forcing.tauwx.b',status='old', action='read')
           read (uoff+931,'(a79)') preambl
           endif !1st tile
           call preambl_print(preambl)
           call zaiopf(flnmfor(1:lgth)//'forcing.tauwy.a','old',932)
           if     (mnproc.eq.1) then  ! .b file from 1st tile only
           open (unit=uoff+932,file=flnmfor(1:lgth)/
     .           /'forcing.tauwy.b',status='old', action='read')
           read (uoff+932,'(a79)') preambl
           endif !1st tile
           call preambl_print(preambl)
        else
           tauwx(:,:,:) = 0.
           tauwy(:,:,:) = 0.
        endif
c    Wind stress to be added ( wove-to-ocean mom. flux)      
       if (stdwom) then
           call zaiopf(flnmfor(1:lgth)//'forcing.twomx.a','old',933)
           if     (mnproc.eq.1) then ! .b file from 1st tile only
           open (unit=uoff+933,file=flnmfor(1:lgth)/
     .             /'forcing.twomx.b',status='old', action='read')
           read (uoff+933,'(a79)') preambl
           endif !1st tile
           call preambl_print(preambl)
           call zaiopf(flnmfor(1:lgth)//'forcing.twomy.a','old',934)
           if     (mnproc.eq.1) then  ! .b file from 1st tile only
           open (unit=uoff+934,file=flnmfor(1:lgth)/
     .           /'forcing.twomy.b',status='old', action='read')
           read (uoff+934,'(a79)') preambl
           endif !1st tile
           call preambl_print(preambl)
        else
           twomx(:,:,:) = 0.
           twomy(:,:,:) = 0.
        endif
c  Wave mean period 
C$$$         if (altrlangmr.and. langmr.gt.0) then
C$$$            call zaiopf(flnmfor(1:lgth)//'forcing.t01.a','old',935)
C$$$            if     (mnproc.eq.1) then ! .b file from 1st tile only
C$$$            open (unit=uoff+935,file=flnmfor(1:lgth)/
C$$$      .             /'forcing.t01.b',status='old', action='read')
C$$$            read (uoff+935,'(a79)') preambl
C$$$            endif !1st tile
C$$$            call preambl_print(preambl)
C$$$         else
C$$$            tm01(:,:,:) = 0.
C$$$         endif
C$2DIM c ---   read in the Stokes Drift interface depths
C$2DIM c
C$2DIM         do k=1,nsdzi
C$2DIM           call zagetc(clinex,ios, uoff+927)
C$2DIM           if     (ios.ne.0) then
C$2DIM             if     (mnproc.eq.1) then
C$2DIM               write(lp,*)
C$2DIM               write(lp,*) 'error in stokes_forfun -' //
C$2DIM      &                    ' hit end of input forcing.stokex'
C$2DIM               write(lp,*) 'while reading the stokes depths'
C$2DIM               write(lp,*)
C$2DIM             endif !1st tile
C$2DIM             call xcstop('(stokes_forfun)')
C$2DIM                    stop '(stokes_forfun)'
C$2DIM           endif !ios
C$2DIM           call zagetc(cliney,ios, uoff+928)
C$2DIM           if     (ios.ne.0) then
C$2DIM             if     (mnproc.eq.1) then
C$2DIM               write(lp,*)
C$2DIM               write(lp,*) 'error in stokes_forfun -' //
C$2DIM      &                    ' hit end of input forcing.stokey'
C$2DIM               write(lp,*) 'while reading the stokes depths'
C$2DIM               write(lp,*)
C$2DIM             endif !1st tile
C$2DIM             call xcstop('(stokes_forfun)')
C$2DIM                    stop '(stokes_forfun)'
C$2DIM           endif !ios
C$2DIM c
C$2DIM           if     (clinex.ne.cliney) then
C$2DIM             if     (mnproc.eq.1) then
C$2DIM               write(lp,*)
C$2DIM               write(lp,*) 'error in stokes_forfun -' //
C$2DIM      &                    ' x and y depths are different'
C$2DIM               write(lp,*) 'for interface number ',k
C$2DIM               write(lp,'(a)')  trim(clinex)
C$2DIM               write(lp,'(a)')  trim(cliney)
C$2DIM               write(lp,*)
C$2DIM             endif !1st tile
C$2DIM             call xcstop('(stokes_forfun)')
C$2DIM                    stop '(stokes_forfun)'
C$2DIM           endif !clinex.ne.cliney
C$2DIM           read (clinex,*,iostat=ios) sdzi(k)
C$2DIM           if     (ios.ne.0) then
C$2DIM             if     (mnproc.eq.1) then
C$2DIM               write(lp,*)
C$2DIM               write(lp,*) 'error in stokes_forfun -' //
C$2DIM      &                    ' hit end of stokes depths in forcing.stokex'
C$2DIM               write(lp,*)
C$2DIM             endif !1st tile
C$2DIM             call xcstop('(stokes_forfun)')
C$2DIM                    stop '(stokes_forfun)'
C$2DIM           endif !ios
C$2DIM           if     (mnproc.eq.1) then
C$2DIM             write(lp,'(a)')  trim(clinex)
C$2DIM           endif !1st tile
C$2DIM           sdzi(k) = sdzi(k)*onem  !pressure units
C$2DIM           if     (sdzi(k).lt.sdzi(max(k-1,1))) then
C$2DIM             if     (mnproc.eq.1) then
C$2DIM             write(lp,*)
C$2DIM             write(lp,*) 'error in stokes_forfun -' //
C$2DIM      &                  ' Stokes depths must be non-decreasing'
C$2DIM             write(lp,*)
C$2DIM             endif !1st tile
C$2DIM             call xcstop('(stokes_forfun)')
C$2DIM                    stop '(stokes_forfun)'
C$2DIM           endif !zi.k<zi.k-1
C$2DIM         enddo  !k
C$2DIM         call xcsync(flush_lp)
C$2DIM c
C$2DIM         dpthin  = 0.001*onemm
C$2DIM         if     (sdzi(1).ne.0.0 .or. sdzi(2).lt.dpthin) then
C$2DIM           if     (mnproc.eq.1) then
C$2DIM             write(lp,*)
C$2DIM             write(lp,*) 'error in stokes_forfun -' //
C$2DIM      &                  ' 1st Stokes depth must be 0 and' //
C$2DIM      &                  ' 2nd must be positive'
C$2DIM             write(lp,*)
C$2DIM           endif !1st tile
C$2DIM           call xcstop('(stokes_forfun)')
C$2DIM                  stop '(stokes_forfun)'
C$2DIM         endif !sdzi(1:2)?
C$2DIM c
C$2DIM c ---   [uv]sdz.1 is exactly at the surface,
C$2DIM c ---   convert it to a very thin layer by stealing some of layer 2.
C$2DIM         sdzi(1) = dpthin
C$2DIM         scl2    = sdzi(2)/(sdzi(2)-dpthin)  !>1.0
C$2DIM         scl1    = scl2-1.0                  !positive, near zero
c
        dpthin  = 0.001*onemm
c ---   skip ahead to the start time.
c
        nrec   = 0
        dtime1 = huge(dtime1)
        do  ! infinite loop, with exit at end
          dtime0 = dtime1
          nrec   = nrec + 1
C$2DIM           do k=1,nsdzi
            call zagetc(clinex,ios, uoff+927)
            if     (ios.ne.0) then
              if     (mnproc.eq.1) then
              write(lp,*)
              write(lp,*) 'error in stokes_forfun -' //
     &                    ' hit end of input forcing.stokex'
              write(lp,*) 'dtime0,dtime1 = ',dtime0,dtime1
              write(lp,*) 'dtime = ',dtime
              write(lp,*)
              endif !1st tile
              call xcstop('(stokes_forfun)')
                     stop '(stokes_forfun)'
            endif !ios
            call zagetc(cliney,ios, uoff+928)
            if     (ios.ne.0) then
              if     (mnproc.eq.1) then
              write(lp,*)
              write(lp,*) 'error in stokes_forfun -' //
     &                    ' hit end of input forcing.stokey'
              write(lp,*) 'dtime0,dtime1 = ',dtime0,dtime1
              write(lp,*) 'dtime = ',dtime
              write(lp,*)
              endif !1st tile
              call xcstop('(stokes_forfun)')
                     stop '(stokes_forfun)'
            endif !ios
C$2DIM           enddo  !k
c
          i = index(clinex,'=')
          read (clinex(i+1:),*) dtime1
          if     (nrec.eq.1 .and. dtime1.lt.1462.0d0) then
c
c ---       must start after wind day 1462.0, 01/01/1905.
            if     (mnproc.eq.1) then
            write(lp,'(a)')  clinex
            write(lp,'(/ a,a / a,g15.6 /)')
     &        'error in stokes_forfun - actual forcing',
     &        ' must start after wind day 1462',
     &        'dtime1 = ',dtime1
            endif !1st tile
            call xcstop('(stokes_forfun)')
                   stop '(stokes_forfun)'
          endif !before wind day 1462.0
          if     (dtime0.le.dtime .and. dtime1.gt.dtime) then
            exit
          endif
        enddo   ! - infinite loop, with exit above
c
        if     (mnproc.eq.1) then  ! .b file from 1st tile only
          rewind(unit=uoff+927)
          read (uoff+927,'(a79)') preambl
C$2DIM           do k=1,nsdzi
C$2DIM             read (uoff+927,'(a)') clinex
C$2DIM           enddo !k
          rewind(unit=uoff+928)
          read (uoff+928,'(a79)') preambl
C$2DIM           do k=1,nsdzi
C$2DIM             read (uoff+928,'(a)') cliney
C$2DIM           enddo !k
        endif
c
        do iunit= 927,928
          do i= 1,nrec-2
C$2DIM             do k=1,nsdzi
              call skmonth(iunit)
C$2DIM             enddo
          enddo
        enddo
C$2DIM         do k=1,nsdzi
        k=1
        sdprnt = .true.
*       sdprnt = mod(nstep,nsdzi).eq.k-1
        call rdpall1(usdz(1-nbdy,1-nbdy,1,k),stime(927),927,sdprnt)
        call rdpall1(vsdz(1-nbdy,1-nbdy,1,k),stime(928),928,sdprnt)
C$2DIM         enddo !k
c ---   Read stokes transport
        do iunit= 929,930
            do i= 1,nrec-2
               call skmonth(iunit)
            enddo
        enddo    
        call rdpall1(tusd(1-nbdy,1-nbdy,1),stime(929),929,sdprnt)
        call rdpall1(tvsd(1-nbdy,1-nbdy,1),stime(930),930,sdprnt)
c --- read tau wave if stdtau is true
        if (stdtau) then
           do iunit= 931,932
              do i= 1,nrec-2
                 call skmonth(iunit)
              enddo
           enddo
           call rdpall1(tauwx(1-nbdy,1-nbdy,1),stime(931),931,sdprnt)
           call rdpall1(tauwy(1-nbdy,1-nbdy,1),stime(932),932,sdprnt)
        endif
c --- read wom wave if stdwom is true
        if (stdwom) then
           do iunit= 933,934
              do i= 1,nrec-2
                 call skmonth(iunit)
              enddo
           enddo
           call rdpall1(twomx(1-nbdy,1-nbdy,1),stime(933),933,sdprnt)
           call rdpall1(twomy(1-nbdy,1-nbdy,1),stime(934),934,sdprnt)
        endif
c --- read wave mean period if altralngmr is true
C$$$         if (altrlangmr.and. langmr.gt.0) then
C$$$            do i= 1,nrec-2
C$$$               call skmonth(935)
C$$$            enddo
C$$$            call rdpall1(tm01(1-nbdy,1-nbdy,1),stime(935),935,sdprnt)
C$$$         endif

        dtime0 = dtime1
        dtime1 = stime(927)
c        
C$2DIM         do k=1,nsdzi
        k=1 !-> Stokes surface(k=1) velocitiesc        
          sdprnt = .true.
*         sdprnt = mod(nstep,nsdzi).eq.k-1
          call rdpall1(usdz(1-nbdy,1-nbdy,1,k),stime(927),927,sdprnt)
          call rdpall1(vsdz(1-nbdy,1-nbdy,1,k),stime(928),928,sdprnt)
          call xctilr( usdz(1-nbdy,1-nbdy,1,k),1,2, nbdy,nbdy, halo_pv)
          call xctilr( vsdz(1-nbdy,1-nbdy,1,k),1,2, nbdy,nbdy, halo_pv)
C$2DIM         enddo !k
C$2DIM c ---   [uv]sdz.1 is exactly at the surface,
C$2DIM c ---   convert it to a very thin layer (modifies [uv]sdz.2).
C$2DIM         usdz(:,:,:,2) = scl2*usdz(:,:,:,2) - scl1*usdz(:,:,:,1)
C$2DIM         vsdz(:,:,:,2) = scl2*vsdz(:,:,:,2) - scl1*vsdz(:,:,:,1)
C$2DIM c
        ! Stokes transport
        call rdpall1(tusd(1-nbdy,1-nbdy,1),stime(929),929,.true.)
        call rdpall1(tvsd(1-nbdy,1-nbdy,1),stime(930),930,.true.)
        call xctilr(tusd(1-nbdy,1-nbdy,1),1,2, nbdy,nbdy, halo_uv)
        call xctilr(tvsd(1-nbdy,1-nbdy,1),1,2, nbdy,nbdy, halo_vv)
        ! tauw tau_wave
        if (stdtau) then
          call rdpall1(tauwx(1-nbdy,1-nbdy,1),stime(931),931,.true.)
          call rdpall1(tauwy(1-nbdy,1-nbdy,1),stime(932),932,.true.)
          call xctilr(tauwx(1-nbdy,1-nbdy,1),1,2, nbdy,nbdy, halo_uv)
          call xctilr(tauwy(1-nbdy,1-nbdy,1),1,2, nbdy,nbdy, halo_vv)
        endif
       ! twom wave-to-ocean
        if (stdwom) then
          call rdpall1(twomx(1-nbdy,1-nbdy,1),stime(933),933,.true.)
          call rdpall1(twomy(1-nbdy,1-nbdy,1),stime(934),934,.true.)
          call xctilr(twomx(1-nbdy,1-nbdy,1),1,2, nbdy,nbdy, halo_uv)
          call xctilr(twomy(1-nbdy,1-nbdy,1),1,2, nbdy,nbdy, halo_vv)
        endif
c
C$$$         if (altrlangmr.and. langmr.gt.0) then
C$$$           call rdpall1(tm01(1-nbdy,1-nbdy,1),stime(935),935,.true.)
C$$$           call xctilr(tm01(1-nbdy,1-nbdy,1),1,2, nbdy,nbdy, halo_uv)
C$$$         endif
c          
        dtime0 = dtime1
        dtime1 = stime(927)
      
        if     (mnproc.eq.1) then
        write (lp,*) 
        write (lp,*) ' dtime,dtime0,dtime1 = ',dtime,dtime0,dtime1
        write (lp,*) 
        write (lp,*) ' ...finished initializing Stokes Drift fields'
        endif !1st tile
        call xcsync(flush_lp)
      endif  ! --->  initialization
c
      if     (dtime.gt.dtime1) then
c
c ---   get the next set of fields.
*           if     (mnproc.eq.1) then
*           write(lp,*) 'enter rdpall - ',time,dtime0,dtime1
*           endif !1st tile
*           call xcsync(flush_lp)
C$2DIM         do k=1,nsdzi
         k=1
c         sdprnt = mod(nstep,nsdzi).eq.k-1
          sdprnt = .true.
          call rdpall1(usdz(1-nbdy,1-nbdy,1,k),stime(927),927,sdprnt)
          call rdpall1(vsdz(1-nbdy,1-nbdy,1,k),stime(928),928,sdprnt)
          call xctilr( usdz(1-nbdy,1-nbdy,2,k),1,1, nbdy,nbdy, halo_pv)
          call xctilr( vsdz(1-nbdy,1-nbdy,2,k),1,1, nbdy,nbdy, halo_pv)
C$2DIM         enddo !k
C$2DIM c ---   [uv]sdz.1 is exactly at the surface,
C$2DIM c ---   convert it to a very thin layer (modifies [uv]sdz.2).
C$2DIM         usdz(:,:,2,2) = scl2*usdz(:,:,2,2) - scl1*usdz(:,:,2,1)
C$2DIM         vsdz(:,:,2,2) = scl2*vsdz(:,:,2,2) - scl1*vsdz(:,:,2,1)
c
        ! Stokes transport
        call rdpall1(tusd(1-nbdy,1-nbdy,1),stime(929),929,.true.)
        call rdpall1(tvsd(1-nbdy,1-nbdy,1),stime(930),930,.true.)
        call xctilr(tusd(1-nbdy,1-nbdy,2),1,1, nbdy,nbdy, halo_uv)
        call xctilr(tvsd(1-nbdy,1-nbdy,2),1,1, nbdy,nbdy, halo_vv)
        
        if (stdtau) then
           call rdpall1(tauwx(1-nbdy,1-nbdy,1),stime(931),931,.true.)
           call rdpall1(tauwy(1-nbdy,1-nbdy,1),stime(932),932,.true.)
           call xctilr(tauwx(1-nbdy,1-nbdy,2),1,1, nbdy,nbdy, halo_uv)
           call xctilr(tauwy(1-nbdy,1-nbdy,2),1,1, nbdy,nbdy, halo_vv)
        endif
c
        if (stdwom) then
           call rdpall1(twomx(1-nbdy,1-nbdy,1),stime(933),933,.true.)
           call rdpall1(twomy(1-nbdy,1-nbdy,1),stime(934),934,.true.)
           call xctilr(twomx(1-nbdy,1-nbdy,2),1,1, nbdy,nbdy, halo_uv)
           call xctilr(twomy(1-nbdy,1-nbdy,2),1,1, nbdy,nbdy, halo_vv)
        endif
c        
C$$$         if (altrlangmr.and. langmr.gt.0) then
C$$$           call rdpall1(tm01(1-nbdy,1-nbdy,1),stime(935),935,.true.)
C$$$           call xctilr(tm01(1-nbdy,1-nbdy,2),1,1, nbdy,nbdy, halo_uv)
C$$$         endif
c          
        dtime0 = dtime1
        dtime1 = stime(927)
*           if     (mnproc.eq.1) then
*           write(lp,*) ' exit rdpall1 - ',time,dtime0,dtime1
*           endif !1st tile
*           call xcsync(flush_lp)
      endif
c
c --- linear interpolation in time.
      ws0 = (dtime1-dtime)/(dtime1-dtime0)
      ws1 = 1.0 - ws0
*           if     (mnproc.eq.1) then
*           write(lp,*) 'stokes_forfun - dtime,ws0,ws1 = ',
*    &                  dtime,ws0,ws1
*           endif !1st tile
*           call xcsync(flush_lp)
c--------------------------------------------------------------------
c  Calculate the p-grid Stokes Drift Arrays
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c   1. First the Surface Velocity arrays for the Mixed Layer Theories
c
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
          usds(i,j) = ws0*usdz(i,j,1,1) +
     &                ws1*usdz(i,j,2,1)
          vsds(i,j) = ws0*vsdz(i,j,1,1) +
     &                ws1*vsdz(i,j,2,1)
        enddo !i
      enddo !j
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c   2. The -k- layer Average Velocity arrays at the -p- points
c
!$OMP   PARALLEL DO PRIVATE(j)
!$OMP&               SHARED(n)
!$OMP&           SCHEDULE(STATIC,jblk)
      do j= 1,jj
        !call stokes_vertical_j(n,j) 
        call stokes_vertical_j_phil(n,j)
      enddo !j
!$OMP END PARALLEL DO
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  All Stokes Drift Velocity Arrays are now on 'private' -p-  grids 
c  Now calculate the U velocity components on -U- grids
c                and V velocity Components on -V- grids
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      call xctilr(usdp(1-nbdy,1-nbdy,1),1,kk, 1,1, halo_pv)
      call xctilr(vsdp(1-nbdy,1-nbdy,1),1,kk, 1,1, halo_pv)
c
c     write(lp,*) 'maxusdp=',maxval(usdp(:,:)),
c     &     'maxvsp=',maxval(vsdp(:,:))
c
C    !    write (lp,*) 'max usdp(:,:,1)=',
C    ! &            maxval(usdz(:,:,1)), mnproc
C    !        write (lp,*) 'min usdp(:,:,1)=',
C    ! &            minval(usdp(:,:,1)), mnproc
C    !        
!$OMP PARALLEL DO PRIVATE(j,i,k,sum_u,sum_v)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1,jj
        do i= 1,ii
          if (SEA_U) then
            sum_u = 0.0
            do k=1,kk
              usd(i,j,k)=0.5*(usdp(i,j,k)+usdp(i-1,j,k))
              sum_u = sum_u + usd(i,j,k)*dpu(i,j,k,n)
            enddo !k
            usdbavg(i,j) = sum_u/depthu(i,j)
          endif !iu
        enddo !i
        do i= 1,ii
          if (SEA_V) then
            sum_v = 0.0
            do k=1,kk
              vsd(i,j,k)=0.5*(vsdp(i,j,k)+vsdp(i,j-1,k))
              sum_v = sum_v + vsd(i,j,k)*dpv(i,j,k,n)
            enddo !k
            vsdbavg(i,j) = sum_v/depthv(i,j)
          endif !iv
        enddo !i
      enddo !j

!$OMP END PARALLEL DO

c!$OMP PARALLEL DO PRIVATE(j,i,k)
c!$OMP&         SCHEDULE(STATIC,jblk)
c
c                              Calculation of dzdx,dzdy needed in momtum
c                              -----------------------------------------
c
c      do j=1-margin, jj+margin
c        do i= 1,ii
c          if (SEA_U) then
c            ztop(i,j) = -(depths(i,j)-depths(i-1,j))
c            zbot(i,j) = -(depths(i,j)-depths(i-1,j))
c          endif !iu
c        enddo !i
c      enddo !j
c
c         do k = kk, 1, -1
c            do j = 1-margin, jj+margin
c             do i= 1,ii
c               if (SEA_U) then
c                  ztop(i,j) = ztop(i,j) +
c     .                   (dp(i,j,k,nm)*util2(i,j)/grho(i,j,k)-
c     .                    dp(i-1,j,k,nm)*util2(i-1,j)/grho(i-1,j,k))
c                  dzdx(i,j,k) = 0.5*(ztop(i,j)+zbot(i,j))
c                  zbot(i,j) = ztop(i,j)
c              endif !iu
c           enddo  !i
c         enddo  !j
c       enddo  !k
c
c        do i= 1,ii
c          if (SEA_V) then
c
c         do j = 1-margin, jj+margin 
c           do i= 1,ii
c               if (SEA_V) then
c                  ztop(i,j) = -(depths(i,j)-depths(i,j-1))
c                  zbot(i,j) = -(depths(i,j)-depths(i,j-1))      
c              endif !iv
c           enddo  !i
c         enddo  !j
c
c         do k = kk, 1, -1
c            do j = 1-margin, jj+margin
c              do i= 1,ii
c               if (SEA_V) then
c               ztop(i,j) = ztop(i,j) +
c     .                   (dp(i,j,k,nm)*util2(i,j)/grho(i,j,k)-
c     .                    dp(i,j-1,k,nm)*util2(i,j-1)/grho(i,j-1,k))
c               dzdy(i,j,k) = 0.5*(ztop(i,j)+zbot(i,j))
c               zbot(i,j) = ztop(i,j)        
c              endif !iv
c           enddo  !i
c         enddo  !j
c       enddo  !k
c!$OMP END PARALLEL DO
c
       if     (debug_stokes) then
 103    format (i9,2i5,a)
 104    format (30x,i3,2f8.4,f9.3,f9.2)
        if (itest.gt.0 .and. jtest.gt.0) then
          write (lp,103) nstep,itest+i0,jtest+j0,
     .    '  stokes_forfun:  usdz    vsdz    thkns     dpth'
          pzb = 0.0
          do k= 1,nsdzi
            pzt = pzb
            pzb = min(sdzi(k),p(itest,jtest,kk+1))*qonem
            write (lp,104)
     .      k,
     .      ws0*usdz(itest,jtest,1,k)+ws1*usdz(itest,jtest,2,k),
     .      ws0*vsdz(itest,jtest,1,k)+ws1*vsdz(itest,jtest,2,k),
     .      pzb-pzt,pzb
            if     (pzt.eq.p(itest,jtest,kk+1)*qonem) then
              exit
            endif
          enddo !k
c
          write (lp,103) nstep,itest+i0,jtest+j0,
     .    '  stokes_forfun:  usdp    vsdp     thkns     dpth'
          pzb = 0.0
          do k= 1,kk
            pzt = pzb
            pzb = min(pzt+dp(itest,jtest,k,n)*qonem,
     &                    p(itest,jtest,kk+1)*qonem)
            write (lp,104)
     .      k,
     .      usdp(itest,jtest,k),vsdp(itest,jtest,k),
     .      pzb-pzt,pzb
            if     (pzt.eq.p(itest,jtest,kk+1)*qonem) then
              exit
            endif
          enddo !k
c
          write (lp,103) nstep,itest+i0,jtest+j0,
     .    '  stokes_forfun:   usd     vsd     thkns     dpth'
          write (lp,104)
     .    0,usdbavg(itest,jtest),vsdbavg(itest,jtest),
     .    0.0,p(itest,jtest,kk+1)*qonem 
          pzb = 0.0
          do k= 1,kk
            pzt = pzb
            pzb = min(pzt+dp(itest,jtest,k,n)*qonem,
     &                    p(itest,jtest,kk+1)*qonem)
            write (lp,104)
     .      k,
     .      usd(itest,jtest,k),vsd(itest,jtest,k),
     .      pzb-pzt,pzb
            if     (pzt.eq.p(itest,jtest,kk+1)*qonem) then
              exit
            endif
          enddo !k
        endif !test
      endif !debug  
c
*      if((i.eq.1).and.(j.eq.1))then
*         write(6,*)'U and V grid Stokes Drift velocities at point (1,1)'
*         write(6,*)'u/v sdbavg = ',usdbavg(1,1),vsdbavg(1,1)
*         write(6,*)'        k     usd(1,1,k)       vsd(1,1,k)'
*         do k=1,kk
*            write(6,*)k,usd(1,1,k),vsd(1,1,k)
*         enddo
*      endif

c
c   Now ensure  all Stokes Drift Velocity Fields are properly defined on Halos
c

      call xctilr(usd(    1-nbdy,1-nbdy,1),1,kk, 6,6, halo_uv)
      call xctilr(vsd(    1-nbdy,1-nbdy,1),1,kk, 6,6, halo_vv)
      call xctilr(usdbavg(1-nbdy,1-nbdy)  ,1, 1, 6,6, halo_uv)
      call xctilr(vsdbavg(1-nbdy,1-nbdy)  ,1, 1, 6,6, halo_vv)
c
c      call xctilr(dzdx(1-nbdy,1-nbdy,1),1, kk, nbdy,nbdy, halo_uv)
c      call xctilr(dzdy(1-nbdy,1-nbdy,1),1, kk, nbdy,nbdy, halo_vv)

      return
      end subroutine stokes_forfun

      subroutine stokes_vertical_j(n,j)
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
c
      integer   n,j
c
c --- --------------------------------------------
c --- interpolate Stokes in vertical, single j-row
c --- --------------------------------------------
c
      integer i,k
      real    dpthin
      logical lcm(nsdzi)      !use PCM for some layers?
      real    s1d(nsdzi,2),   !input Stokes fields
     &        f1d(kdm,  2),   !model Stokes fields
     &        c1d(nsdzi,2,3), !interpolation coefficients
     &        dpi( nsdzi),    !input layer thicknesses, >= dpthin
     &        dprs(nsdzi),    !input layer thicknesses
     &        pres(nsdzi+1),  !input layer interfaces
     &        prsf(kdm+1)     !model layer interfaces
c
      dpthin = 0.001*onemm
      lcm(:) = .false.
c
      do i= 1,ii
        if (SEA_P) then
c ---     1-d existing arrays
          pres(1)=0.0
          do k=1,nsdzi
            s1d(k,1) = ws0*usdz(i,j,1,k) +
     &                 ws1*usdz(i,j,2,k)
            s1d(k,2) = ws0*vsdz(i,j,1,k) +
     &                 ws1*vsdz(i,j,2,k)
            pres(k+1)=min(sdzi(k),p(i,j,kk+1))
            dprs(k)  =pres(k+1)-pres(k)
            dpi( k)  =max(dprs(k),dpthin)
          enddo !k 1:nsdzi
          prsf(1)=0.0
          do k=1,kk
            prsf(k+1) = p(i,j,k+1)
          enddo !k 1:kk
c ---     remap
          call hybgen_weno_coefs(s1d,     dpi, lcm,c1d,
     &                                    nsdzi,   2,dpthin)
          call hybgen_weno_remap(s1d,pres,dprs,    c1d,
     &                           f1d,prsf,nsdzi,kk,2,dpthin)
c ---     1-d new arrays
          do k=1,kk
            usdp(i,j,k) = f1d(k,1)
            vsdp(i,j,k) = f1d(k,2)
          enddo !k 1:kk
        endif !ip
      enddo !i
c
      return
      end subroutine stokes_vertical_j
c     ><<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<
C   !  Construct stokes profile using Phillips approx.,
C      Breivik et al., 2016. 
C      Need on the Stokes transport and surface velocity
c     ><<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<< 
      subroutine stokes_vertical_j_phil(n,j)
      use mod_xc                ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
c
      integer   n,j
c
      integer i,k,l
      real    knm(2), sum0   !stokes depth  
      real    ,parameter :: betaa=1.0
      real    zbot,ztop, Vzbot,Vztop    !stokes depth  
      logical hycom_isnaninf
      real    ts1d(2)     !input Stokes trans. 
      real    f1dphil(kdm,2),       !Phillips actua. prof. 
     &        f1dphilavg(kdm,2),    !Phillips averaged prof.
     &        prsf(kdm+1)           !model layer interfaces
c
      real  v0sp,tvsp,thetatv0,kbar,v0spzi,
     $     Vztopsp,Vzbotsp,v0spzilayer 
      
c --- ---------------------------------------------------------
c --- Construct Stokes in vertical using Phillips, single j-row
c --- ---------------------------------------------------------
c
         do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))   
c     --- hycom layers
C$2DIM              do k=1,nsdzi   ! k looooooooop
C$2DIM                 s1d(k,1) = ws0*usdz(i,j,1,k) +
C$2DIM      &                     ws1*usdz(i,j,2,k)
C$2DIM                 s1d(k,2) = ws0*vsdz(i,j,1,k) +
C$2DIM      &                     ws1*vsdz(i,j,2,k)
C$2DIM               enddo             !k 1:nsdzi  
C    --- input transport at i,j (time interp.)
              ts1d(1) = ws0*tusd(i,j,1) +
     &                  ws1*tusd(i,j,2)
              ts1d(2) = ws0*tvsd(i,j,1) +
     &                  ws1*tvsd(i,j,2)
C    ---  comp speed and angle of the transport
C    ---  comp stokes surface speed              
              v0sp=sqrt(usds(i,j)**2+vsds(i,j)**2) 
              tvsp=sqrt(ts1d(1)**2+ts1d(2)**2)
              thetatv0=atan2(ts1d(2),ts1d(1))
C    ---  comp Stokes depth  
              kbar=v0sp/(2.0*max(tvsp,epsil))
     $            *(1.0-2.0*betaa/3.0)
c    --  Now comp layer averaged profile
         if     (debug_stokes) then    
          if((i+i0).eq. 430 .and. (j+j0).eq. 500)then
  103      format (i9,2i5,a)
            write (lp,103) nstep,itest+i0,jtest+j0,
     .      ' stokes_forfun: Phlavgsp  Phlsp  usdz      dpth'
          endif  
         endif
         sum0=0.0
         prsf(1)=0.0
         do k=1,kk
c            prsf(k+1) = p(i,j,k+1)
            prsf(k+1)=min(prsf(k)+dp(i,j,k,n),p(i,j,kk+1))
            !   layer k:  zbot,ztop
            ztop=prsf(k)/onem    ! top    
            zbot=prsf(k+1)/onem  ! bot
            ! 
c   --   actual Stokes profile (not really needed, just for debugging)                   
            v0spzi=phillips_profile(v0sp,
     $          kbar,-1*ztop,betaa)                   
            f1dphil(k,1)=v0spzi*cos(thetatv0)
            f1dphil(k,2)=v0spzi*sin(thetatv0)
c   --      layer_averaged Stokes profile               
            Vztopsp=phillips_transportz(v0sp,
     $                              kbar,-1*ztop,betaa)
            Vzbotsp=phillips_transportz(v0sp,
     $                              kbar,-1*zbot,betaa)
            v0spzilayer=(Vzbotsp-Vztopsp)/max(zbot-ztop,epsil)
            f1dphilavg(k,1)=v0spzilayer*cos(thetatv0)
            f1dphilavg(k,2)=v0spzilayer*sin(thetatv0)
C     ------ Finish   computing layeravg stokes
C            now update assign stokes at p points
            usdp(i,j,k) = f1dphilavg(k,1)
            vsdp(i,j,k) = f1dphilavg(k,2)
c   >---------------------------------------------------------
C$$        if     (debug_stokes) then
C$$  103    format (i9,2i5,a)
C$$  104    format (30x,i3,2f8.4,f9.3,f9.2)
C$$         if (itest.gt.0 .and. jtest.gt.0) then
C$$           write (lp,103) nstep,itest+i0,jtest+j0,
C$$      .    '  stokes_forfun:  usdz    vsdz    thkns     dpth'
C$$           pzb = 0.0
C$$           do k= 1,nsdzi
C$$             pzt = pzb
C$$             pzb = min(sdzi(k),p(itest,jtest,kk+1))*qonem
C$$             write (lp,104)
C$$      .      k,
C$$      .      ws0*usdz(itest,jtest,1,k)+ws1*usdz(itest,jtest,2,k),
C$$      .      ws0*vsdz(itest,jtest,1,k)+ws1*vsdz(itest,jtest,2,k),
C$$      .      pzb-pzt,pzb
C$$             if     (pzt.eq.p(itest,jtest,kk+1)*qonem) then
C$$               exit
C$$             endif
C$$ C$$           enddo !k
C$$       endif
        if     (debug_stokes) then        
        if((i+i0).eq. 430 .and. (j+j0).eq. 500)then
c  103    format (i9,2i5,a)
  104    format (30x,i3,3f8.4,f9.2)
  105    format (30x,i3,f8.4,f9.2)         
c           write (lp,103) nstep,itest+i0,jtest+j0,
c     .    '  stokes_forfun:  Philavgsp   Philsp  usdz    dpth'
             sum0 = sum0 
     $              + sqrt(f1dphilavg(k,1)**2
     $                     +f1dphilavg(k,2)**2)*dp(i,j,k,n)/onem
           if (k.le.nsdzi)then 
             write (lp,104)
     .       k,
     $       sqrt(f1dphilavg(k,1)**2+ f1dphilavg(k,2)**2),
     $       sqrt(f1dphil(k,1)**2+ f1dphil(k,2)**2),  
     $       v0sp,
     $       ztop  
           else
              write (lp,105)
     .        k,
     $        sqrt(f1dphilavg(k,1)**2+ f1dphilavg(k,2)**2),
     $        ztop
              if     (ztop.eq.p(itest,jtest,kk+1)*qonem) then
                  exit
              endif
             endif
           endif
        endif !debugstokes
c     >-------------------------------------------------------            
          enddo  !k 1:kk
          if  (debug_stokes) then                
            if((i+i0).eq. 430 .and. (j+j0).eq. 500)then
              write(lp,'(a,2E16.8)')
     $                     '-->>> T_phil_avg=,T_in  ', sum0
     $           ,sqrt(ts1d(1)**2 +ts1d(2)**2 )         
            endif
           endif
          enddo  !i
       enddo   !l

C$$$c        !##############################################Phillips
C$$$cy c     !##############################################Phillips     
      return
      end subroutine stokes_vertical_j_phil
c      
      real function phillips_transportz(v0sur,knm,zlvl,bta)
C  ---   comp stokes transport between surface and zlevel 
        IMPLICIT NONE
        REAL :: v0sur, knm, zlvl, bta
        real,    parameter :: epsil=1.e-20
C     REAL             :: s
        phillips_transportz=v0sur/(2*max(knm,epsil))
     $       *(1 - exp(2*knm*zlvl) 
     $       - (2*bta/3)*(1+sqrt(pi)*(-2*knm*zlvl)**(1.5)
     $       *get_erfc(sqrt(-2*knm*zlvl)) -
     $       (1 -2*knm*zlvl)*exp(2*knm*zlvl)))
        
      END FUNCTION phillips_transportz
      
      real function phillips_profile(v0surp,knmp,zlvlp,btap)
C  ---   comp stokes profile based on Philips       
      IMPLICIT NONE
      REAL :: v0surp, knmp, zlvlp, btap      
      phillips_profile=v0surp*(exp(2*knmp*zlvlp)
     $     - btap*sqrt(-2*pi*knmp*zlvlp)
     $     *get_erfc(sqrt(-2*knmp*zlvlp)))
!     print *,'z, phill_prof =',zlvlp,phillips_profile     
      END FUNCTION phillips_profile
c
c      
      end module mod_stokes
c> May  2014 - use land/sea masks (e.g. ip) to skip land

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
#if defined(NERSC_stokes)
      use m_get_erfc
#endif
!
      implicit none
!
! --- HYCOM Stokes Drift from external data files
!
      logical, parameter, private :: debug_stokes=.false.  !usually .false.
#if defined(NERSC_stokes)
      logical, parameter, private :: debug_stokesPhl=.false.  !usually .false.
#endif
!
!     Logical flags for Stokes Drift Effects
!     The modification of the turbulent viscosity near the surface 
!     has been implemented  following:
!          - Craig & Banner + Ardhuin-Filipot-Perenne
!          - Uchiyama-McWilliams-Shchepetkin
!
      logical, save, public :: &
       stdflg,    & ! Stokes Drift Velocities:           TRUE/FALSE
       stdsur,    & ! Stokes Drift Surface Stresses:     TRUE/FALSE (dissipation due to breaking and rolling)
       stdbot,    & ! Stokes Drift Bottom Stress:        TRUE/FALSE (dissipation due to bottom friction)
#if defined(NERSC_stokes)
       stdtau,    & ! Removing from the total wind stress the part transfered to the wave field subjected to dissipation: TRUE/FALSE
       stdwom,    & ! Adding the wave-to-ocean momentum flux due to wave dissipation/breaking: TRUE/FALSE 
#endif
       stdarc       ! Stokes Drift Velocities in Archive TRUE/FALSE
!
!    Constant values
!     Number of fixed interface depths for Stokes input
!
      integer,   save, public   :: &
        nsdzi,    & ! Number of fixed interface depths for Stokes input 
        langmr      ! Langmuir turb enhancement (KPP) 0: None 1:McWilliams-Sulliva 2:Smyth 3:McWilliams-Harcourt 4:Takaya

!     Arrays holding Surface Stokes Velocities

      real,    allocatable, dimension(:,:), &
               save, public  :: &
       usds,     & ! Surface Stokes Drift U Velocity, p-grid
       vsds     ! Surface Stokes Drift V Velocity, p-grid

!     -U- & -V- grid Arrays holding Vertically averaged Stokes Velocities

      real,    allocatable, dimension(:,:), &
               save, public  :: &
       usdbavg,     & ! Vertical average Stokes Drift U Velocity, u-grid
       vsdbavg     ! Vertical average Stokes Drift V Velocity, v-grid

!     Arrays holding Stokes Drift Velocities on U and V grids

      real,    allocatable, dimension(:,:,:), &
               save, public  :: &
       usd,     & ! Stokes Drift U Velocity, u-grid
       vsd     ! Stokes Drift V Velocity, v-grid

!     Arrays holding Stokes Drift Velocities on -p- grids

      real,    allocatable, dimension(:,:,:), &
               save, public  :: &
       usdp,     & ! Stokes Drift U Velocity, p-grid
       vsdp     ! Stokes Drift V Velocity, p-grid


!     Arrays to hold Input Stokes Drift Velocity

      real,    allocatable, dimension(:), &
               save, private :: &
       sdzi       ! Input Stokes Drift layer depths

      real,    allocatable, dimension(:,:,:,:), &
               save, private :: &
       usdz,     & ! Stokes Drift U Velocity for fixed layers, p-grid
       vsdz     ! Stokes Drift V Velocity for fixed layers, p-grid

!     Array holding the wave-induced mean pressure

      real,    allocatable, dimension(:,:,:), &
               save, public  :: &
       sj     ! wave-induced mean pressure on pressure grid
!     Arrays holding the horizontal derivatives of the vertical position of the middle of the layer
!
!      real,    allocatable, dimension(:,:,:),
!     &         save, public  ::
!     & dzdx,    ! derivative of z with respect to x
!     & dzdy     ! derivative of z with respect to y

!     Arrays used to calculate the dissipations

      real,    allocatable, dimension(:,:,:), &
               save, public  :: &
       sdka,         & ! frequency as a function ot time and space (=sdk)
       wave_dir,     & ! main wave propagation direction
       h_sig,        & ! significant wave height
       eps_brk      ! wave breaking dissipation rate

!     Dissipations (bottom and breakin)

      real,    allocatable, dimension(:,:,:), &
               save, public  :: &
       wave_bdx(:,:,:),         & ! Bottom dissipation  
       wave_bdy(:,:,:),         & ! Bottom dissipation  
       wave_brkx(:,:,:),        & ! Dissipation  by breaking
       wave_brky(:,:,:)           ! Dissipation  by breaking
#if defined(NERSC_stokes)
!     Arrays used to calculate transport

      real,    allocatable, dimension(:,:,:),
           save, private  :: &
        tusd,       & ! Stokes transport u direction, p-grid
        tvsd       ! Stokes transport v direction, p-grid
#endif

!     Arrays used for the parameterization of the vertical mixing

      real,    allocatable, dimension(:,:,:), &
               save, public  :: &
       phi_ocw,   & ! Turbulent kinetic energy flux (notation from WW3)
       z0topw,    & ! wind sea significant height(if stz0tp=1) or surface roughness
       ustarw3w,  & ! Velocity scale from WW3
       difx_wave ! Vertical diffusion enhancement due to wave breaking (Uchiyama et al.)
      real,    allocatable, dimension(:,:), &
               save, public  :: &
       phi_oc,    & ! Turbulent kinetic energy flux (notation from WW3)
       alpha,     & ! Wave energy factor from WW3
       z0top,     & ! Surface roughness from WW3
       ustarw3   ! Velocity scale from WW3
!

!    Weights for holding values during a barotropic time step.

      real,    save, public :: &
        ws0, &
        ws1


!    Stokes indexes
      integer,    save, public :: &
        ls0, &
        ls1

      contains 
!*********************************************************************************************
     

      subroutine stokes_set(dtime)
      implicit none
!
      real*8 dtime
!
! --- Stokes Drift velocity setup

      if     (mnproc.eq.1) then
      write(6,*)'=================================================='
      write(6,*)'  In mod_stokes.F  allocating arrays!'
      write(6,*)'nbdy,idm,jdm,kdm = ',nbdy,idm,jdm,kdm
      write(6,*)'=================================================='
      endif !1st tile

!
! --- Allocate Stokes Drift Velocity arrays used in rest of HYCOM
!
      allocate( usd(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm), &
                vsd(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm))
      call mem_stat_add( 2*(idm+2*nbdy)*(jdm+2*nbdy)*kdm )
      usd(:,:,:) = 0.0
      vsd(:,:,:) = 0.0

      allocate( usdp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm), &
                vsdp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm))
      call mem_stat_add( 2*(idm+2*nbdy)*(jdm+2*nbdy)*kdm )
      usdp(:,:,:) = 0.0
      vsdp(:,:,:) = 0.0

      allocate( usdbavg(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                vsdbavg(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      call mem_stat_add( 2*(idm+2*nbdy)*(jdm+2*nbdy) )
      usdbavg(:,:) = 0.0
      vsdbavg(:,:) = 0.0

      allocate( usds(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                vsds(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      call mem_stat_add( 2*(idm+2*nbdy)*(jdm+2*nbdy) )
      usds(:,:) = 0.0
      vsds(:,:) = 0.0

#if defined(NERSC_stokes)
!     Arrar for transport
       allocate( tusd(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
                 tvsd(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2))
       call mem_stat_add( 2*(idm+2*nbdy)*(jdm+2*nbdy)*2 )
       tusd(:,:,:) = 0.0
       tvsd(:,:,:) = 0.0
!
#endif
!      allocate( dzdx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm),
!     &          dzdy(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm))
!      call mem_stat_add( 2*(idm+2*nbdy)*(jdm+2*nbdy)*kdm )
!      dzdx(:,:,:) = 0.0
!      dzdy(:,:,:) = 0.0

!
      ws0 = -99.0
      ws1 = -99.0

      if     (stdflg) then
        ls0 = 1
        ls1 = 2
!
! ---   Read in Stokes Drift profile arrays for current model start time
!
        allocate( sdzi(nsdzi) )

        allocate( usdz(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2,nsdzi), &
                  vsdz(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2,nsdzi))
        call mem_stat_add( 2*(idm+2*nbdy)*(jdm+2*nbdy)*2*nsdzi )
        usdz(:,:,:,:) = 0.0
        vsdz(:,:,:,:) = 0.0

        call stokes_forfun(dtime,1)

        if     (mnproc.eq.1) then
          write (lp,*) '...finished initializing Stokes Drift ', &
                       'velocity Fields'
        endif !1st tile
        call xcsync(flush_lp)
      else
        nsdzi=0
        if     (mnproc.eq.1) then
          write (lp,*)'Stokes drift version of HYCOM called with ', &
                      'stdflg == 0'
          write (lp,*)'All Stokes Drift fields are set to zero!'
          write (lp,*)'No attempt is made to read WW3 Stokes Drift ', &
                      'Data'
        endif !1st tile
        call xcsync(flush_lp)
      endif     !stdflg:else
      return
      end subroutine stokes_set

      subroutine stokes_forfun(dtime,n)
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
!
      real*8    dtime
      integer   n
!
!     Stokes Drift Velocity Fields 
!
! --- input for fixed interface depths on the p-grid.
!
! --- I/O and array I/O units 927 and 928 are reserved for the entire run.
!
! --- all input fields much be defined at all grid points
!
      real*8    stime(927:928)
!
      real*8    dtime0,dtime1
      save      dtime0,dtime1
!
      real      scl1,scl2
      save      scl1,scl2
!
      character preambl(5)*79,clinex*80,cliney*80
      integer   i,ios,iunit,j,lgth,nrec,k,f,fmax
      logical   sdprnt
      real      dpthin,sum_u,sum_v,pzb,pzt
!
! --- ws0 negative on first call only.
      if     (ws0.lt.-1.0) then
!
! ---   initialize Stokes fields
!
! ---   open all stokes files.
!
        if     (mnproc.eq.1) then
        write (lp,*) ' now initializing Stokes Drift fields ...'
        endif !1st tile
        call xcsync(flush_lp)
!
        lgth = len_trim(flnmfor)
!
        call zaiopf(flnmfor(1:lgth)//'forcing.stokex.a', 'old', 927)
        if     (mnproc.eq.1) then  ! .b file from 1st tile only
        open (unit=uoff+927,file=flnmfor(1:lgth)//'forcing.stokex.b', &
              status='old', action='read')
        read (uoff+927,'(a79)') preambl
        endif !1st tile
        call preambl_print(preambl)
!
        call zaiopf(flnmfor(1:lgth)//'forcing.stokey.a', 'old', 928)
        if     (mnproc.eq.1) then  ! .b file from 1st tile only
        open (unit=uoff+928,file=flnmfor(1:lgth)//'forcing.stokey.b', &
              status='old', action='read')
        read (uoff+928,'(a79)') preambl
        endif !1st tile
        call preambl_print(preambl)
!
! ---   read in the Stokes Drift interface depths
!
        do k=1,nsdzi
          call zagetc(clinex,ios, uoff+927)
          if     (ios.ne.0) then
            if     (mnproc.eq.1) then
              write(lp,*)
              write(lp,*) 'error in stokes_forfun -' // &
                          ' hit end of input forcing.stokex'
              write(lp,*) 'while reading the stokes depths'
              write(lp,*)
            endif !1st tile
            call xcstop('(stokes_forfun)')
                   stop '(stokes_forfun)'
          endif !ios
          call zagetc(cliney,ios, uoff+928)
          if     (ios.ne.0) then
            if     (mnproc.eq.1) then
              write(lp,*)
              write(lp,*) 'error in stokes_forfun -' // &
                          ' hit end of input forcing.stokey'
              write(lp,*) 'while reading the stokes depths'
              write(lp,*)
            endif !1st tile
            call xcstop('(stokes_forfun)')
                   stop '(stokes_forfun)'
          endif !ios
!
          if     (clinex.ne.cliney) then
            if     (mnproc.eq.1) then
              write(lp,*)
              write(lp,*) 'error in stokes_forfun -' // &
                          ' x and y depths are different'
              write(lp,*) 'for interface number ',k
              write(lp,'(a)')  trim(clinex)
              write(lp,'(a)')  trim(cliney)
              write(lp,*)
            endif !1st tile
            call xcstop('(stokes_forfun)')
                   stop '(stokes_forfun)'
          endif !clinex.ne.cliney
          read (clinex,*,iostat=ios) sdzi(k)
          if     (ios.ne.0) then
            if     (mnproc.eq.1) then
              write(lp,*)
              write(lp,*) 'error in stokes_forfun -' // &
                          ' hit end of stokes depths in forcing.stokex'
              write(lp,*)
            endif !1st tile
            call xcstop('(stokes_forfun)')
                   stop '(stokes_forfun)'
          endif !ios
          if     (mnproc.eq.1) then
            write(lp,'(a)')  trim(clinex)
          endif !1st tile
          sdzi(k) = sdzi(k)*onem  !pressure units
          if     (sdzi(k).lt.sdzi(max(k-1,1))) then
            if     (mnproc.eq.1) then
            write(lp,*)
            write(lp,*) 'error in stokes_forfun -' // &
                        ' Stokes depths must be non-decreasing'
            write(lp,*)
            endif !1st tile
            call xcstop('(stokes_forfun)')
                   stop '(stokes_forfun)'
          endif !zi.k<zi.k-1
        enddo  !k
        call xcsync(flush_lp)
!
        dpthin  = 0.001*onemm
        if     (sdzi(1).ne.0.0 .or. sdzi(2).lt.dpthin) then
          if     (mnproc.eq.1) then
            write(lp,*)
            write(lp,*) 'error in stokes_forfun -' // &
                        ' 1st Stokes depth must be 0 and' // &
                        ' 2nd must be positive'
            write(lp,*)
          endif !1st tile
          call xcstop('(stokes_forfun)')
                 stop '(stokes_forfun)'
        endif !sdzi(1:2)?
!
! ---   [uv]sdz.1 is exactly at the surface,
! ---   convert it to a very thin layer by stealing some of layer 2.
        sdzi(1) = dpthin
        scl2    = sdzi(2)/(sdzi(2)-dpthin)  !>1.0
        scl1    = scl2-1.0                  !positive, near zero
!
! ---   skip ahead to the start time.
!
        nrec   = 0
        dtime1 = huge(dtime1)
        do  ! infinite loop, with exit at end
          dtime0 = dtime1
          nrec   = nrec + 1
          do k=1,nsdzi
            call zagetc(clinex,ios, uoff+927)
            if     (ios.ne.0) then
              if     (mnproc.eq.1) then
              write(lp,*)
              write(lp,*) 'error in stokes_forfun -' // &
                          ' hit end of input forcing.stokex'
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
              write(lp,*) 'error in stokes_forfun -' // &
                          ' hit end of input forcing.stokey'
              write(lp,*) 'dtime0,dtime1 = ',dtime0,dtime1
              write(lp,*) 'dtime = ',dtime
              write(lp,*)
              endif !1st tile
              call xcstop('(stokes_forfun)')
                     stop '(stokes_forfun)'
            endif !ios
          enddo  !k
!
          i = index(clinex,'=')
          read (clinex(i+1:),*) dtime1
          if     (nrec.eq.1 .and. dtime1.lt.1462.0d0) then
!
! ---       must start after wind day 1462.0, 01/01/1905.
            if     (mnproc.eq.1) then
            write(lp,'(a)')  clinex
            write(lp,'(/ a,a / a,g15.6 /)') &
              'error in stokes_forfun - actual forcing', &
              ' must start after wind day 1462', &
              'dtime1 = ',dtime1
            endif !1st tile
            call xcstop('(stokes_forfun)')
                   stop '(stokes_forfun)'
          endif !before wind day 1462.0
          if     (dtime0.le.dtime .and. dtime1.gt.dtime) then
            exit
          endif
        enddo   ! infinite loop, with exit above
        if     (mnproc.eq.1) then  ! .b file from 1st tile only
          rewind(unit=uoff+927)
          read (uoff+927,'(a79)') preambl
          do k=1,nsdzi
            read (uoff+927,'(a)') clinex
          enddo !k
          rewind(unit=uoff+928)
          read (uoff+928,'(a79)') preambl
          do k=1,nsdzi
            read (uoff+928,'(a)') cliney
          enddo !k
        endif
        do iunit= 927,928
          do i= 1,nrec-2
            do k=1,nsdzi
              call skmonth(iunit)
            enddo
          enddo
        enddo
        do k=1,nsdzi
          sdprnt = .true.
!         sdprnt = mod(nstep,nsdzi).eq.k-1
          call rdpall1(usdz(1-nbdy,1-nbdy,1,k),stime(927),927,sdprnt)
          call rdpall1(vsdz(1-nbdy,1-nbdy,1,k),stime(928),928,sdprnt)
        enddo !k
        dtime0 = dtime1
        dtime1 = stime(927)
        do k=1,nsdzi
          sdprnt = .true.
!         sdprnt = mod(nstep,nsdzi).eq.k-1
          call rdpall1(usdz(1-nbdy,1-nbdy,1,k),stime(927),927,sdprnt)
          call rdpall1(vsdz(1-nbdy,1-nbdy,1,k),stime(928),928,sdprnt)
          call xctilr( usdz(1-nbdy,1-nbdy,1,k),1,2, nbdy,nbdy, halo_pv)
          call xctilr( vsdz(1-nbdy,1-nbdy,1,k),1,2, nbdy,nbdy, halo_pv)
        enddo !k
! ---   [uv]sdz.1 is exactly at the surface,
! ---   convert it to a very thin layer (modifies [uv]sdz.2).
        usdz(:,:,:,2) = scl2*usdz(:,:,:,2) - scl1*usdz(:,:,:,1)
        vsdz(:,:,:,2) = scl2*vsdz(:,:,:,2) - scl1*vsdz(:,:,:,1)
!
        dtime0 = dtime1
        dtime1 = stime(927)
      
        if     (mnproc.eq.1) then
        write (lp,*) 
        write (lp,*) ' dtime,dtime0,dtime1 = ',dtime,dtime0,dtime1
        write (lp,*) 
        write (lp,*) ' ...finished initializing Stokes Drift fields'
        endif !1st tile
        call xcsync(flush_lp)
      endif  ! initialization
!
      if     (dtime.gt.dtime1) then
!
! ---   get the next set of fields.
!diag           if     (mnproc.eq.1) then
!diag           write(lp,*) 'enter rdpall - ',time,dtime0,dtime1
!diag           endif !1st tile
!diag           call xcsync(flush_lp)
        do k=1,nsdzi
          sdprnt = mod(nstep,nsdzi).eq.k-1
          call rdpall1(usdz(1-nbdy,1-nbdy,1,k),stime(927),927,sdprnt)
          call rdpall1(vsdz(1-nbdy,1-nbdy,1,k),stime(928),928,sdprnt)
          call xctilr( usdz(1-nbdy,1-nbdy,2,k),1,1, nbdy,nbdy, halo_pv)
          call xctilr( vsdz(1-nbdy,1-nbdy,2,k),1,1, nbdy,nbdy, halo_pv)
        enddo !k
! ---   [uv]sdz.1 is exactly at the surface,
! ---   convert it to a very thin layer (modifies [uv]sdz.2).
        usdz(:,:,2,2) = scl2*usdz(:,:,2,2) - scl1*usdz(:,:,2,1)
        vsdz(:,:,2,2) = scl2*vsdz(:,:,2,2) - scl1*vsdz(:,:,2,1)
!
        dtime0 = dtime1
        dtime1 = stime(927)
!diag           if     (mnproc.eq.1) then
!diag           write(lp,*) ' exit rdpall1 - ',time,dtime0,dtime1
!diag           endif !1st tile
!diag           call xcsync(flush_lp)
      endif
!
! --- linear interpolation in time.
      ws0 = (dtime1-dtime)/(dtime1-dtime0)
      ws1 = 1.0 - ws0
      if     (debug_stokes) then
        if     (mnproc.eq.1) then
        write(lp,'(a,i2,f20.10,2f10.7)') &
                  'stokes_forfun - n,dtime,ws0,ws1 = ', &
                  n,dtime,ws0,ws1
        endif !1st tile
        call xcsync(flush_lp)
      endif !debug
!--------------------------------------------------------------------
!  Calculate the p-grid Stokes Drift Arrays
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!
!   1. First the Surface Velocity arrays for the Mixed Layer Theories
!
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
          usds(i,j) = ws0*usdz(i,j,1,1) + &
                      ws1*usdz(i,j,2,1)
          vsds(i,j) = ws0*vsdz(i,j,1,1) + &
                      ws1*vsdz(i,j,2,1)
        enddo !i
      enddo !j
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!
!   2. The -k- layer Average Velocity arrays at the -p- points
!
!$OMP   PARALLEL DO PRIVATE(j) &
!$OMP                SHARED(n) &
!$OMP            SCHEDULE(STATIC,jblk)
      do j= 1,jj
        call stokes_vertical_j(n,j)
      enddo !j
!$OMP END PARALLEL DO
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  All Stokes Drift Velocity Arrays are now on 'private' -p-  grids 
!  Now calculate the U velocity components on -U- grids
!                and V velocity Components on -V- grids
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      call xctilr(usdp(1-nbdy,1-nbdy,1),1,kk, 1,1, halo_pv)
      call xctilr(vsdp(1-nbdy,1-nbdy,1),1,kk, 1,1, halo_pv)
!
!$OMP PARALLEL DO PRIVATE(j,i,k,sum_u,sum_v) &
!$OMP          SCHEDULE(STATIC,jblk)
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
!
      if     (debug_stokes) then
 103    format (i9,2i5,a)
 104    format (30x,i3,2f8.4,f9.3,f9.2)
        if (itest.gt.0 .and. jtest.gt.0) then
          write (lp,103) nstep,itest+i0,jtest+j0, &
          '  stokes_forfun:  usdz    vsdz    thkns     dpth'
          pzb = 0.0
          do k= 1,nsdzi
            pzt = pzb
            pzb = min(sdzi(k),p(itest,jtest,kk+1))*qonem
            write (lp,104) &
            k, &
            ws0*usdz(itest,jtest,1,k)+ws1*usdz(itest,jtest,2,k), &
            ws0*vsdz(itest,jtest,1,k)+ws1*vsdz(itest,jtest,2,k), &
            pzb-pzt,pzb
            if     (pzt.eq.p(itest,jtest,kk+1)*qonem) then
              exit
            endif
          enddo !k
!
          write (lp,103) nstep,itest+i0,jtest+j0, &
          '  stokes_forfun:  usdp    vsdp     thkns     dpth'
          pzb = 0.0
          do k= 1,kk
            pzt = pzb
            pzb = min(pzt+dp(itest,jtest,k,n)*qonem, &
                          p(itest,jtest,kk+1)*qonem)
            write (lp,104) &
            k, &
            usdp(itest,jtest,k),vsdp(itest,jtest,k), &
            pzb-pzt,pzb
            if     (pzt.eq.p(itest,jtest,kk+1)*qonem) then
              exit
            endif
          enddo !k
!
          write (lp,103) nstep,itest+i0,jtest+j0, &
          '  stokes_forfun:   usd     vsd     thkns     dpth'
          write (lp,104) &
          0,usdbavg(itest,jtest),vsdbavg(itest,jtest), &
          0.0,p(itest,jtest,kk+1)*qonem 
          pzb = 0.0
          do k= 1,kk
            pzt = pzb
            pzb = min(pzt+dp(itest,jtest,k,n)*qonem, &
                          p(itest,jtest,kk+1)*qonem)
            write (lp,104) &
            k, &
            usd(itest,jtest,k),vsd(itest,jtest,k), &
            pzb-pzt,pzb
            if     (pzt.eq.p(itest,jtest,kk+1)*qonem) then
              exit
            endif
          enddo !k
        endif !test
      endif !debug
!
!   Now ensure  all Stokes Drift Velocity Fields are properly defined on Halos
!
      call xctilr(usd(    1-nbdy,1-nbdy,1),1,kk, 6,6, halo_uv)
      call xctilr(vsd(    1-nbdy,1-nbdy,1),1,kk, 6,6, halo_vv)
      call xctilr(usdbavg(1-nbdy,1-nbdy)  ,1, 1, 6,6, halo_uv)
      call xctilr(vsdbavg(1-nbdy,1-nbdy)  ,1, 1, 6,6, halo_vv)
      return
      end subroutine stokes_forfun

      subroutine stokes_vertical_j(n,j)
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
!
      integer   n,j
!
! --- --------------------------------------------
! --- interpolate Stokes in vertical, single j-row
! --- --------------------------------------------
!
      integer i,k
      real    dpthin
      logical lcm(nsdzi)      !use PCM for some layers?
      real    s1d(nsdzi,2),    & !input Stokes fields
              f1d(kdm,  2),    & !model Stokes fields
              c1d(nsdzi,2),    & !interpolation coefficients
              dpi( nsdzi),     & !input layer thicknesses, >= dpthin
              dprs(nsdzi),     & !input layer thicknesses
              pres(nsdzi+1),   & !input layer interfaces
              prsf(kdm+1)     !model layer interfaces
!
      dpthin = 0.001*onemm
      lcm(:) = .false.
!
      do i= 1,ii
        if (SEA_P) then
! ---     1-d existing arrays
          pres(1)=0.0
          do k=1,nsdzi
            s1d(k,1) = ws0*usdz(i,j,1,k) + &
                       ws1*usdz(i,j,2,k)
            s1d(k,2) = ws0*vsdz(i,j,1,k) + &
                       ws1*vsdz(i,j,2,k)
            pres(k+1)=min(sdzi(k),p(i,j,kk+1))
            dprs(k)  =pres(k+1)-pres(k)
            dpi( k)  =max(dprs(k),dpthin)
          enddo !k 1:nsdzi
          prsf(1)=0.0
          do k=1,kk
            prsf(k+1) = min(prsf(k)+dp(i,j,k,n),p(i,j,kk+1))
          enddo !k 1:kk
! ---     remap
          call hybgen_plm_coefs(s1d,     dpi, lcm,c1d, &
                                         nsdzi,   2,dpthin)
          call hybgen_plm_remap(s1d,pres,dprs,    c1d, &
                                f1d,prsf,nsdzi,kk,2,dpthin)
! ---     1-d new arrays
          do k=1,kk
            usdp(i,j,k) = f1d(k,1)
            vsdp(i,j,k) = f1d(k,2)
          enddo !k 1:kk
        endif !ip
      enddo !i
      return
      end subroutine stokes_vertical_j
#if defined(NERSC_stokes)
!     ><<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<
!   !  Construct stokes profile using Phillips approx.,
!      Breivik et al., 2016. 
!      Need on the Stokes transport and surface velocity
!     ><<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<< 
      subroutine stokes_vertical_j_phil(n,j)
      use mod_xc                ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
!
      integer   n,j
!
      integer i,k,l
      real    knm(2), sum0   !stokes depth  
      real    ,parameter :: betaa=1.0
      real    zbot,ztop, Vzbot,Vztop    !stokes depth  
      logical hycom_isnaninf
      real    ts1d(2)     !input Stokes trans. 
      real    f1dphil(kdm,2),    &  !Phillips actua. prof. 
              f1dphilavg(kdm,2), &  !Phillips averaged prof.
              prsf(kdm+1)           !model layer interfaces
!
      real  v0sp,tvsp,thetatv0,kbar,v0spzi, &
           Vztopsp,Vzbotsp,v0spzilayer 
      
! --- ---------------------------------------------------------
! --- Construct Stokes in vertical using Phillips, single j-row
! --- ---------------------------------------------------------
!
         do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))   
!     --- hycom layers
!$2DIM              do k=1,nsdzi   ! k looooooooop
!$2DIM                 s1d(k,1) = ws0*usdz(i,j,1,k) +
!$2DIM      &                     ws1*usdz(i,j,2,k)
!$2DIM                 s1d(k,2) = ws0*vsdz(i,j,1,k) +
!$2DIM      &                     ws1*vsdz(i,j,2,k)
!$2DIM               enddo             !k 1:nsdzi  
!    --- input transport at i,j (time interp.)
              ts1d(1) = ws0*tusd(i,j,1) + &
                        ws1*tusd(i,j,2)
              ts1d(2) = ws0*tvsd(i,j,1) + &
                        ws1*tvsd(i,j,2)
!    ---  comp speed and angle of the transport
!    ---  comp stokes surface speed              
              v0sp=sqrt(usds(i,j)**2+vsds(i,j)**2) 
              tvsp=sqrt(ts1d(1)**2+ts1d(2)**2)
              thetatv0=atan2(ts1d(2),ts1d(1))
!    ---  comp Stokes depth  
              kbar=v0sp/(2.0*max(tvsp,epsil)) &
                  *(1.0-2.0*betaa/3.0)
!    --  Now comp layer averaged profile
         if     (debug_stokesPhl) then    
          if( itest.gt. 0 .and. jtest.gt. 0)then
  103      format (i9,2i5,a)
            write (lp,103) nstep,itest+i0,jtest+j0, &
            ' stokes_forfun: Phlavgsp  Phlsp  usdz      dpth'
          endif  
         endif
         sum0=0.0
         prsf(1)=0.0
         do k=1,kk
!            prsf(k+1) = p(i,j,k+1)
            prsf(k+1)=min(prsf(k)+dp(i,j,k,n),p(i,j,kk+1))
            !   layer k:  zbot,ztop
            ztop=prsf(k)/onem    ! top    
            zbot=prsf(k+1)/onem  ! bot
            ! 
!   --   actual Stokes profile (not really needed, just for debugging)                   
            v0spzi=phillips_profile(v0sp, &
                kbar,-1*ztop,betaa)                   
            f1dphil(k,1)=v0spzi*cos(thetatv0)
            f1dphil(k,2)=v0spzi*sin(thetatv0)
!   --      layer_averaged Stokes profile               
            Vztopsp=phillips_transportz(v0sp, &
                                    kbar,-1*ztop,betaa)
            Vzbotsp=phillips_transportz(v0sp, &
                                    kbar,-1*zbot,betaa)
            v0spzilayer=(Vzbotsp-Vztopsp)/max(zbot-ztop,epsil)
            f1dphilavg(k,1)=v0spzilayer*cos(thetatv0)
            f1dphilavg(k,2)=v0spzilayer*sin(thetatv0)
!     ------ Finish   computing layeravg stokes
!            now update assign stokes at p points
            usdp(i,j,k) = f1dphilavg(k,1)
            vsdp(i,j,k) = f1dphilavg(k,2)
!   >---------------------------------------------------------
!$$        if     (debug_stokes) then
!$$  103    format (i9,2i5,a)
!$$  104    format (30x,i3,2f8.4,f9.3,f9.2)
!$$         if (itest.gt.0 .and. jtest.gt.0) then
!$$           write (lp,103) nstep,itest+i0,jtest+j0, &
!$$           '  stokes_forfun:  usdz    vsdz    thkns     dpth'
!$$           pzb = 0.0
!$$           do k= 1,nsdzi
!$$             pzt = pzb
!$$             pzb = min(sdzi(k),p(itest,jtest,kk+1))*qonem
!$$             write (lp,104)                                       &
!$$             k,                                                   &
!$$             ws0*usdz(itest,jtest,1,k)+ws1*usdz(itest,jtest,2,k), &
!$$             ws0*vsdz(itest,jtest,1,k)+ws1*vsdz(itest,jtest,2,k), &
!$$             pzb-pzt,pzb
!$$             if     (pzt.eq.p(itest,jtest,kk+1)*qonem) then
!$$               exit
!$$             endif
!$$ C$$           enddo !k
!$$       endif
        if     (debug_stokesPhl) then        
        if((i+i0).eq. itest .and. (j+j0).eq. jtest)then
!  103    format (i9,2i5,a)
  104    format (30x,i3,3f8.4,f9.2)
  105    format (30x,i3,f8.4,f9.2)         
!           write (lp,103) nstep,itest+i0,jtest+j0,
!          '  stokes_forfun:  Philavgsp   Philsp  usdz    dpth'
             sum0 = sum0                                        &
                   + sqrt(f1dphilavg(k,1)**2                    &
                          +f1dphilavg(k,2)**2)*dp(i,j,k,n)/onem
           if (k.le.nsdzi)then 
             write (lp,104)                                &
             k,                                            &
             sqrt(f1dphilavg(k,1)**2+ f1dphilavg(k,2)**2), &
             sqrt(f1dphil(k,1)**2+ f1dphil(k,2)**2),       &
             v0sp,                                         &
             ztop  
           else
              write (lp,105)                               &
              k,                                           &
              sqrt(f1dphilavg(k,1)**2+ f1dphilavg(k,2)**2),&
              ztop 
              if     (ztop.eq.p(itest,jtest,kk+1)*qonem) then
                  exit
              endif
             endif
           endif
        endif !debugstokes
!     >-------------------------------------------------------            
          enddo  !k 1:kk
          if  (debug_stokesPhl) then                
            if((i+i0).eq. itest .and. (j+j0).eq. jtest)then
              write(lp,'(a,2E16.8)')                        &
                           '-->>> T_phil_avg=,T_in  ', sum0 &
                 ,sqrt(ts1d(1)**2 +ts1d(2)**2 )         
            endif
           endif
          enddo  !i
       enddo   !l

!$$$c        !##############################################Phillips
!$$$cy c     !##############################################Phillips     
      return
      end subroutine stokes_vertical_j_phil
!      
      real function phillips_transportz(v0sur,knm,zlvl,bta)
!  ---   comp stokes transport between surface and zlevel 
        IMPLICIT NONE
        REAL :: v0sur, knm, zlvl, bta
        real,    parameter :: epsil=1.e-20
!     REAL             :: s
        phillips_transportz=v0sur/(2*max(knm,epsil))     &
            *(1 - exp(2*knm*zlvl)                        &
            - (2*bta/3)*(1+sqrt(pi)*(-2*knm*zlvl)**(1.5) &
            *get_erfc(sqrt(-2*knm*zlvl)) -               &
            (1 -2*knm*zlvl)*exp(2*knm*zlvl)))
        
      END FUNCTION phillips_transportz
      
      real function phillips_profile(v0surp,knmp,zlvlp,btap)
!  ---   comp stokes profile based on Philips       
      IMPLICIT NONE
      REAL :: v0surp, knmp, zlvlp, btap      
      phillips_profile=v0surp*(exp(2*knmp*zlvlp) &
           - btap*sqrt(-2*pi*knmp*zlvlp)         &
           *get_erfc(sqrt(-2*knmp*zlvlp)))
!     print *,'z, phill_prof =',zlvlp,phillips_profile     
      END FUNCTION phillips_profile
#endif
!
!      
      end module mod_stokes
!> May  2014 - use land/sea masks (e.g. ip) to skip land
!> Aug  2015 - added stdarc

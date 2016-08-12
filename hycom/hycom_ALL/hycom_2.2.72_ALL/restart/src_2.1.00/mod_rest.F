      module mod_rest
      use mod_xc
      implicit none
c
      integer, save :: kdm
c
c-----------------------------------------------------------------------------
c
      real,    save, allocatable, dimension (:,:) :: 
     & depths         ! water depth

      real,    save, allocatable, dimension (:,:,:,:) :: 
     & u,v,           ! velocity components
     & dp,            ! layer thickness
     & temp,          ! temperature
     & saln,          ! salinity
     & th3d           ! potential density

      real,    save, allocatable, dimension (:,:,:) :: 
     & ubavg,vbavg,   ! barotropic velocity
     & pbavg          ! barotropic pressure

      real,    save, allocatable, dimension (:,:) :: 
     & pbot           ! bottom pressure at t=0
     & psikk,         ! montg.pot. in bottom layer
     & thkk           ! virtual potential density in bottom layer

      real,    save, allocatable, dimension (:,:,:) :: 
     & dpmixl,        ! mixed layer depth
     & thmix          ! mixed layer density (MICOM only)

      real,    save, allocatable, dimension (:,:) :: 
     & thkice,        ! grid-cell avg. ice thknss (m)
     & covice,        ! ice coverage (rel.units)
     & temice         ! ice surf.temp.

      real,    save, allocatable, dimension (:,:,:) :: 
     & tracer         ! inert tracer (optional)

      integer, save, allocatable, dimension (:,:) :: 
     & ip,         ! pressure   grid land/sea mask
     & iu,         ! u-velocity grid land/sea mask
     & iv          ! v-velocity grid land/sea mask

      real,    save, allocatable, dimension (:) :: 
     & sigma    ! layer densities (MICOM only)

      real ::
     & thbase,  ! reference density (MICOM only)
     & baclin,  ! baroclinic time step
     & huge     ! large number used to indicate land points

      common/parms1r/ thbase,baclin,huge
      save  /parms1r/

      integer ::
     & yrflag,  ! days in year flag (0=360,1=366)
     & iversn,  ! hycom version number x10
     & iexpt    ! experiment number x10

      common/parms1i/ yrflag,iversn,iexpt
      save  /parms1i/

      logical ::
     & icegln,  ! use energy loan ice model
     & trcout   ! advect tracer and save results in history/restart file

      common/swtchs/ icegln,trcout
      save  /swtchs/

      character*48 ::
     & flnmdep,  !  input bathy   filename
     & flnmrsa,  !  input archive filename
     & flnmrsi,  !  input restart filename
     & flnmrso   ! output restart filename

      common/iovars/ flnmdep,flnmrsa,flnmrsi,flnmrso
      save  /iovars/
c-----------------------------------------------------------------------------
c
c --- module subroutines
c
      contains

      subroutine rest_alloc
      implicit none
c
c --- initialize allocatable arrays.
c
      allocate( depths(idm,jdm) )
c
      allocate(      u(idm,jdm,kdm,2) )
      allocate(      v(idm,jdm,kdm,2) )
      allocate(     dp(idm,jdm,kdm,2) )
      allocate(   temp(idm,jdm,kdm,2) )
      allocate(   saln(idm,jdm,kdm,2) )
      allocate(   th3d(idm,jdm,kdm,2) )
      allocate( tracer(idm,jdm,kdm,2) )
c
      allocate(  ubavg(idm,jdm,3) )
      allocate(  vbavg(idm,jdm,3) )
      allocate(  pbavg(idm,jdm,3) )
c
      allocate(   pbot(idm,jdm) )
      allocate(  psikk(idm,jdm) )
      allocate(   thkk(idm,jdm) )
c
      allocate( dpmixl(idm,jdm,2) )
      allocate(  thmix(idm,jdm,2) )
c
      allocate( thkice(idm,jdm) )
      allocate( covice(idm,jdm) )
      allocate( temice(idm,jdm) )
c
      allocate(     ip(idm,jdm) )
      allocate(     iu(idm,jdm) )
      allocate(     iv(idm,jdm) )
c
      allocate(  sigma(kdm) )

      end subroutine rest_alloc

      end module mod_rest

      program grid
      use mod_xc 
      use mod_za
      use mod_sigma
      use m_parse_blkdat
      ! Planed to use this as illustration of idealized forcing setup. 
      ! TODO: Make a illustration of how to proceed
      implicit none

      real, parameter :: onem=9806.0
      real, parameter :: radian= 57.2957795

      character(len=80) gridid, cline, flnmdep, flnmrso
      character(len=11) tag7

      real, allocatable, dimension(:,:)   ::  depths
      real, allocatable, dimension(:)     ::  theta

      integer, allocatable, dimension(:,:) ::  ip

      real :: centerlon, centerlat, dx, dy
      real :: pi, rad, deg
      real :: xmin, xmax, rdummy
      integer :: idummy
      integer :: i,j,k,l

      integer :: kapflg, thflag
      real    :: thbase, dtime
      integer :: kdm
      integer :: yrflag, iexpt, iversn, nstep=0


      flnmdep='regional.depth'
      dtime=0.0
      iexpt=1
      iversn=21
      nstep=0




      ! First tests -- hardcoded grid size etc
      idm=400
      jdm=200
      call xcspmd()
      call zaiost()

      if (idm/=idm .or. jdm/=jdm) then
         print *,'Grid size inconsistency - i: ',idm,idm
         print *,'Grid size inconsistency - j: ',jdm,jdm
         stop '()'
      end if


      allocate(depths(idm,jdm))
      allocate(ip    (idm,jdm))
      ip=0

      ! Get model bathymetry
      call zaiopf(flnmdep(1:len_trim(flnmdep))//'.a','old', 12)
      call zaiord(depths, ip,.false., xmin,xmax,  12)
      call zaiocl(12)


      !TODO Bring in examples from synoptic-forcing here...

      end

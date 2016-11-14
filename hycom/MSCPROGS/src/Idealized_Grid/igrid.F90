      program grid
      use mod_xc 
      use mod_za
      use m_griddump
      use m_parse_blkdat
! --- this program can be used to create simple idealized grid
! --- configurations. Currently it takes two configurations, seamount
! --- or flat bottom. Grid size is read from blkdat.input. Output is
! --- topography files (regional.grid/depth) which can be understood by 
! --- hycom

! --- For other configurations than flat bottom/sea mount, you must modify
! --- this file ( see depths generation below ) . Make sure you close the 
! --- southern and northern boundaries (j=1 and j=jdm). The eastern/western
! --- boundary can be open ( the grid will be interpreted as periodic in 
! --- that case )

! --- NB: For hycom 2.2 and standard hycom, remember to set depths to huge where 
! ---     depths=0 (before saving to regional.depth). This is done in routine 
! ---     m_griddump
      implicit none

      character(len=80) gridid, bathyfile, tmparg, mode
      character(len=11) tag7

      real, allocatable, dimension(:,:) ::  modlat, modlon, &
         plat,plon, ulat,ulon, vlat,vlon, qlat,qlon, depths, &
         scpx, scpy, scux, scuy, scvx, scvy, scqx, scqy, asp, pang, &
         corio
      real*8, allocatable, dimension(:,:) ::  iodepths,iolon,iolat

      real :: centerlon, centerlat, dx, dy
      real :: pi, rad, deg, rdummy
      integer :: mapflg
      integer :: i,j
#if defined (IARGC)
      integer*4, external :: iargc
#endif



      if (iargc()==4) then
         print *,'Routine for creating simplified hycom topography, grid files'
         print *,'Usage:'
         print *,' igrid <topo option> <grid size> <centerlon> <centerlat>'
         print *
         print *,'First argument:'
         print *,'seamount  for seamount bathymetry placed mid-basin'
         print *,'flat      for flat bottom'

         print *,'Second argument is grid size (uniform) in meters' 
         print *,'Third  argument is longitude (also uniform) '
         print *,'Fourth argument is latitude (also uniform) '
         print *
         print *,'Grid size will be retrieved from blkdat.input'
         print *,'Routine produces regional.grid/depth files which'
         print *,'are understood by hycom'
         stop '()'
      end if


      ! First check if topo files are already present



      call getarg(1,mode)
      call getarg(2,tmparg); read(tmparg,*) dx; dy=dx
      call getarg(3,tmparg); read(tmparg,*) centerlon
      call getarg(4,tmparg); read(tmparg,*) centerlat



      ! Frst tests -- hardcoded grid size etc
      print *,'NB: grid size retrieved from blkdat.input'
      call parse_blkdat('idm   ','integer',rdummy,idm)
      call parse_blkdat('jdm   ','integer',rdummy,jdm)
      call zaiost()

      pi=4.*atan(1.)
      rad=pi/180.
      deg=180./pi
      mapflg=4
      allocate (modlat(idm,jdm),modlon(idm,jdm))
      allocate (plat(idm,jdm),plon(idm,jdm))
      allocate (qlat(idm,jdm),qlon(idm,jdm))
      allocate (ulat(idm,jdm),ulon(idm,jdm))
      allocate (vlat(idm,jdm),vlon(idm,jdm))

      allocate(scpx(idm,jdm),scpy(idm,jdm))
      allocate(scux(idm,jdm),scuy(idm,jdm))
      allocate(scvx(idm,jdm),scvy(idm,jdm))
      allocate(scqx(idm,jdm),scqy(idm,jdm))

      allocate (asp (idm,jdm),pang(idm,jdm))
      allocate (depths(idm,jdm))
      allocate (corio (idm,jdm))
      allocate (iodepths(idm,jdm))
      allocate (iolon(idm,jdm))
      allocate (iolat(idm,jdm))

      ! Constant lon/lat and grid size
      ulon=centerlon
      vlon=centerlon
      plon=centerlon
      qlon=centerlon

      ulat=centerlat
      vlat=centerlat
      plat=centerlat
      qlat=centerlat

      corio=sin(centerlat)*4.*pi/86400.
      print *,rad,pi,maxval(corio)



      scpx=dx;
      scpy=dx;
      scux=dx;
      scuy=dx;
      scvx=dx;
      scvy=dx;
      scqx=dx;
      scqy=dx;



      modlat(1:idm,1:jdm)=plat(1:idm,1:jdm)
      modlon(1:idm,1:jdm)=plon(1:idm,1:jdm)
      iolon=modlon
      iolat=modlat
      open(10,file='newpos.uf',form='unformatted',status='unknown')
         write(10)iolat,iolon
      close(10)


      ! Added by Gisle
      ! Seamount case according to Berntsen et al. (2006).
      if (trim(mode)=='seamount') then
         depths(:,1:2)=0.
         depths(:,jdm-1:jdm)=0.
         do i=1,idm
         do j=3,jdm-2
             depths(i,j)=5000.*(1.-0.9*EXP((-(((idm/2.)-i)**2+((jdm/2.)-j)**2))/(40000./dx)**2))
         end do
         end do
         depths(1:2,:)=0.
         depths(idm-1:idm,:)=0.
      elseif (trim(mode)=='flat') then

         ! depths - flat bottom
         depths=500.
         depths(:,1:3   )=0.
         depths(:,jdm-3:jdm)=0.

      else

         print *,'Input is '//trim(mode)
         print *,'seamount  for seamount bathymetry placed mid-basin'
         print *,'flat      for flat bottom'
         stop '()'
      end if





      print *,'Creating regional.depth.[ab]'
      print *,'Creating regional.grid.[ab]'
      call griddump(plon,plat,ulon,ulat,vlon,vlat,qlon,qlat, &
         scpx,scpy,scux,scuy,scvx,scvy,scqx,scqy,corio,asp,pang,depths, &
         mapflg)

      iodepths=depths
      if (idm>999 .or. jdm > 999) then
         write(tag7,'(i5.5,a,i5.5)')idm,'x',jdm
      else
         write(tag7,'(i3.3,a,i3.3)')idm,'x',jdm
      end if
      bathyfile='depths'//trim(tag7)//'.uf'
      open(unit=10,file=trim(bathyfile),status='replace',form='unformatted')
      write(10)iodepths
      close(10)





! Write latlon.dat
!      gridid='CONFORMAL GRID'
!      open(10,file='latlon.dat',form='formatted')
!      write(10,'(2i5)')idm,jdm
!      write(10,'(a80)')gridid
!      write(10,'(15e14.7)')((qlat(i,j),i=0,idm+1),j=0,jdm+1)
!      write(10,'(15e14.7)')((qlon(i,j),i=0,idm+1),j=0,jdm+1)
!      write(10,'(15e14.7)')((plat(i,j),i=0,idm+1),j=0,jdm+1)
!      write(10,'(15e14.7)')((plon(i,j),i=0,idm+1),j=0,jdm+1)
!      write(10,'(15e14.7)')((ulat(i,j),i=0,idm+1),j=0,jdm+1)
!      write(10,'(15e14.7)')((ulon(i,j),i=0,idm+1),j=0,jdm+1)
!      write(10,'(15e14.7)')((vlat(i,j),i=0,idm+1),j=0,jdm+1)
!      write(10,'(15e14.7)')((vlon(i,j),i=0,idm+1),j=0,jdm+1)
!      close(10)




      print *,'Max dx = ',dx
      print *,'Min dx = ',dx


      stop '(normal)'

      end

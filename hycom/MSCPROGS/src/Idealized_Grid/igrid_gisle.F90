      program grid
      use mod_xc 
      use mod_za
      use mod_hycom_fileio
! --- this program create a grid as specified in the file grid.info

! --- --------------------------------------------------------------------------------
! --- KAL: 20060616 -- Changed p/q/u/v arrays from real to real*8. A lot
! ---                  of the steps in the transformation can lead to 
! ---                  precision loss, so it is best to keep real 8 for these
! ---                  Also, most variables in oldtonew and newtoold, and
! ---                  in confmap are now real*8 for the same reason.
! --- 
! ---                  NB: Before this, the precision loss mainly occured when calculating
! ---                  grid sizes. lon/lat values looked ok, but miniscule
! ---                  differences affected distance between lon/lat points.
! --- --------------------------------------------------------------------------------

      implicit none

      character(len=80) gridid, bathyfile
      character(len=11) tag7

      real, allocatable, dimension(:,:) ::  modlat, modlon, &
         plat,plon, ulat,ulon, vlat,vlon, qlat,qlon, depths, &
         scpx, scpy, scux, scuy, scvx, scvy, scqx, scqy, asp, pang, &
         corio
      real*8, allocatable, dimension(:,:) ::  iodepths,iolon,iolat

      real :: centerlon, centerlat, dx, dy
      real :: pi, rad, deg
      integer :: mapflg
      integer :: i,j

      !real, parameter :: radian= 57.2957795
      !real, parameter :: pi    = 3.1415926536




      ! Frst tests -- hardcoded grid size etc
      idm=50
      jdm=54
      call xcspmd()
      call zaiost()

      if (idm/=idm .or. jdm/=jdm) then
         print *,'Grid size inconsistency - i: ',idm,idm
         print *,'Grid size inconsistency - j: ',jdm,jdm
         stop '()'
      end if


      centerlon=-18.
      centerlat=-45.
      dx=20000./3.;
      dy=20000./3.;
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

      corio=sind(centerlat)*4.*pi/86400.
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

      !depths when using flat bottom
      depths(:,1:2)=0.
      depths(:,3:52)=5000.
      depths(:,53:54)=0.

      ! Sloping topography (Slørdal PhD thesis 1995), closed western/eastern boundaries
      !depths(1:2,:)=0.
      !depths(3,:)=3000.
      !do i=4,idm-2
      !   depths(i,:)=210.-(3000.-210.)/(idm-1)*(i-(idm-2))
      !end do
      !depths(idm-1:idm,:)=0.


      ! Sloping topography (Slørdal PhD thesis 1995), closed northern/southern boundaries
      !depths(:,1:2)=0.
      !depths(:,3)=3000.
      !do j=4,jdm-2
      !   depths(:,j)=210.-(3000.-210.)/(jdm-1)*(j-(jdm-2))
      !end do
      !depths(:,jdm-1:jdm)=0.

      ! Seamount case according to Berntsen et al. (2006).
      !depths(:,1:2)=0.
      !do i=1,idm
      !  do j=3,jdm-2
      !    depths(i,j)=5000.*(1.-0.9*EXP((-((25.-i)**2+(27.-j)**2))/(40000./dx)**2))
      !end do
      !  end do
      !depths(:,jdm-1:jdm)=0.
      !print*,'depths',depths 

      print *,'Creating regional.depth.[ab]'
      print *,'Creating regional.grid.[ab]'
      call dump_regional(plon,plat,ulon,ulat,vlon,vlat,qlon,qlat, &
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

      !call checktopo(depths,idm,jdm)

      ! Check if depth is non zero along boundaries
      if (any(depths(1,:) > .5)  .or. any(depths(idm,:) > .5)  ) then
         print *
         print *
         print *,'NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB'
         print *,'Depths along eastern/western boundaries are non-zero'
         print *,'This is ok for a periodic model domain, but not else..'
         print *,'Insert the following lines in grid.topofix to fix this:'
         print '(4i5,f4.1)',1,1,1,jdm,0.
         print '(4i5,f4.1)',idm,idm,1,jdm,0.
      end if

      if (any(depths(:,1)>.5)  .or. any(depths(:,jdm)>.5)  ) then
         print *
         print *
         print *,'NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB'
         print *,'Depths along northern/southern boundaries are non-zero'
         print *,'This is NOT ok'
         print *,'Insert the following lines in grid.topofix to fix this:'
         print '(4i5,f4.1)',1,idm,1,1,0.
         print '(4i5,f4.1)',1,idm,jdm,jdm,0.
      end if

      stop '(normal)'

      end

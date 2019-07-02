   program grid
   use mod_xc , only : xcspmd, idm , jdm
   use mod_confmap
   use m_bigrid
   use m_topography
   use m_topofix
   use m_tecconfgrid
   use m_spherdist
   use m_checktopo
   use m_grid_consistency
   use m_grid_to_hycom
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

   real maxlat,maxlon,minlon,minlat
#if defined (IARGC)
   integer*4, external :: iargc
#endif
   real lat_n,lon_n,lat_o,lon_o,lat_r,lon_r
   logical smooth
   integer ish,ifilt
   integer, parameter :: shdim=16
   real sh(0:shdim)
   integer i,j

   character(len=80) gridid
   character(len=11) tag7

   integer, parameter :: ms=20

   real, allocatable ::  modlat(:,:),modlon(:,:)
   real*8, allocatable ::  plat(:,:),plon(:,:)
   real*8, allocatable ::  ulat(:,:),ulon(:,:)
   real*8, allocatable ::  vlat(:,:),vlon(:,:)
   real*8, allocatable ::  qlat(:,:),qlon(:,:)
   real, allocatable ::  dbodc(:,:)
   real, allocatable ::  dibcao(:,:)
   real, allocatable ::  detopo5(:,:)
   real, allocatable ::  dgebco(:,:)
   real, allocatable ::  dconman(:,:)
   real, allocatable ::  depths(:,:)
   real*8, allocatable ::  iodepths(:,:),iolon(:,:),iolat(:,:)
   real, allocatable ::  dx(:,:)
   real, allocatable ::  dy(:,:)
   logical, allocatable ::  mask(:,:)

   integer, allocatable :: ip(:,:)
   integer, allocatable :: ifp(:,:),ilp(:,:),isp(:)
   integer, allocatable :: jfp(:,:),jlp(:,:),jsp(:)
   character(len=99) :: path



   ! initialize conformal map. When create=true, the two first
   ! arguments are RETURNED from the routine
   ! NB: setting idm and jdm directly removes the need for
   ! calling xcspmd()
   call initconfmap(idm,jdm,createin=.true.)

   ! Dumping dimension.h. Not really needed anymore , but just in case
   open(10, file='dimension.h')
   write(10,*)'     INTEGER idm,jdm,izdm,iwmax,nrmonths'
   write(10,'(a,i5,a,i5,a)')'      PARAMETER(idm=',idm,', jdm=', &
      jdm,', izdm=12, iwmax=18)'
   write(10,'(a)')'      PARAMETER(nrmonths=izdm) '// &
                  ' ! do not edit this line (nz is 12 months)'
   close(10)

   allocate (modlat(idm,jdm),modlon(idm,jdm))
   allocate (plat(0:idm+1,0:jdm+1),plon(0:idm+1,0:jdm+1))
   allocate (qlat(0:idm+1,0:jdm+1),qlon(0:idm+1,0:jdm+1))
   allocate (ulat(0:idm+1,0:jdm+1),ulon(0:idm+1,0:jdm+1))
   allocate (vlat(0:idm+1,0:jdm+1),vlon(0:idm+1,0:jdm+1))
   allocate (dbodc(idm,jdm))
   allocate (dibcao(idm,jdm))
   allocate (dgebco(idm,jdm))
   allocate (dconman(idm,jdm))
   allocate (depths(idm,jdm))
   allocate (iodepths(idm,jdm))
   allocate (iolon(idm,jdm))
   allocate (iolat(idm,jdm))
   allocate (detopo5(idm,jdm))
   allocate (dx(idm,jdm))
   allocate (dy(idm,jdm))
   allocate (mask(idm,jdm))
   allocate(ip(idm,jdm))
   allocate(ifp(jdm,ms),ilp(jdm,ms),isp(jdm))
   allocate(jfp(idm,ms),jlp(idm,ms),jsp(idm))

! --- calculate positions for gridpoint in new projection using
! --- specified gridsize

   do j=0,jdm+1
   do i=0,idm+1
!         lon_n=wlim+(i-1)*di
!         lat_n=slim+(j-1)*dj
!         if (mercator) then
!            lat_n=slim+(j-1)*dm
!            !print *,i,j,'lat_n',lat_n
!            lat_n=(2.*atan(exp((lat_n*rad)))-pi*.5)*deg
!            !print *,i,j,'lat_n',lat_n
!         endif
!         call newtoold(lat_n,lon_n,lat_o,lon_o)
      call getgridP(i,j,lat_o,lon_o)
      plat(i,j)=lat_o
      plon(i,j)=lon_o


! --- --- for testing of the transformation

!         call oldtonew(lat_o,lon_o,lat_r,lon_r)
!         if ((abs(lat_r-lat_n).gt.0.001).or. &
!             (abs(lon_r-lon_n).gt.0.001)) then
!           write(2,*) i,j
!           write(2,*) lat_n,lon_n
!           write(2,*) lat_r,lon_r
!         endif

      if (.true.) then
! U point
!            lon_n=wlim+(float(i)-1.0-0.5)*di
!            lat_n=slim+(float(j)-1.0    )*dj
!            if (mercator) then
!               lat_n=slim+(float(j)-1.0   )*dm
!               lat_n=(2.*atan(exp((lat_n*rad)))-pi*.5)*deg
!            endif
!            call newtoold(lat_n,lon_n,lat_o,lon_o)
         call getgridU(i,j,lat_o,lon_o)
         ulat(i,j)=lat_o
         ulon(i,j)=lon_o
! V point
!            lon_n=wlim+(float(i)-1.0)*di
!            lat_n=slim+(float(j)-1.0-0.5)*dj
!            if (mercator) then
!               lat_n=slim+(float(j)-1.0-0.5)*dm
!               lat_n=(2.*atan(exp((lat_n*rad)))-pi*.5)*deg
!            endif
!            call newtoold(lat_n,lon_n,lat_o,lon_o)
         call getgridV(i,j,lat_o,lon_o)
         vlat(i,j)=lat_o
         vlon(i,j)=lon_o
! q point
!            lon_n=wlim+(float(i)-1.0-0.5)*di
!            lat_n=slim+(float(j)-1.0-0.5)*dj
!            if (mercator) then
!               lat_n=slim+(float(j)-1.0-0.5)*dm
!               lat_n=(2.*atan(exp((lat_n*rad)))-pi*.5)*deg
!            endif
!            call newtoold(lat_n,lon_n,lat_o,lon_o)
         call getgridQ(i,j,lat_o,lon_o)
         qlat(i,j)=lat_o
         qlon(i,j)=lon_o
      endif
   enddo
   enddo

   do j=1,jdm
      do i=1,idm
        dy(i,j)=0.5*0.001*(&
        spherdist8(plon(i,j-1),plat(i,j-1),plon(i,j),plat(i,j))      &
       +spherdist8(plon(i-1,j-1),plat(i-1,j-1),plon(i-1,j),plat(i-1,j)))

        dx(i,j)=0.5*0.001*(&
        spherdist8(plon(i-1,j),plat(i-1,j),plon(i,j),plat(i,j))      &
       +spherdist8(plon(i-1,j-1),plat(i-1,j-1),plon(i,j-1),plat(i,j-1)))
      enddo
   enddo

! modlat modlon generation
   modlat(1:idm,1:jdm)=plat(1:idm,1:jdm)
   modlon(1:idm,1:jdm)=plon(1:idm,1:jdm)
   iolon=modlon
   iolat=modlat
   open(10,file='newpos.uf',form='unformatted',status='unknown')
      write(10)iolat,iolon
   close(10)


! Saving the new lon,lat coordinates
!      open(10,file='lonN.dat')
!         do i=0,idm+1
!            write(10,'(f13.5,i5,f13.5)')di,i,wlim+(i-1)*di
!         enddo
!      close(10)
!      open(10,file='latN.dat')
!         do j=0,jdm+1
!            if (mercator) then
!               lat_n=slim+(j-1)*dm
!               lat_n=(2.*atan(exp((lat_n*rad)))-pi*.5)*deg
!            else
!               lat_n=slim+(j-1)*dj
!            endif
!            write(10,'(f13.5,i5,f13.5)')dm,j,lat_n
!         enddo
!      close(10)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Layout of model grid
   minlon=minval(modlon); print *,'Minimum longitude is: ',minlon
   maxlon=maxval(modlon); print *,'Maximum longitude is: ',maxlon
   minlat=minval(modlat); print *,'Minimum latitude  is: ',minlat
   maxlat=maxval(modlat); print *,'Maximum latitude  is: ',maxlat
   print '(a,f8.2,a)'    ,'     #--------- ',maxlat,' -----------#'
   print '(a)'           ,'     |                              |'
   print '(a)'           ,'     |                              |'
   print '(a)'           ,'     |                              |'
   print '(a)'           ,'     |                              |'
   print '(a)'           ,'     |                              |'
   print '(f8.2,a,f8.2)' ,minlon,'                       ',maxlon
   print '(a)'           ,'     |                              |'
   print '(a)'           ,'     |                              |'
   print '(a)'           ,'     |                              |'
   print '(a)'           ,'     |                              |'
   print '(a)'           ,'     |                              |'
   print '(a,f8.2,a)'    ,'     #--------- ',minlat,' -----------#'


! --- may skip the topography
   depths=0.0
   if (confmap_calctopo()) then 
      call topography(depths,modlon,modlat,dbodc,detopo5,dibcao,dgebco,dconman, &
                      idm,jdm,trim(path))
      print *,'topography done'
      call bigrid(depths,idm,jdm)
      print *,'bigrid done'
      !call grid_consistency(depths,idm,jdm,80)
      !print *,'grid_consistency done'
      call topofix(depths,idm,jdm)
      print *,'topofix done'
      call bigrid(depths,idm,jdm)
      print *,'bigrid done'

!         if (smooth) then
!            call index(etopo5,idm,jdm,ip,ifp,ilp,isp,jfp,jlp,jsp,ms)
!            call shapiro_ini(ish,sh,shdim)
!            call shapiro_filt(ish,sh,etopo5,idm,jdm,shdim,ifilt,&
!                              ifp,ilp,isp,jfp,jlp,jsp,ms)
!         endif

      ! Knut - simple fix for huge grids (idm or jdm > 999 )
      if (idm>999 .or. jdm > 999 ) then
         write(tag7,'(i5.5,a,i5.5)')idm,'x',jdm
      else
         write(tag7,'(i3.3,a,i3.3)')idm,'x',jdm
      end if
      iodepths=depths
      open(unit=10,file='depths'//trim(tag7)//'.uf',status='unknown',form='unformatted')
         write(10)iodepths
      close(10)

      call  tecconfgrid(depths,dx,dy,modlon,modlat,dbodc,detopo5,dibcao,dgebco,dconman,idm,jdm)

      ! mask setup
      mask=.false.
      where (depths > 0.0) mask=.true.
      open(10,file='mask.uf',status='unknown',form='unformatted')
         write(10)mask
      close(10)
      open(10,file='mask.asc',form='formatted')
      write(10,'(66i1)') ((min(1,int(depths(i,j))),j=1,jdm),i=1,idm)
      close(10)


      !KAL -- new - dump the regional.depths and regional.grid files
      print *,'Creating regional.depth.[ab]'
      print *,'Creating regional.grid.[ab]'
      call grid_to_hycom(plon,plat,qlon,qlat,ulon,ulat,vlon,vlat,&
                         depths)
   else
      print *,'topo flag false in grid.info - not dumping bathymetry'
   endif

    
! Dump boundary lines
   open(10,file='boundary.dat')
   j=1
   do i=1,idm
      write(10,'(2f10.4)')modlon(i,j),modlat(i,j)
   enddo
   i=idm
   do j=2,jdm
      write(10,'(2f10.4)')modlon(i,j),modlat(i,j)
   enddo
   j=jdm
   do i=idm-1,1,-1
      write(10,'(2f10.4)')modlon(i,j),modlat(i,j)
   enddo
   i=1
   do j=jdm-1,2,-1
      write(10,'(2f10.4)')modlon(i,j),modlat(i,j)
   enddo
   close(10)


   ! Write latlon.dat
   gridid='CONFORMAL GRID'
   open(10,file='latlon.dat',form='formatted')
   write(10,'(2i5)')idm,jdm
   write(10,'(a80)')gridid
   write(10,'(15e14.7)')((qlat(i,j),i=0,idm+1),j=0,jdm+1)
   write(10,'(15e14.7)')((qlon(i,j),i=0,idm+1),j=0,jdm+1)
   write(10,'(15e14.7)')((plat(i,j),i=0,idm+1),j=0,jdm+1)
   write(10,'(15e14.7)')((plon(i,j),i=0,idm+1),j=0,jdm+1)
   write(10,'(15e14.7)')((ulat(i,j),i=0,idm+1),j=0,jdm+1)
   write(10,'(15e14.7)')((ulon(i,j),i=0,idm+1),j=0,jdm+1)
   write(10,'(15e14.7)')((vlat(i,j),i=0,idm+1),j=0,jdm+1)
   write(10,'(15e14.7)')((vlon(i,j),i=0,idm+1),j=0,jdm+1)
   close(10)
   print *,'Max dx = ',max(maxval(dx),maxval(dy))
   print *,'Min dx = ',min(minval(dx),minval(dy))

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

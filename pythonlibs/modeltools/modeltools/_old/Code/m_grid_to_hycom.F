      module m_grid_to_hycom
      contains
      subroutine grid_to_hycom(plon,plat,qlon,qlat,ulon,ulat,vlon,vlat,
     &                         depths)
      use mod_xc
      use mod_za
      use m_spherdist
      implicit none


! August 2000 : modified for hycom 1.0.08
! July   2003 : modified for hycom 2.1.03
! Aug    2005 : part of conf_grid 
! --- --------------------------------------------------------
! --- Computes/load the Coriolis parameter and estimates the
! --- grid size parameters.
! --- 
! --- The Coriolis parameter is stored in corio(udm,jdm)
! --- scux:    zonal grid distance at u-point
! --- scuy:    meridional grid distance at u-point
! --- scvx:    zonal grid distance at v-point
! --- scvy:    meridional grid distance at v-point
! --- scpx:    zonal grid distance at p-point
! --- scpy:    meridional grid distance at p-point
! --- scqx:    zonal grid distance at q-point (vorticity)
! --- scqy:    meridional grid distance at q-point
! --- scXXi:   Inverse of scXX
! --- scp2:    scpx*scpy
! --- scp2i:   1/scp2
! --- scu2 :   scux*scuy
! --- scv2 :   scvx*scvy
! --- scq2i:=  1/(scqx*scqy) !Only used once in this form
! ---
! --- --------------------------------------------------------


      !integer, intent(in) :: idm,jdm
      real*8, intent(inout), dimension(0:idm+1,0:jdm+1) ::   
     &   plon ,      ! Longitude of q-points
     &   plat ,      ! Latitude  of q-points
     &   qlon ,      ! Longitude of q-points
     &   qlat ,      ! Latitude  of q-points
     &   ulon ,      ! Longitude of u-points
     &   ulat ,      ! Latitude  of u-points
     &   vlon ,      ! Longitude of v-points
     &   vlat        ! Latitude  of v-points

      real, intent(inout), dimension(idm,jdm) ::   
     & depths         ! water depth

      real, dimension (idm,jdm) ::
     & util1,util2,   ! arrays for temporary storage
     & util3,util4,   ! arrays for temporary storage
     & util5,util6,   ! arrays for temporary storage
     & scux, scuy,    ! mesh size at u pts in x,y dir.
     & scvx, scvy,    ! mesh size at v pts in x,y dir.
     & scpx, scpy,    ! mesh size at p pts in x,y dir.
     & scqx, scqy,    ! mesh size at q pts in x,y dir.
     & scu2, scv2,    ! grid box area at u,v pts
     & scp2, scq2,    ! grid box area at p,q pts
     & scp2i,scq2i,   ! inverses of scp2,scq2
     & scuxi,scvyi,   ! inverses of scux,scvy
     &   pang ,      ! Longitude of u-points
     & asp        ,   ! v-grid aspect ratios for diffusion
     & pgfx, pgfy,    ! horiz. presssure gradient
     & gradx,grady,   ! horiz. presssure gradient
     & depthu,depthv, ! bottom pres. at u,v points
     & pvtrop,        ! pot.vort. of barotropic flow
     & drag,          ! bottom drag
     & corio, rlat,
     & iofld


      integer, dimension (idm,jdm) :: ip

      integer :: nerr,iter,nneigh ! For sanity check of depths matrix


      real realat, hmin, hmax,maxscux
      real :: dlon,dlat
      integer :: ios
      integer im,jm,iread,jread,i,j
      integer im1,ip1,jm1,jp1
      logical tdiag,ex
      character*80 gridid,cline
      character(len=11) :: tag7
      character*48  flnmdep,flnmgrd,flnmrsi,flnmrso,flnmflx,
     &              flnmarc,flnmovr,flnmfor,flnmforw,tmpchar


      ! Names of depths files
      character(len=25) :: bathyfile, newbathyfile

      real :: pi,radian
      integer :: mapflg

      call zaiost()



      tdiag=.TRUE.     !Dump in Tecplot format
      flnmdep='regional.depth'
      flnmgrd='regional.grid'
      pi=2*acos(0.)
      radian=180./pi


! --- ------------------------------------------------------------------
! --- Calc scuy and scvx from qlat, qlon
! ---
! --- stencil:      u--p
! ---               |  |  
! --- We have q   ->q--v
! --- ------------------------------------------------------------------
c$OMP  PARALLEL DO PRIVATE(j,i,im1,jm1,ip1,jp1,realat,dlon,dlat)
c$OMP& SCHEDULE(STATIC,jblk)
      do j=1,jdm
      jm1=j-1
      jp1=j+1
      do i=1,idm
         im1=i-1
         ip1=i+1

         ! Lon, lat is qlon, qlat
         scuy(i,j)=spherdist8(qlon(i,jp1),qlat(i,jp1),   
     &                       qlon(i,j)  ,qlat(i,j))
         scvx(i,j)=spherdist8(qlon(ip1,j),qlat(ip1,j),  
     &                       qlon(i,j)  ,qlat(i,j))
      end do
      end do



! --- ------------------------------------------------------------------
! --- Calc scvy and scux from plat, plon
! ---
! --- stencil:      u--p<- We have p-point
! ---               |  |  
! ---               q--v
! --- ------------------------------------------------------------------
c$OMP  PARALLEL DO PRIVATE(j,i,im1,jm1,ip1,jp1,realat,dlon,dlat)
c$OMP& SCHEDULE(STATIC,jblk)
      do j=1,jdm
      jm1=j-1
      jp1=j+1
      do i=1,idm
         im1=i-1
         ip1=i+1

         ! Lon, lat is plon,plat
         scux(i,j)=spherdist8(plon(i,j)  ,plat(i,j),     
     &                        plon(im1,j),plat(im1,j))
         scvy(i,j)=spherdist8(plon(i,j)  ,plat(i,j),    
     &                        plon(i,jm1),plat(i,jm1))

         ! Rotation angle of x-grid direction rel to local lat line
         dlon=plon(ip1,j)-plon(i,j)
         dlat=plat(ip1,j)-plat(i,j)
         do while (dlon>180. .or.dlon < -180.)
            if (dlon >  180.) dlon = dlon -360.
            if (dlon < -180.) dlon = dlon +360.
         end do
         pang(i,j)=atan2(dlat,dlon)*radian

      end do
      end do


! --- ------------------------------------------------------------------
! --- Calc scpx and scqy from ulat, ulon
! ---
! --- stencil:    ->u--p   We have u-point
! ---               |  |  
! ---               q--v
! --- ------------------------------------------------------------------
c$OMP  PARALLEL DO PRIVATE(j,i,im1,jm1,ip1,jp1,realat,dlon,dlat)
c$OMP& SCHEDULE(STATIC,jblk)
      do j=1,jdm
      jm1=j-1
      jp1=j+1
      do i=1,idm
         im1=i-1
         ip1=i+1

         ! Lon, lat is ulon, ulat
         scpx(i,j)=spherdist8(ulon(ip1,j),ulat(ip1,j),
     &                        ulon(i,j)  ,ulat(i,j))
         scqy(i,j)=spherdist8(ulon(i,j)  ,ulat(i,j),    
     &                        ulon(i,jm1),ulat(i,jm1))
      end do
      end do

! --- ------------------------------------------------------------------
! --- Calc scpy and scqy from vlat,vlon
! ---
! --- stencil:      u--p
! ---               |  |  
! ---               q--v <- We have this point
! --- ------------------------------------------------------------------
c$OMP  PARALLEL DO PRIVATE(j,i,im1,jm1,ip1,jp1,realat,dlon,dlat)
c$OMP& SCHEDULE(STATIC,jblk)
      do j=1,jdm
      jm1=j-1
      jp1=j+1
      do i=1,idm
         im1=i-1
         ip1=i+1

         ! Lon, lat is vlon, vlat
         scqx(i,j)=spherdist8(vlon(i,j)  ,vlat(i,j),    
     &                        vlon(im1,j),vlat(im1,j))
         scpy(i,j)=spherdist8(vlon(i,jp1),vlat(i,jp1),   
     &                        vlon(i,j)  ,vlat(i,j))
      end do
      end do



! --- ------------------------------------------------------------------
! --- define gridsize at u, v, p, and q points
! ---
! --- stencil:    - o    or    u p      
! ---             x |          q v
! ---
! --- ------------------------------------------------------------------
      !print *,mnproc,' -- calc grid props '
      maxscux=0.
c$OMP  PARALLEL DO PRIVATE(j,i,realat)
c$OMP& SCHEDULE(STATIC,jblk)
      do j=1,jdm
      do i=1,idm

         ! Aspect ratio of grid
         asp(i,j)=scpx(i,j)/scpy(i,j)

         ! Inverse grid dimensions
         scuxi(i,j)=1./scux(i,j)
         scvyi(i,j)=1./scvy(i,j)

         ! Grid areas
         scu2(i,j)=scux(i,j)*scuy(i,j)
         scv2(i,j)=scvx(i,j)*scvy(i,j)
         scp2(i,j)=scpx(i,j)*scpy(i,j)
         scq2(i,j)=scqx(i,j)*scqy(i,j)

         ! Inverse grid areas
         scp2i(i,j)=1./scp2(i,j)
         scq2i(i,j)=1./scq2(i,j)

         !latitudes (in rad) for export 
         rlat(i,j)=plat(i,j)/radian               

         ! Coriolis parameter
         realat=qlat(i,j)/radian
         corio(i,j)=sin(realat)*4.*pi/86400.

      ENDDO
      ENDDO
!$OMP END PARALLEL DO 


      ! Map flag:
      ! This is a mess, it says 4= input, but it is really an f-plane.
      ! Lets use mapflg=-1.
      mapflg=-1

      ip=0
      where(depths>.5) ip=1

      ! Dump bathymetry
      open (unit=9,file=flnmdep(1:len_trim(flnmdep))//'.b',
     &   status='unknown',position='rewind')
      write(9,'(a)') 'Depth Grid generated by conformal mapping'
      write(9,'(a)') ' '
      write(9,'(a)') ' '
      write(9,'(a)') ' '
      write(9,'(a)') ' '
      write(9,'(a,2f11.3)') 'min,max depth = ',
     &minval(depths),maxval(depths)
      close(9)
      

      ! Dump header
      write(lp,'(a)') 'Dumping bathymetry To new format'

      ! Dump to grid header file
      open (unit=9,file=flnmgrd(1:len_trim(flnmgrd))//'.b',
     &      status='unknown')
      write(cline,'(i6,a6)') idm,  'idm   ' ;  
      write(9,'(a80)') cline
      write(cline,'(i6,a6)') jdm,  'jdm   ' ;  
      write(9,'(a80)') cline
      write(cline,'(i6,a6)') mapflg,'mapflg' ;  
      write(9,'(a80)') cline
      close(9)
      call zaiopf(flnmdep(1:len_trim(flnmdep))//'.a','replace', 9)
      call zaiowr(depths,ip,.false.,hmin,hmax,9,.true.)
      call zaiocl(9)




      ! Having obtained grids, dump them in a format the "new" hycom will understand
      write(lp,'(a)') 'Converting grid to new' //
     &' hycom format (from version 2.1.03)'

      ! Dump to grid header file
      !open (unit=9,file=flnmgrd(1:len_trim(flnmgrd))//'.b',status='new')
      open (unit=9,file=flnmgrd(1:len_trim(flnmgrd))//'.b',
     &        status='replace')
      write(cline,'(i5,4x,a8)') idm,  "'idm   '" ;  
      write(9,'(a80)') cline
      write(cline,'(i5,4x,a8)') jdm,  "'jdm   '" ;  
      write(9,'(a80)') cline
      write(cline,'(i5,4x,a8)') mapflg,  "'mapflg'" ;  
      write(9,'(a80)') cline
      call zaiopf(flnmgrd(1:len_trim(flnmgrd))//'.a','replace', 9)


      iofld=plon(1:idm,1:jdm)
      !call zaiowr(plon(1:idm,1:jdm),ip,.false.,hmin,hmax,9,.true.)
      call zaiowr(iofld,ip,.false.,hmin,hmax,9,.true.)
      write(9,'(a4,":  min,max =",2f16.5)') "plon",hmin,hmax

      iofld=plat(1:idm,1:jdm)
      !call zaiowr(plat(1:idm,1:jdm),ip,.false.,hmin,hmax,9,.true.)
      call zaiowr(iofld,ip,.false.,hmin,hmax,9,.true.)
      write(9,'(a4,":  min,max =",2f16.5)') "plat",hmin,hmax

      iofld=qlon(1:idm,1:jdm)
      !call zaiowr(qlon(1:idm,1:jdm),ip,.false.,hmin,hmax,9,.true.)
      call zaiowr(iofld,ip,.false.,hmin,hmax,9,.true.)
      write(9,'(a4,":  min,max =",2f16.5)') "qlon",hmin,hmax

      iofld=qlat(1:idm,1:jdm)
      !call zaiowr(qlat(1:idm,1:jdm),ip,.false.,hmin,hmax,9,.true.)
      call zaiowr(iofld,ip,.false.,hmin,hmax,9,.true.)
      write(9,'(a4,":  min,max =",2f16.5)') "qlat",hmin,hmax

      iofld=ulon(1:idm,1:jdm)
      !call zaiowr(ulon(1:idm,1:jdm),ip,.false.,hmin,hmax,9,.true.)
      call zaiowr(iofld,ip,.false.,hmin,hmax,9,.true.)
      write(9,'(a4,":  min,max =",2f16.5)') "ulon",hmin,hmax

      iofld=ulat(1:idm,1:jdm)
      !call zaiowr(ulat(1:idm,1:jdm),ip,.false.,hmin,hmax,9,.true.)
      call zaiowr(iofld,ip,.false.,hmin,hmax,9,.true.)
      write(9,'(a4,":  min,max =",2f16.5)') "ulat",hmin,hmax

      iofld=vlon(1:idm,1:jdm)
      !call zaiowr(vlon(1:idm,1:jdm),ip,.false.,hmin,hmax,9,.true.)
      call zaiowr(iofld,ip,.false.,hmin,hmax,9,.true.)
      write(9,'(a4,":  min,max =",2f16.5)') "vlon",hmin,hmax

      iofld=vlat(1:idm,1:jdm)
      !call zaiowr(vlat(1:idm,1:jdm),ip,.false.,hmin,hmax,9,.true.)
      call zaiowr(iofld,ip,.false.,hmin,hmax,9,.true.)
      write(9,'(a4,":  min,max =",2f16.5)') "vlat",hmin,hmax

      call zaiowr(pang,ip,.false.,hmin,hmax,9,.true.)
      write(9,'(a4,":  min,max =",2f16.5)') "pang",hmin,hmax

      call zaiowr(scpx,ip,.false.,hmin,hmax,9,.true.)
      write(9,'(a4,":  min,max =",2f16.5)') "scpx",hmin,hmax

      call zaiowr(scpy,ip,.false.,hmin,hmax,9,.true.)
      write(9,'(a4,":  min,max =",2f16.5)') "scpy",hmin,hmax

      call zaiowr(scqx,ip,.false.,hmin,hmax,9,.true.)
      write(9,'(a4,":  min,max =",2f16.5)') "scqx",hmin,hmax

      call zaiowr(scqy,ip,.false.,hmin,hmax,9,.true.)
      write(9,'(a4,":  min,max =",2f16.5)') "scqy",hmin,hmax

      call zaiowr(scux,ip,.false.,hmin,hmax,9,.true.)
      write(9,'(a4,":  min,max =",2f16.5)') "scux",hmin,hmax

      call zaiowr(scuy,ip,.false.,hmin,hmax,9,.true.)
      write(9,'(a4,":  min,max =",2f16.5)') "scuy",hmin,hmax

      call zaiowr(scvx,ip,.false.,hmin,hmax,9,.true.)
      write(9,'(a4,":  min,max =",2f16.5)') "scvx",hmin,hmax

      call zaiowr(scvy,ip,.false.,hmin,hmax,9,.true.)
      write(9,'(a4,":  min,max =",2f16.5)') "scvy",hmin,hmax

      call zaiowr(corio,ip,.false.,hmin,hmax,9,.true.)
      write(9,'(a4,":  min,max =",2f16.10)') "cori",hmin,hmax

      call zaiowr(asp,ip,.false.,hmin,hmax,9,.true.)
      write(9,'(a4,":  min,max =",2f16.5)') "pasp",hmin,hmax

      call zaiocl(9)
      close(9)



      ! print a latlon.dat file which can be read from the 
      ! newf program that interpolates forcing functions.
      ! The point is that the program uses -180:180 in longitudes
      ! while micom uses 0:360

      !  For tecplot 
      IF (tdiag) THEN
         open(10,file='teclatlon.dat',form='formatted',status='unknown')
         WRITE(10,*)'TITLE=""'
         WRITE(10,*)'VARIABLES=i,j,modlon,modlat,ulon,ulat,vlon,vlat,'//
     &      'qlon,qlat,depths,corio,scux,scuy,scvx,scvy,scpx,scpy,'//
     &      'scqx,scqy,scp2,scq2,scp2i,scq2i'
         WRITE(10,*)'ZONE I=',idm,',J=',jdm,',F=BLOCK'
         WRITE(10,101)((i,i=1,idm),j=1,jdm)
         WRITE(10,101)((j,i=1,idm),j=1,jdm)

         WRITE(10,100)((plon(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((plat(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((ulon(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((ulat(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((vlon(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((vlat(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((qlon(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((qlat(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((depths(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((corio(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((scux(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((scuy(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((scvx(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((scvy(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((scpx(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((scpy(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((scqx(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((scqy(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((scp2(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((scq2(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((scp2i(i,j),i=1,idm),j=1,jdm)
         WRITE(10,100)((scq2i(i,j),i=1,idm),j=1,jdm)
         close(10)
      end if

100   FORMAT(10(1x,e12.6)) 
101   FORMAT(30i4)
      end subroutine grid_to_hycom
      end module m_grid_to_hycom

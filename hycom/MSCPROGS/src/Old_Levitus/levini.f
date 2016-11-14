      PROGRAM levini
      use mod_xc
      use mod_za
      use mod_grid
      IMPLICIT NONE
c -------------------------------------------------------------
c --- Author Knud Simonsen, NERSC, 10/12 1993
c --- Modified:
c ---   Knut Liseter, model independent version (allocatable vars)
c ---                 Not really used anymore, but kept just
c ---                 in case ...
c -------------------------------------------------------------
c --- The purpose oif this program is to read the Levitus
c --- annual mean temperature and salinity, interpolate them
c --- in a model grid and calculate the depth of the bottom
c --- of the isopycnal layers given by the sigma values in
c --- the file 'gridspec.inp'.
c -------------------------------------------------------------
c --- External subroutines:       External files
c ---      rgrp.f              >>   gridspec.inp
c ---      readsigma.f         >>   gridspec.inp
c ---      newold.f
c ---      rdepth.f            >>   depth.dat  
c ---      levread.f           >>   SALTEMP
c ---      intpol.f
c ---      fillup.f
c ---      eqstat.f
c ---      verpol.f  
c ---      teclev.f
c ---      tec3d.f
c --- --------------------------------------------------------
c --- External files:
c ---      gridspec.inp          Dimension info, sigmavalues
c ---      flags.dat             Landmask
c ---      SALTEMP               Levitus salinity and temp.
c --- --------------------------------------------------------
c --- WARNING:
c ---  Check the dimensions in the routines : fillup.f,rdepth.f
c -------------------------------------------------------------
      INTEGER nxd, nyd, nzd
      PARAMETER(nxd=360,nyd=180,nzd=12) !Data dim. NODC

C     INTEGER iflg(idm,jdm),dpos(10,2),np
C     REAL ypivo,xpivn,ypivn,gridn,  !Grid specifications
C    +     mlat(idm,jdm),              !Model latitudes at pressure point
C    +     mlon(idm,jdm),              !Model longitudes at pressure point
C    +     qlat(0:idm,0:jdm),          !Model latitudes at q- point
C    +     qlon(0:idm,0:jdm),          !Model longitudes at q- point
C    +     ulat(idm,0:jdm),            !Model latitudes at u- point
C    +     ulon(idm,0:jdm),            !Model longitudes at u- point
C    +     vlat(0:idm,jdm),            !Model latitudes at v- point
C    +     vlon(0:idm,jdm),            !Model longitudes at v- point
C    +     mxlon(idm,jdm),mxlat(idm,jdm),!Model (0,1) in geogr. coord. 
C    +     dm(idm,jdm),                !Depth matrix
C    +     mlsalt(idm,jdm,nzd),        !Lev. salt., model grid in hor.  
C    +     mlsig(idm,jdm,nzd),         !Lev. sigm.,model grid in hor.
C    +     salt(nxd,nyd,nzd),        !Levitus saltinity
C    +     zl(nzd),                  !Levitus level depths.  
C    +     dlat(nyd),dlon(nxd)       !Levitus grid positions.
C
C     ! KAL newpos is hardcoded real*8
C     real*8 iolon(idm,jdm),iolat(idm,jdm)

      REAL ypivo,xpivn,ypivn,gridn  !Grid specifications
      integer,allocatable,dimension(:,:) :: 
     +     iflg
      real,allocatable,dimension(:,:) :: 
     +     mlat,                       !Model latitudes at pressure point
     +     mlon,                       !Model longitudes at pressure point
     +     dm,                         !Depth matrix
     +     tmp                         !Depth matrix
      real,allocatable,dimension(:,:,:) :: 
     +     mlsalt,                     !Lev. salt., model grid in hor.  
     +     mlsig,                      !Lev. sigm.,model grid in hor.
     +     salt                        !Levitus saltinity
      real,allocatable,dimension(:) :: 
     +   dlon,
     +   dlat,
     +   zl

      INTEGER :: i,j,k,n
      REAL :: dpos(10,2)


      ! Retrieve model grid and bathymetry
      call xcspmd()
      call zaiost()
      call get_grid()

      ! generic version
      allocate(iflg(idm,jdm))
      allocate(mlat(idm,jdm))
      allocate(mlon(idm,jdm))
      allocate(dm  (idm,jdm))
      allocate(tmp (idm,jdm))
      allocate(mlsalt(idm,jdm,nzd))
      allocate(mlsig(idm,jdm,nzd))
      allocate(salt(nxd,nyd,nzd))
      allocate(dlon(nxd))
      allocate(dlat(nyd))
      allocate(zl(nzd))

      mlat=plat
      mlon=plon
      dm=depths



c
      do 1000 n=1,2
c
      DO i=1,10
       dpos(i,1) = 0
       dpos(i,2) = 0
      ENDDO  

      if (n.eq.1) then
          CALL rNODC(salt,dlon,dlat,zl,nxd,nyd,nzd,1)
      elseif (n.eq.2) then
          CALL rNODC(salt,dlon,dlat,zl,nxd,nyd,nzd,2)
      endif

                                           !Dump data on model
                                           !grid to tecplot file
      CALL tecle2(nxd,nyd,nzd,zl,salt,
     +              'Initial','lev1.dat',
     +              1,nxd,120,nyd,1,nzd)        
                                          
                                           !Load positions from file
C      CALL rlatlon(mlat,mlon)
C     open(10,file='newpos.uf',form='unformatted',status='unknown')
C        !read(10)mlat,mlon
C        read(10)iolat,iolon
C     close(10)
C     mlat=iolat
C     mlon=iolon

                                      !Interpolate horizontally
      CALL intpol(salt,dlat,dlon,nxd,nyd,   !Saln.
     +        mlsalt,mlat,mlon,idm,jdm,nzd)

      WRITE(*,*)'After Intpol salt '

                                      !Load depth matrix 
!      CALL reddat(idm,jdm,1,dm,'Data/depth.forc') 
!      WRITE(*,*)'After reddat'

       DO j=1,jdm                       !Create land-sea mask: 
        DO i=1,idm
!         iflg(i,j)= INT((.5+SIGN(.5,dm(i,j)-1.) ))
          iflg(i,j)= 1   !fill entire model domain
        ENDDO
       ENDDO

      DO k=1,nzd
                                      !Fillup missing values
                                      !in all grid defined as
                                      !seapoints (iflg), also
                                      !if they are below the
                                      !bottom depth.
       CALL fillup(mlsalt,idm,jdm,iflg,nzd,k,10.)
       WRITE(*,*)'After fillup',k
      ENDDO

      if (n.eq.1) then
      OPEN(11,FILE='Data/sss_nodc.dat',form='formatted')
      elseif (n.eq.2) then
      close(11)
      OPEN(11,FILE='Data/sst_nodc.dat',form='formatted')
      endif
c
      DO k=1,nzd                      !Fillup remaining areas
       DO i=1,idm
        DO j=1,jdm 
CKAL     mlati,j)=mlsalt(i,j,k) ! Nice "gotcha"
         tmp(i,j)=mlsalt(i,j,k)
        ENDDO
       ENDDO
CKAL   WRITE(11,'(10f9.4)')mlat
       WRITE(11,'(10f9.4)')tmp
      ENDDO
      CLOSE(11)
         
                                      !Dump to tecplot file
      if (n.eq.1) then
      CALL tecle2(idm,jdm,nzd,zl,mlsalt,
     +              'Initial','Data/tec_nodc_sss.dat',
     +              1,idm,1,jdm,1,nzd)        
      elseif (n.eq.2) then
      CALL tecle2(idm,jdm,nzd,zl,mlsalt,
     +              'Initial','Data/tec_nodc_sst.dat',
     +              1,idm,1,jdm,1,nzd)        
      endif
c
      print *
 1000 continue
c
      STOP 'Finito'
      END  


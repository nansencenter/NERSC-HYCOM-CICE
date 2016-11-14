      program newforc
      use mod_xc
      use mod_za
      use mod_grid , only : depths, mlon => plon, mlat => plat, get_grid
      IMPLICIT NONE
c
c --- Author:Knud Simonsen, NERSC.
c --- First version: 1/12 1993
c --- Revision:     13/1 1994 Fixed some bugs.
c --- Revision:      1/4 1996 Made it general.........
c ---
c ----------------------------------------------------------------
c --- This program creates new model forcing fields 
c ---
c -----------------------------------------------------------------
c --- External datafiles:
c ---         flags.dat: The land/sea mask.
c ---                    See rmask.f
c ---           fort.62: Contains positions of the model grid.  
c ---                    See rlatlong.f
c ---       tmpmean.dat: Contains coastal st. monthly mean  temp. data. 
c ---                    See rstat.f
c ---   Other datafiles: See file fields.inp
c ---
c --- Datafiles.
c ------------------------------------------------------------------
                                           !Size of the model grid

      REAL ypivo,xpivn,ypivn,gridn,       !Model grid specifications
     +     a,b                            !Dummy constants
      INTEGER, PARAMETER:: nrmonths=12    !Nr months (used to be nz?)
      INTEGER, PARAMETER:: IWMAX   =18    !Input fields max len
      INTEGER iw(iwmax)                   !Input parameters
CKAL  INCLUDE 'dimension.h'
CKAL  INTEGER iflg(idm,jdm)               !Land mask
CKAL  REAL       mlat(idm,jdm),           !Model latitudes at pressure point
CKAL +           mlon(idm,jdm),           !Model longitudes at pressure point
CKAL +           mxlat(idm,jdm),          !Unit vector (0,1) in 
CKAL +           mxlon(idm,jdm),          !geogr. components  
CKAL +           mdata(idm,jdm,nrmonths),       !Array for model grid data
CKAL +           cloud(idm,jdm,nrmonths),
CKAL +           humid(idm,jdm,nrmonths),
CKAL +           sst(idm,jdm,nrmonths),
CKAL +           airt(idm,jdm,nrmonths),
CKAL +           wind(idm,jdm,nrmonths),
CKAL +           windev(idm,jdm,nrmonths),
CKAL +           uw(idm,jdm,nrmonths),
CKAL +           vw(idm,jdm,nrmonths),
CKAL +           prcp(idm,jdm,nrmonths),
CKAL +           ice(idm,jdm,nrmonths),
CKAL +           depth(idm,jdm),          !Depth matrixe
CKAL +           dd(0:idm+1,0:jdm+1),
CKAL +           qlat(0:idm,0:jdm)
      INTEGER, ALLOCATABLE:: iflg(:,:)               !Land mask
      REAL, ALLOCATABLE, DIMENSION(:,:) ::
     +           mxlat,          !Unit vector (0,1) in 
     +           mxlon           !geogr. components  

      REAL, ALLOCATABLE, DIMENSION(:,:,:) ::
     +           mdata,       !Array for model grid data
     +           cloud,
     +           humid,
     +           sst,
     +           airt,
     +           wind,
     +           windev,
     +           uw,
     +           vw,
     +           prcp,
     +           ice
      
      

      INTEGER I,J,K
      character*20 filen
  
c --------------------user-supplied statements----------------------
      real radian, pi
      data radian/57.29578/,pi/3.14159265/
c      gridn = 1.                      !Gridresolution
c      ypivo = -10.!-40.               !Model equator meridiane
c      xpivn = 10.!  39.               !i-index for the ME
c      ypivn = -50.!  8.               !j-index for true eq.
                                      !Make your choice such you
                                      !do not hit the pole. 

c ----------------- Init Grid sizes-------------------------------
      call xcspmd()
      call zaiost()
      call get_grid()
      allocate(mxlat (idm,jdm))
      allocate(mxlon (idm,jdm))
      allocate(iflg  (idm,jdm))
      allocate(mdata (idm,jdm,nrmonths))
      allocate(cloud (idm,jdm,nrmonths))
      allocate(humid (idm,jdm,nrmonths))
      allocate(sst   (idm,jdm,nrmonths))
      allocate(airt  (idm,jdm,nrmonths))
      allocate(wind  (idm,jdm,nrmonths))
      allocate(windev(idm,jdm,nrmonths))
      allocate(uw    (idm,jdm,nrmonths))
      allocate(vw    (idm,jdm,nrmonths))
      allocate(prcp  (idm,jdm,nrmonths))
      allocate(ice   (idm,jdm,nrmonths))
      


c ----------------- Init Grid sizes-------------------------------
C .. whats .. the .. point
C     j=1
C     k=1                              !Initialize fields.
C     DO i=1,idm*jdm*nrmonths
C      mdata(i,j,k) = 0.
C      cloud(i,j,k) = 0.8
C      humid(i,j,k) = 0.8
C      sst(i,j,k) = 0.
C      airt(i,j,k) = 0.
C      wind(i,j,k) = 0.
C      windev(i,j,k) = 0.
C      prcp(i,j,k) = 0.
C      uw(i,j,k) = 0.
C      vw(i,j,k) = 0. 
C     ENDDO

      j=1
      k=1                              !Initialize fields.
      DO k=1,nrmonths
      DO j=1,jdm
      DO i=1,idm
       mdata(i,j,k) = 0.
       cloud(i,j,k) = 0.8
       humid(i,j,k) = 0.8
       sst(i,j,k) = 0.
       airt(i,j,k) = 0.
       wind(i,j,k) = 0.
       windev(i,j,k) = 0.
       prcp(i,j,k) = 0.
       uw(i,j,k) = 0.
       vw(i,j,k) = 0. 
      ENDDO
      ENDDO
      ENDDO
c --- ------------------------------- Load commands from fields.inp ---

 200  OPEN(10,FILE='fields.inp',STATUS='old')
        READ(10,'(18I3)')(iw(i),i=1,iwmax)
      CLOSE(10)

c --- ------------------------------- Grid positions ---------------------
! Load pos from file.
!      CALL rlatlon(mlat,mlon)
CKAL  open(10,file='newpos.uf',form='unformatted',status='unknown')
CKAL     !read(10)mlat,mlon
CKAL     read(10)iolat,iolon
CKAL  close(10)
CKAL  mlat=iolat
CKAL  mlon=iolon


c --- ------------------------------- SST field  ----
      IF(iw(5).GT.0)THEN
      CALL csst(sst,iw(5),idm,jdm,12,iflg,mlon,mlat)
      CALL dmpdat(idm,jdm,nrmonths,sst,10000,
     +                                 'Data/sstmp.forc')  
      ENDIF

c --- ------------------------------- Set iflg=1 everywhere  ----
      DO i=1,idm
       DO j=1,jdm
        iflg(i,j)=1
       ENDDO
      ENDDO
c --- ------------------------------- Wind velocity field  ----
      IF(iw(1).GT.0)THEN
      CALL cuv(uw,vw,iw(1),idm,jdm,nrmonths,         
     +        iflg,mlon,mlat,mxlon,mxlat)
      CALL dmpdat(idm,jdm,nrmonths,uw,10000,'Data/uwind.forc')  
      CALL dmpdat(idm,jdm,nrmonths,vw,10000,'Data/vwind.forc')  
      ENDIF

c --- ------------------------------- St Wind deviation  ----
      IF(iw(2).GT.0)THEN
      CALL cwstde(windev,iw(2),idm,jdm,nrmonths,         
     +        iflg,mlon,mlat)
      CALL dmpdat(idm,jdm,nrmonths,windev,10000,
     +                                'Data/wnddv.forc')  
      ENDIF

c --- ------------------------------- Abs. Wind field  ----
      IF(iw(3).GT.0)THEN
      CALL cwind(wind,windev,uw,vw,iw(3),idm,jdm,12,         
     +        iflg,mlon,mlat)
      CALL dmpdat(idm,jdm,nrmonths,wind,10000,
     +                                'Data/wndab.forc')  
      ENDIF

c --- ------------------------------- Air temp field  ----
      IF(iw(4).GT.0)THEN
      CALL cairt(airt,sst,iw(4),idm,jdm,12,iflg,mlon,mlat)
      CALL dmpdat(idm,jdm,nrmonths,airt,10000,
     +                                 'Data/airtp.forc')  
                             !Note: the sst array is used as
                             !dummy array in this routine
      ENDIF
 
c --- ------------------------------- Cloud field  ----
      IF(iw(6).GT.0)THEN
      CALL ccloud(cloud,iw(6),idm,jdm,12,iflg,mlon,mlat)
      CALL dmpdat(idm,jdm,nrmonths,cloud,10000,
     +                                 'Data/cloud.forc')  
      ENDIF

c --- ------------------------------- Relative humidity  ----
      IF(iw(7).GT.0)THEN
      CALL chumid(humid,iw(7),idm,jdm,12,iflg,mlon,mlat)
      CALL dmpdat(idm,jdm,nrmonths,humid,10000,
     +                                  'Data/humid.forc')  
      ENDIF
c --- ------------------------------- Precipitation  ----
      IF(iw(8).GT.0)THEN
      CALL crain(prcp,iw(8),idm,jdm,12,iflg,mlon,mlat)
      CALL dmpdat(idm,jdm,nrmonths,prcp,1000,
     +                              'Data/precp.forc')  
      ENDIF
c --- ------------------------------- Sea ice concentration ----
      IF(iw(11).GT.0)THEN
      CALL cice(ice,iw(11),idm,jdm,12,iflg,mlon,mlat)
      CALL dmpdat(idm,jdm,nrmonths,ice,10,'Data/icec.forc')
      ENDIF

c --- ------------------------------- Dump to tecplot files
c      filen='c.dat'
c      CALL tec1(idm,jdm,12,filen,cloud)
       filen = 'Data/forc.dat'
      CALL tecall(idm,jdm,12,filen,iflg,
     +            mlon,mlat,wind,ice,prcp,airt,cloud,humid,uw,vw)
  99  FORMAT(20I4)
 100  FORMAT(10(1x,e12.6))

 900  STOP      
      END 


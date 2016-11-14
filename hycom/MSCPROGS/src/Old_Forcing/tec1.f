      SUBROUTINE tec1(nx,ny,mo,filen,mlon,mlat,cc,mean)
c ---------------------------------------------------------------------
c --- Author: Knud Simonsen, NERSC 6/12 1993
c ---------------------------------------------------------------------
c --- Dump the vector desriped by (u,v) on file in a format accepted
c --- by tecplot.
c --- Input:
c ---        cc : cloud cover
c ---------------------------------------------------------------------
      CHARACTER*10 filen
      CHARACTER*7 month(1:12)
      INTEGER nx,ny,mo
      REAL cc(nx,ny,mo),
     &    mlon(nx,ny),mlat(nx,ny),mean
   
      month(1)='Jan. '
      month(2)='Feb. '
      month(3)='March'
      month(4)='April'
      month(5)='Mai  '
      month(6)='June '
      month(7)='July '
      month(8)='Aug. '
      month(9)='Sept.'
      month(10)='Oct.'
      month(11)='Nov.'
      month(12)='Dec.'

      OPEN(10,FILE=filen,STATUS='UNKNOWN')
      WRITE(10,101)filen
      WRITE(10,102)
      DO k=1,mo
       WRITE(10,103)month(k),nx,ny
        WRITE(10,100)((mlon(i,j),i=1,nx),j=1,ny)
        WRITE(10,100)((mlat(i,j),i=1,nx),j=1,ny)
        WRITE(10,100)(((cc(i,j,k)-mean),i=1,nx),j=1,ny)
      ENDDO
      CLOSE(10)
      WRITE(*,*)'Field written to ',filen
   99 FORMAT(20I4)  
  100 FORMAT(10(1x,e10.3))
  101 FORMAT('TITLE="',A10,'"')
  102 FORMAT('VARIABLES="lon","lat","CC"')
  103 FORMAT('ZONE T="',A7,'", I=',I3,', J=',I3,', F=BLOCK')
      RETURN
      END    




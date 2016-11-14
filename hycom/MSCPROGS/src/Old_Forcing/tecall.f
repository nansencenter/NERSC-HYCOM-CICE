      SUBROUTINE tecall(nx,ny,mo,filen,iflg,
     +           lon,lat,w,ice,prcp,ta,cc,rh,u,v)
c ---------------------------------------------------------------------
c --- Author: Knud Simonsen, NERSC 6/12 1993
c ---------------------------------------------------------------------
c --- Dump the vector desriped by (u,v) on file in a format accepted
c --- by tecplot.
c --- Input:
c ---        cc : cloud cover
c ---------------------------------------------------------------------
      CHARACTER*20 filen
      CHARACTER*7 month(1:12)
      INTEGER nx,ny,mo,iflg(nx,ny)
      REAL cc(nx,ny,mo),w(nx,ny,mo),
     +     rh(nx,ny,mo),prcp(nx,ny,mo),
     +     ta(nx,ny,mo),u(nx,ny,mo),
     +     ice(nx,ny,mo),
     +     v(nx,ny,mo),ti,
     +     lat(nx,ny),lon(nx,ny)
   
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
      ti = 273.15

      OPEN(10,FILE=trim(filen),STATUS='UNKNOWN')
      WRITE(*,*)'Open ',filen
      WRITE(10,101)filen
      WRITE(10,102)
      DO k=1,mo
       WRITE(10,103)month(k),nx,ny
        if (k==1) then
           WRITE(10,100)((lon(i,j),i=1,nx),j=1,ny)
           WRITE(10,100)((lat(i,j),i=1,nx),j=1,ny)
        else
           write(10,'(a)') 'D=(1,2)'
        end if
        WRITE(10,100)((w(i,j,k),i=1,nx),j=1,ny)
        WRITE(10,100)((ice(i,j,k),i=1,nx),j=1,ny)
        WRITE(10,100)(((prcp(i,j,k)),i=1,nx),j=1,ny)
        WRITE(10,100)(((ta(i,j,k)-ti),i=1,nx),j=1,ny)
        WRITE(10,100)((cc(i,j,k),i=1,nx),j=1,ny)
        WRITE(10,100)((rh(i,j,k),i=1,nx),j=1,ny)
        WRITE(10,100)((u(i,j,k),i=1,nx),j=1,ny)
        WRITE(10,100)((v(i,j,k),i=1,nx),j=1,ny)
        WRITE(10,100)((FLOAT(iflg(i,j)),i=1,nx),j=1,ny)
      ENDDO
      CLOSE(10)
   99 FORMAT(20I4)  
  100 FORMAT(10(1x,e11.4))
  101 FORMAT('TITLE= "',A10,'"')
  102 FORMAT('VARIABLES="lon","lat","W","Ice",',
     +        '"P","Ta","CC","RH","U","V","fl"')
C 102 FORMAT('VARIABLES="lon","lat","Ta"')
  103 FORMAT('ZONE T="',A7,'", I=',I3,', J=',I3,', F=BLOCK')
      RETURN
      END    




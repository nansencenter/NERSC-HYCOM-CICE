      SUBROUTINE rwalsh(dd,dlat,dlon,nxd,nyd,
     +               filen,a,b,vmin,vmax)
c -----------------------------------------------------------------
c --- Author Knud Simonsen, NERSC,25/2 1997.
c --- -------------------------------------------------------------
c --- Reads the Walsh ice concentration data (output from
c --- programme read1.f written by HD)
c -----------------------------------------------------------------
c --- Input: filen: filaname 
c ---        nxd,nyd: size of data file. 
c ---        a,b: const. to convert.
c ---        vmin,vmax: min and max accepted values
c --- Ouput: dd: data array in physical units.
c ---        dlon: contains the longitudes.
c ---        dlat: contains the lattitudes
c ----------------------------------------------------------------
      REAL dd(nxd,nyd,12),
     +     dlat(nxd,nyd),dlon(nxd,nyd)
      CHARACTER filen*20,dum*80
      REAL a,b,vmin,vmax,val,td
      INTEGER x,y,k,mo ,id(180)

      IF(nxd.NE.80)STOP 'Check rwalsh'

      WRITE(*,*)'Read datafile: ',filen
      OPEN(89,FILE=filen,ACCESS='SEQUENTIAL',FORM='FORMATTED',ERR=602) 
 
      DO i=1,nxd
       DO j=1,nyd
        READ(89,'(2f9.3,12f7.3)')
     &    dlat(i,j),dlon(i,j),(dd(i,j,m),m=1,12)
        DO m=1,12
         IF(dd(i,j,m).LT.0.0) dd(i,j,m)=-999.9999
        ENDDO
       ENDDO
      ENDDO 
      CLOSE(89)
      GOTO 603       
  602 STOP 'File Error, rwalsh.f'
  603 RETURN
      END

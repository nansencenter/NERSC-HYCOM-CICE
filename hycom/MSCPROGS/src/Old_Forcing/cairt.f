      SUBROUTINE cairt(airt,th,iw,nx,ny,mo,
     +        iflg,mlon,mlat)
c --- ----------------------------------------------------
c --- Call the routines nessesary to interpolate data from
c --- data file to model grid.
c --- Positions: mlon,mlat
c --- Land mask: iflg
c --- Data source: spec by iw

      INTEGER iw,nx,ny,mo,iflg(nx,ny)
      REAL airt(nx,ny,mo),mlon(nx,ny),mlat(nx,ny),
     +     th(nx,ny,mo),
     +    a,b,vmin,vmax
      CHARACTER*15 filen
     
      IF(iw.EQ.1) THEN 
       filen= 'Data/AIRGLOB'    !Read COADS air temp. In the
       a=100.                      !Arctic the data from Shea(1986)
       b = 273.15-60.              !and the Antarctic the data
       vmin = 211.
       vmax = 321.
       CALL lcoads(airt,iflg,nx,ny,mlon,mlat,filen,a,b,
     +               vmin,vmax)
                                     !Fillup missing data
                                     !and land mask.
       CALL fillup(airt,nx,ny,iflg,12,271.0)
      ELSE IF(iw.EQ.2) THEN
       filen= 'Data/TSC'            !Read ECMWF temp
       a=10000./120.
       b= 273.16-75.
       vmin = 211.
       vmax = 321.
       CALL lecmw2(airt,iflg,nx,ny,mlon,mlat,filen,a,b,
     +            vmin,vmax)
                                    !Fillup missing data
                                     !and land mask.
       CALL fillup(airt,nx,ny,iflg,12,271.0)
      ELSE IF(iw.EQ.3)THEN
                                    !Load/generates air tempdata,
                                    !which is a mixture of COADS
                                    !EXMWF and coastal weather stations  
       CALL mixairt(airt,th,nx,ny,iflg,mlon,mlat)

      ELSE
       WRITE(*,*)'Airt field does not exist'
      ENDIF
c.diag      filen='airt.dat'                    !Dump to tec.
c.diag      CALL tec1(nx,ny,12,filen,mlon,mlat,airt,273.15)
c.diag      stop'test'

      RETURN
      END




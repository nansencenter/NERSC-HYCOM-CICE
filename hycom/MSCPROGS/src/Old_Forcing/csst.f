      SUBROUTINE csst(sst,iw,nx,ny,mo,
     +        iflg,mlon,mlat)
c --- ----------------------------------------------------
c --- Call the routines nessesary to interpolate data from
c --- data file to model grid.
c --- Positions: mlon,mlat
c --- Land mask: iflg
c --- Data source: spec by iw
c --- Fiels: Sea Surface Temperature
c --- ----------------------------------------------------

      INTEGER iw,nx,ny,mo,iflg(nx,ny)
      REAL sst(nx,ny,mo),mlon(nx,ny),mlat(nx,ny),
     +    a,b,vmin,vmax
      CHARACTER*15 filen
     
      IF(iw.EQ.1) THEN 
       filen= 'Data/SSTGLOB'       !Read COADS SST.
       a=100.                      !Where COADS data not defined
       b= 273.15                   !the Reynolds data are applied
       vmin = 271.34
       vmax = 321.
       CALL lcoads(sst,iflg,nx,ny,mlon,mlat,filen,a,b,
     +               vmin,vmax)
                                      !Fillup missing data
                                      !and land mask.
       CALL fillup(sst,nx,ny,iflg,12,273.1)
      ELSE IF(iw.EQ.2) THEN

       filen= 'Data/AMICLIM'            !Read AMIP SST
       a=100.
       b = 273.15
       vmin = 271.34
       vmax =321.
       CALL lamip(sst,iflg,nx,ny,mlon,mlat,filen,a,b,
     +               vmin,vmax)
                                      !Fillup missing data
                                      !and land mask.
       CALL fillup(sst,nx,ny,iflg,12,273.1)
      ELSE
       WRITE(*,*)'SST field does not exist'
      ENDIF
      RETURN
      END




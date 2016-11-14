      SUBROUTINE cwstde(wind,iw,nx,ny,mo,
     +        iflg,mlon,mlat)

      INTEGER iw,nx,ny,mo,iflg(nx,ny)
      REAL wind(nx,ny,mo),mlon(nx,ny),mlat(nx,ny),
     +    a,b,vmin,vmax
      CHARACTER*15 filen
       
     
      IF(iw.EQ.1) THEN 
       filen= 'Data/UVDEV'          !Read COADS wind speed.
       a= 100.                 !st. dev.
       b= 0.
       vmin = 0.
       vmax = 20.
       CALL lcoads(wind,iflg,nx,ny,mlon,mlat,filen,a,b,
     +               vmin,vmax)
                                      !Fillup missing data
                                      !and land mask.
       CALL fillup(wind,nx,ny,iflg,12,1.0)

      ELSE IF(iw.EQ.2) THEN
       filen= 'Data/STDV_USC'       !Read ECMWF wind speed
       a=500.                  !st. dev.    
       b = 0.
       vmin = 0.
       vmax =20.
       CALL lecmwf(wind,iflg,nx,ny,mlon,mlat,filen,a,b,
     +              vmin,vmax)
                                      !Fillup missing data
                                      !and land mask.
       CALL fillup(wind,nx,ny,iflg,12,1.0)
      ELSE
       WRITE(*,*)'Wind speed st. dev. is not loaded'
      ENDIF
      RETURN
      END




      SUBROUTINE crain(prcp,iw,nx,ny,mo,
     +        iflg,mlon,mlat)

      INTEGER iw,nx,ny,mo,iflg(nx,ny)
      REAL prcp(nx,ny,mo),mlon(nx,ny),mlat(nx,ny),
     +    a,b,vmin,vmax
      CHARACTER*15 filen
       
     
      IF(iw.EQ.1) THEN 
       filen= 'Data/LEGATES'         !Read Relative Humidity
c       a= 1000.*2592000.               !m/s
       WRITE(*,*)'Precipitation is given in mm/month (crain.f)'
       a=1.                             !mm/month
       b= 0.
       vmin=-0.1
       vmax=1000.0
       CALL lrain(prcp,iflg,nx,ny,mlon,mlat,filen,a,b,
     +            vmin,vmax)
                                      !Fillup missing data
                                      !and land mask.
       CALL fillup(prcp,nx,ny,iflg,12,100.)


   
      ELSE
       WRITE(*,*)'Rel.Humid. does not exist'
      ENDIF
      RETURN
      END




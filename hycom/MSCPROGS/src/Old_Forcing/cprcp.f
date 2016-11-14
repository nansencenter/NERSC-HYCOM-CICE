      SUBROUTINE cprcp(prcp,iw,nx,ny,mo,
     +        iflg,mlon,mlat)

      INTEGER iw,nx,ny,mo,iflg(nx,ny)
      REAL humid(nx,ny,mo),mlon(nx,ny),mlat(nx,ny),
     +    a,b,vmin,vmax,ncc(12)
      CHARACTER*15 filen
       
     
      IF(iw.EQ.1) THEN 
       filen= 'Data/WETNESS'         !Read Relative Humidity
       a= 1000.
       b= 0.
       vmin=0.1
       vmax=1.0
       CALL lcoads(prcp,iflg,nx,ny,mlon,mlat,filen,a,b,
     +            vmin,vmax)
                                      !Fillup missing data
                                      !and land mask.
       CALL fillup(prcp,nx,ny,iflg,12,0.8)

      ELSE
       WRITE(*,*)'Rel.Humid. does not exist'
      ENDIF
      RETURN
      END




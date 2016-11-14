      SUBROUTINE ccloud(cloud,iw,nx,ny,mo,
     +        iflg,mlon,mlat)

      INTEGER iw,nx,ny,mo,iflg(nx,ny)
      REAL cloud(nx,ny,mo),mlon(nx,ny),mlat(nx,ny),
     +    a,b,vmin,vmax
      CHARACTER*15 filen
       
     
      IF(iw.EQ.1.OR.iw.EQ.2) THEN 
       filen= 'Data/COVER'            !Read Cloud cover
       a= 80.
       b= 0.
       vmin=0.1
       vmax=1.0
       CALL lcoads(cloud,iflg,nx,ny,mlon,mlat,filen,a,b,
     +            vmin,vmax)
      ENDIF
      IF (iw.EQ.2) THEN               !Include some values in no 
       CALL inclcc(cloud,nx,ny,       !data point in the Arctic.
     +      mlon,mlat)
      ENDIF
                                      !Fillup missing data
                                      !and land mask. 

       CALL fillup(cloud,nx,ny,iflg,12,0.7)

   
      IF(iw.LT.1.OR.iw.GT.2) THEN 
       WRITE(*,*)'Cloud Cover does not exist'
      ENDIF
      RETURN
      END




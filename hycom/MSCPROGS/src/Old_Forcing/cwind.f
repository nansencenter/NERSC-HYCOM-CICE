      SUBROUTINE cwind(wind,stw,uw,vw,iw,nx,ny,mo,
     +        iflg,mlon,mlat)

      INTEGER iw,nx,ny,mo,iflg(nx,ny)
      REAL wind(nx,ny,mo),mlon(nx,ny),mlat(nx,ny),
     +    stw(nx,ny,mo),uw(nx,ny,mo),vw(nx,ny,mo),
     +    a,b,vmin,vmax,uvabs,rflg,fl
      CHARACTER*15 filen
       
      write(*,*)'Choiche of abs. wind: ',iw 
      IF(iw.EQ.1.OR.iw.EQ.3.OR.iw.EQ.5)THEN 
       filen= 'Data/UVABS'       !Read COADS wind speed.
       a= 10.
       b= 0.
       vmin = 0.
       vmax = 20.
       CALL lcoads(wind,iflg,nx,ny,mlon,mlat,filen,a,b,
     +               vmin,vmax)
       
      ENDIF

      IF(iw.EQ.2.OR.iw.EQ.4.OR.iw.EQ.6)THEN
       filen= 'Data/USC'            !Read ECMWF wind
       a=500.                      
       b = 0.
       vmin = 0.
       vmax =20.
       CALL lecmwf(wind,iflg,nx,ny,mlon,mlat,filen,a,b,
     +              vmin,vmax)
      ENDIF 

      IF(iw.EQ.3.OR.iw.EQ.4) THEN
       k=1                     !Add st. dev to the wind field
       j=1
       DO i=1,nx*ny*mo
        wind(i,j,k)=wind(i,j,k)+MAX(0.,stw(i,j,k))
       ENDDO
      ENDIF
 
      IF(iw.EQ.5.OR.iw.EQ.6) THEN 
                                   !Fillup with sqrt(u2+v2)
       DO j=1,ny                   !if no data exist 
        DO i=1,nx
         rflg=FLOAT(iflg(i,j))     !Land/sea flag 
         DO k=1,mo 
                                   !abs. wind
                                   ! -999 if land  
          uvabs = uw(i,j,k)*uw(i,j,k)+vw(i,j,k)*vw(i,j,k)
          uvabs = SQRT(uvabs)*rflg - (1.-rflg)*999.
                                   !=1 data exist
                                   !=0 no data or land 
          fl=(0.5 + SIGN(0.5,wind(i,j,k)-.01) )*rflg

          wind(i,j,k)=wind(i,j,k)*fl+(1.-fl)*uvabs
         ENDDO       
        ENDDO
       ENDDO
      ENDIF
                                   !Fillup missing data
                                   !and land mask.
      CALL fillup(wind,nx,ny,iflg,12,5.)

      IF(iw.GT.6.OR.iw.LT.1) THEN
       WRITE(*,*)'Abs. wind field does not exist'
      ENDIF
      RETURN
      END




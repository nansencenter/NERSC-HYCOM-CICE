      SUBROUTINE cuv(uw,vw,iw,nx,ny,mo,
     +        iflg,mlon,mlat,mxlon,mxlat)

      INTEGER iw,nx,ny,mo,iflg(nx,ny)
      REAL uw(nx,ny,mo),vw(nx,ny,mo),mlon(nx,ny),mlat(nx,ny),
     +    mxlon(nx,ny,mo),mxlat(nx,ny,mo),
     +    a,b,vmin,vmax
      CHARACTER*15 filen
       
     
      IF(iw.EQ.1) THEN 
                                        !Read He. and Ro wind field
       filen = 'Data/UVGLOB'
       CALL lherov(uw,vw,mxlat,mxlon,mlat,mlon,
     +    iflg,nx,ny,gridn,0.0)

      ELSE IF(iw.EQ.2) THEN
       filen= 'Data/TAUXY'            !Read ECMWF wind stress
       a= 5000.
       b= -1.
       vmin = -30.
       vmax = +30.
       CALL lecmwv(uw,vw,iflg,nx,ny,mlon,mlat,
     +               mxlon,mxlat,filen,a,b,
     +               vmin,vmax,0.0)

      ELSE
       WRITE(*,*)'Wind velocity field does not exist'
      ENDIF
      RETURN
      END




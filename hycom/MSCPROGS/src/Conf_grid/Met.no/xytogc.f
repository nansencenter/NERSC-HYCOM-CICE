C
C ********************************************************************
C
      SUBROUTINE XYTOGC(NPOS,X,Y,GRID)
C
C        From polarstereographic to geographic coordinates.
C        NB: NO OTHER CONVERSION POSSIBLE
C            Use 'xyconvert.f' for more general interpolation 
C
C  Input:
C        NPOS:       no. of positions
C        X(NPOS):    polarstereographic X coordinate
C        Y(NPOS):    polarstereographic Y coordinate
C        GRID(1-6):  XP,YP,AN,FI,projection latitude,not used
C
C  Output:
C        X(NPOS):    geographic latitude  (decimal, +/- = N/S)
C        Y(NPOS):    geographic longitude (decimal, +/- = E/W)
C
C----------------------------------------------------------------------
C  DNMI/FoU   30.11.1992   Anstein Foss
C  DNMI/FoU   14.11.1998   Anstein Foss
C  DNMI/FoU   24.01.2003   Harald Engedahl
C-------------------------------------------------------------------
C
      implicit none
c
      INTEGER NPOS,N
      REAL    X(NPOS),Y(NPOS),GRID(6)
      real    XP,YP,AN,FI,DEG,DEG2,DXP,DYP,RR,GLAT,GLON
C
      XP=GRID(1)
      YP=GRID(2)
      AN=GRID(3)
      FI=GRID(4)
C
      DEG=180./3.141592654
      DEG2=DEG*2.
C
      DO N=1,NPOS
        DXP=X(N)-XP
        DYP=YP-Y(N)
        RR=SQRT(DXP*DXP+DYP*DYP)
        GLAT=90.-DEG2*ATAN(RR/AN)
        GLON=0.
        IF(RR.GT.1.E-20)  GLON=FI+DEG*ATAN2(DXP,DYP)
        IF(GLON.LE.-180.) GLON=GLON+360.
        IF(GLON.GT.+180.) GLON=GLON-360.
        X(N)=GLAT
        Y(N)=GLON
      END DO
C
      RETURN
      END

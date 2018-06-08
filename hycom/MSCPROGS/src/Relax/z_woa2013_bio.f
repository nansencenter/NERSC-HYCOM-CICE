      PROGRAM TSZINT
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     DEFINE INPUT CLIMATOLOGY GRID.
C
C     SETUP FOR 1.0 DEGREE LEVITUS GLOBAL CLIMATOLOGY.
C
      INTEGER    IWI,JWI,KWI
      REAL       XFIN,YFIN,DXIN,DYIN
      PARAMETER (IWI=360, JWI=180, KWI=102)
      PARAMETER (XFIN=-179.5, YFIN=-89.5, DXIN=1.0, DYIN=1.0)
C
C     CLIM ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL,    ALLOCATABLE :: PLON(:,:),PLAT(:,:),XAF(:,:),YAF(:,:)
      REAL,    ALLOCATABLE :: NITM(:,:),PHOM(:,:),SILM(:,:),OXYM(:,:)
C
      REAL      XAMAX,XAMIN,YAMAX,YAMIN
      REAL*4    ZLEV(KWI)
      REAL*4    NITSEAIN(IWI,JWI),PHOSEAIN(IWI,JWI),
     +          SILSEAIN(IWI,JWI),OXYSEAIN(IWI,JWI)
      CHARACTER PREAMBL(5)*79
C
C     INTERPOLATION ARRAYS.
C
      INTEGER IBD(4)
      REAL    NITSEAI(IWI+4,JWI+4),PHOSEAI(IWI+4,JWI+4),
     +        SILSEAI(IWI+4,JWI+4),OXYSEAI(IWI+4,JWI+4)
      REAL    FXI(IWI+4),FYI(JWI+4),
     +        WQSEA3(IWI+4,JWI+4,3),WK(3*(IWI+JWI+8)+1)
C
C     NAMELIST.
C
      INTEGER          JPR
      COMMON/NPROCS/   JPR
      SAVE  /NPROCS/
C
      CHARACTER*40     CTITLE
      NAMELIST/AFTITL/ CTITLE
      INTEGER          ICTYPE,IFRZ,KSIGMA,INTERP,MONTH,ITEST,JTEST
      NAMELIST/AFFLAG/ ICTYPE,IFRZ,KSIGMA,INTERP,MONTH,ITEST,JTEST,JPR
C
C**********
C*
C 1)  FROM UNFORMATTED Sig0-T OR Sig0-S OR T-S DATA ON ITS NATIVE GRID,
C      CREATE A FORMATTED MODEL GRID CLIM FILE SUITABLE FOR INPUT TO THE
C      HYCOM ISOPYCNAL CLIMATOLOGY GENERATOR OVER THE GIVEN REGION.
C
C     INTERPOLATION IS EITHER PIECEWISE BILINEAR OR CUBIC SPLINE.
C
C 2)  PARAMETERS:
C
C     NATIVE CLIM GRID SPECIFICATION, SEE (4):
C
C        IWI    = 1ST DIMENSION OF CLIM GRID
C        JWI    = 2ND DIMENSION OF CLIM GRID
C        JWI    = 3RD DIMENSION OF CLIM GRID (NUMBER OF Z-LEVELS)
C        XFIN   = LONGITUDE OF 1ST CLIM GRID POINT
C        YFIN   = LATITUDE  OF 1ST CLIM GRID POINT
C        DXIN   = CLIM LONGITUDINAL GRID SPACING
C        DYIN   = CLIM LATITUDINAL  GRID SPACING
C
C 3)  NAMELIST INPUT:
C
C     /AFTITL/
C        CTITLE - ONE (40-CHARACTER) LINE TITLE.
C
C     /AFFLAG/
C        ICTYPE -  INPUT FILE TYPE
C                   =1; Sigma-0 AND POTENTIAL TEMPERATURE (DEAFULT)
C                   =2; Sigma-0 AND SALINITY
C                   =3; POTENTIAL TEMPERATURE AND SALINITY
C                   =4; POTENTIAL TEMPERATURE AND SALINITY, UNSTABLE
C        IFRZ   -  FREEZING POINT TYPE
C                   =0; NONE
C                   =1; CONSTANT (DEFAULT)
C                   =2; SALINITY DEPENDENT
C        KSIGMA - OUTPUT FILE TYPE
C                   =0; Sigma-0, POTENTIAL TEMPERATURE AND SALINITY
C                   =2; Sigma-2, POTENTIAL TEMPERATURE AND SALINITY
C        INTERP - INTERPOLATION FLAG.
C                   =0; PIECEWISE LINEAR
C                   =1; CUBIC SPLINE (DEFAULT)
C        MONTH  - MONTH OF CLIMATOLOGY (1 TO 12)
C        ITEST  - 1ST ARRAY INDEX FOR DIAGNOSTIC PRINTOUT
C                   =0; NO DIAGNOSTIC PRINTOUT
C        JTEST  - 2ND ARRAY INDEX FOR DIAGNOSTIC PRINTOUT
C                   =0; NO DIAGNOSTIC PRINTOUT
C
C 4)  INPUT:
C        ON UNIT  5: NAMELIST /AFTITL/, /AFTIME/
C        ON UNIT 71: UNFORMATTED NATIVE Sig0 CLIM FILE, SEE (5).
C        ON UNIT 72: UNFORMATTED NATIVE PotT CLIM FILE, SEE (5).
C        ON UNIT 73: UNFORMATTED NATIVE S    CLIM FILE, SEE (5).
C     OUTPUT:
C        ON UNIT 11:    UNFORMATTED MODEL    CLIM FILE, SEE (6).
C        ON UNIT 12:    UNFORMATTED MODEL    CLIM FILE, SEE (6).
C        ON UNIT 13:    UNFORMATTED MODEL    CLIM FILE, SEE (6).
C
C 5)  THE INPUT CLIM FIELDS ON UNITS 71-73, ARE ON THE 'NATIVE' LAT-LON
C      GRID, STARTING AT THE POINT 'XFIN' EAST AND 'YFIN' NORTH WITH 
C      'YFIN' NORTH WITH GRIDSIZE 'DXIN' BY 'DYIN' DEGREES.  THE
C      INPUT ARRAY SIZE IS 'IWI' BY 'JWI', AND THERE ARE NO REPEATED
C      NODES (EVEN FOR GLOBAL DATA SETS).  THE CONTENTS OF EACH INPUT
C      FILE IS AS FOLLOWS:
C       RECORD 1:   A 40-CHARACTER TITLE
C       RECORD 2:   IWI,JWI,KWI,XFIN,YFIN,DXIN,DYIN,ZLEV
C       RECORD 2+N: CLIM Z-LEVEL N, N=1...KWI.
C
C     ALL CLIMATOLOGY FIELDS MUST BE DEFINED AT EVERY GRID POINT,
C      INCLUDING LAND AND BELOW THE OCEAN FLOOR.
C     THE POTENTIAL DENSITY FIELDS MUST BE STABALY STRATIFIED
C      IN THE VERTICAL EVERYWHERE, AND THEREFORE SHOULD NOT BE
C      INPUT FOR KSIGMA=2 (Sigma2 OUTPUT), I.E. THEN USE ICTYPE=3,4.
C
C     IF ICTYPE=4, DON'T INFORCE STABLE DENSITY STRATIFICATION.
C
C 6)  THE OUTPUT CLIMS ARE AT EVERY GRID POINT OF THE MODEL'S 'P' GRID.
C     ARRAY SIZE IS 'IDM' BY 'JDM'.
C
C 7)  SEVERAL STATISTICS ARE WRITTEN OUT IN ORDER TO CHECK THE 
C      INTERPOLATION BETWEEN VARIOUS MACHINES.  MIN, MAX, MEAN AND RMS 
C      OF THE ENTIRE BASIN ARE OUTPUT FOR (TM, SM, RM) FOR EACH 
C      RECORD.  NOTE HOWEVER THAT THESE VALUES MAY NOT REPRESENT THE
C      STATISTICS OF THE CLIMS AS SEEN BY THE MODEL, IF THE INPUT CLIM 
C      DATA HAS NON-REALISTIC VALUES OVER LAND.  IT IS UP TO THE USER 
C      TO CHECK THE LOG FILES FOR CONSISTENCY BETWEEN MACHINES.
C
C 8)  ALAN J. WALLCRAFT,  PLANNING SYSTEMS INC.,  OCTOBER 1995.
C      BASED ON EARILER VERSIONS BY SEVERAL AUTHORS.
C*
C**********
C
      EXTERNAL LINEAR, AVERMS,MINMAX
C
      REAL       ZERO,RADIAN
      PARAMETER (ZERO=0.0, RADIAN=57.2957795)
      REAL*8     DZERO,DTHIRD
      PARAMETER (DZERO=0.D0,DTHIRD=1.D0/3.D0)
C
      CHARACTER*80 CLINE
      CHARACTER*40 CCNAME
      INTEGER IUNIT,IWIT,JWIT,KWIT
      REAL*4  XFINT,YFINT,DXINT,DYINT
      REAL    HMINA,HMINB,HMAXA,HMAXB
C
      INTEGER I,IWIX,J,KREC
      REAL    XFD,YFD,DXD,DYD
      REAL    XLIN,XFDX,XOV,YOU,
     +        XMIN,XMAX,XAVE,XRMS,TFRZ
C
      INTEGER      LEN_TRIM
      CHARACTER*11 CMONTH(12)
      DATA CMONTH /  ', January',
     +               ', February',
     +               ', March',
     +               ', April',
     +               ', May',
     +               ', June',
     +               ', July',
     +               ', August',
     +               ', September',
     +               ', October',
     +               ', November',
     +               ', December'  /
C
C     STATEMENT FUNCTIONS.
C
      REAL*8 C1,C2,C3,C4,C5,C6,C7
C --- coefficients for sigma-0 (based on Brydon & Sun fit)
      DATA C1,C2,C3,C4,C5,C6,C7/
     . -1.36471E-01, 4.68181E-02, 8.07004E-01,-7.45353E-03,-2.94418E-03,
     .  3.43570E-05, 3.48658E-05/
C
      REAL*4  R4
      REAL*8  R8
      REAL*8  R,S,T
      REAL*8  A0,A1,A2,CUBQ,CUBR,CUBAN,CUBRL,CUBIM,TOFSIG,SOFSIG,SIG
C
C --- auxiliary statement for real*4 to real*8 conversion
      R8(R4)=R4
C
C --- auxiliary statements for finding root of 3rd degree polynomial
      A0(S)=(C1+C3*S)/C6
      A1(S)=(C2+C5*S)/C6
      A2(S)=(C4+C7*S)/C6
      CUBQ(S)=DTHIRD*A1(S)-(DTHIRD*A2(S))**2
      CUBR(R,S)=DTHIRD*(0.5D0*A1(S)*A2(S)-1.5D0*(A0(S)-R/C6))
     .   -(DTHIRD*A2(S))**3
C --- if q**3+r**2>0, water is too dense to yield real root at given
C --- salinitiy. setting q**3+r**2=0 in that case is equivalent to
C --- lowering sigma until a double real root is obtained.
      CUBAN(R,S)=DTHIRD*ATAN2(SQRT(MAX(DZERO,
     .   -(CUBQ(S)**3+CUBR(R,S)**2))),CUBR(R,S))
      CUBRL(R,S)=SQRT(-CUBQ(S))*COS(CUBAN(R,S))
      CUBIM(R,S)=SQRT(-CUBQ(S))*SIN(CUBAN(R,S))
C
C --- temp (deg c) as a function of sigma and salinity (mil)
      TOFSIG(R,S)=-CUBRL(R,S)+SQRT(3.)*CUBIM(R,S)-DTHIRD*A2(S)
C
C --- salinity (mil) as a function of sigma and temperature (deg c)
      SOFSIG(R,T)=(R-C1-T*(C2+T*(C4+C6*T)))/(C3+T*(C5+C7*T))
C
C --- sigma-theta as a function of temp (deg c) and salinity (mil)
C --- (friedrich-levitus 3rd degree polynomial fit)
      SIG(T,S)=(C1+C3*S+T*(C2+C5*S+T*(C4+C7*S+C6*T)))
C
C --- MODEL ARRAYS.
C
      CALL XCSPMD  !define idm,jdm
      ALLOCATE(  MSK(IDM,JDM) )
      ALLOCATE( PLON(IDM,JDM) )
      ALLOCATE( PLAT(IDM,JDM) )
      ALLOCATE(  XAF(IDM,JDM) )
      ALLOCATE(  YAF(IDM,JDM) )
      ALLOCATE(   NITM(IDM,JDM) )
      ALLOCATE(   PHOM(IDM,JDM) )
      ALLOCATE(   SILM(IDM,JDM) )
      ALLOCATE(   OXYM(IDM,JDM) )
C
C     NAMELIST INPUT.
C
      CALL ZHOPEN(6, 'FORMATTED', 'UNKNOWN', 0)
C
      CTITLE = ' '
      WRITE(6,*) 'READING /AFTITL/'
      CALL ZHFLSH(6)
      READ( 5,AFTITL)
      WRITE(6,AFTITL)
C
      ICTYPE = 1
      IFRZ   = 1
      KSIGMA = 0
      INTERP = 1
      MONTH  = 1
      ITEST  = 0
      JTEST  = 0
      JPR    = 8
      WRITE(6,*) 'READING /AFFLAG/'
      CALL ZHFLSH(6)
      READ( 5,AFFLAG)
      WRITE(6,AFFLAG)
      WRITE(6,*)
      CALL ZHFLSH(6)
C
C     GRID INPUT.
C
      CALL ZAIOST
C
      CALL ZHOPNC(21, 'regional.grid.b', 'FORMATTED', 'OLD', 0)
      CALL ZAIOPF('regional.grid.a', 'OLD', 21)
C
      READ(21,*) ! skip idm
      READ(21,*) ! skip jdm
      READ(21,*) ! skip mapflg
      READ(21,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
      CALL ZAIORD(PLON,MSK,.FALSE., HMINA,HMAXA, 21)
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (plon):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      READ(21,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
      CALL ZAIORD(PLAT,MSK,.FALSE., HMINA,HMAXA, 21)
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (plat):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      CLOSE(UNIT=21)
      CALL ZAIOCL(21)
C
C     INITIALIZE OUTPUT.
C
      CTITLE = CTITLE(1:LEN_TRIM(CTITLE)) // CMONTH(MONTH)
      WRITE(6,6000) 'OUTPUT:',CTITLE
      CALL ZHFLSH(6)
C
      CALL ZAIOPN('NEW', 10)
      CALL ZAIOPN('NEW', 11)
      CALL ZAIOPN('NEW', 12)
      CALL ZAIOPN('NEW', 13)
C
      CALL ZHOPEN(10, 'FORMATTED', 'NEW', 0)
      CALL ZHOPEN(11, 'FORMATTED', 'NEW', 0)
      CALL ZHOPEN(12, 'FORMATTED', 'NEW', 0)
      CALL ZHOPEN(13, 'FORMATTED', 'NEW', 0)
C
      PREAMBL(1) = CTITLE
      PREAMBL(2) = ' '
      PREAMBL(3) = ' '
      PREAMBL(4) = ' '
      WRITE(PREAMBL(5),'(A,2I5)')
     +        'i/jdm =',
     +       IDM,JDM
C
      PREAMBL(2) = 'Nitrate'
      WRITE(10,4101) PREAMBL
C
      PREAMBL(2) = 'Phophate'
      WRITE(11,4101) PREAMBL
C
      PREAMBL(2) = 'Silicate'
      WRITE(12,4101) PREAMBL
C
      PREAMBL(2) = 'Oxygen'
      WRITE(13,4101) PREAMBL

      IF     (KSIGMA.EQ.0) THEN
      ELSEIF (KSIGMA.EQ.2) THEN
      ELSE
        WRITE(6,*)
        WRITE(6,*) 'ERROR - KSIGMA MUST BE 0 OR 2'
        WRITE(6,*)
        STOP
      ENDIF
C
      WRITE(6,*)
      WRITE(6, 4101) PREAMBL
      WRITE(6,*)
C
C     INITIALIZE CLIMS.
C
      IF     (ICTYPE.EQ.1) THEN
        IUNIT=71
        CALL ZHOPEN(IUNIT, 'UNFORMATTED', 'OLD', 0)
        READ(IUNIT)
        CCNAME=' '
        READ(IUNIT) IWIT,JWIT,KWIT,XFINT,YFINT,DXINT,DYINT,ZLEV
        IF     (IWIT.NE.IWI .OR.
     +          JWIT.NE.JWI .OR.
     +          KWIT.NE.KWI .OR.
     +          ABS(XFINT-XFIN).GT.0.01 .OR.
     +          ABS(YFINT-YFIN).GT.0.01 .OR.
     +          ABS(DXINT-DXIN).GT.0.01 .OR.
     +          ABS(DYINT-DYIN).GT.0.01     ) THEN
          WRITE(6,9000) IUNIT,IWI,JWI,KWI,XFIN,YFIN,DXIN,DYIN
          CALL ZHFLSH(6)
          STOP
        ENDIF
        IUNIT=72
        CALL ZHOPEN(IUNIT, 'UNFORMATTED', 'OLD', 0)
        READ(IUNIT)
        CCNAME=' '
        READ(IUNIT) IWIT,JWIT,KWIT,XFINT,YFINT,DXINT,DYINT,ZLEV
        IF     (IWIT.NE.IWI .OR.
     +          JWIT.NE.JWI .OR.
     +          KWIT.NE.KWI .OR.
     +          ABS(XFINT-XFIN).GT.0.01 .OR.
     +          ABS(YFINT-YFIN).GT.0.01 .OR.
     +          ABS(DXINT-DXIN).GT.0.01 .OR.
     +          ABS(DYINT-DYIN).GT.0.01     ) THEN
          WRITE(6,9000) IUNIT,IWI,JWI,KWI,XFIN,YFIN,DXIN,DYIN
          CALL ZHFLSH(6)
          STOP
        ENDIF
        IUNIT=73
        CALL ZHOPEN(IUNIT, 'UNFORMATTED', 'OLD', 0)
        READ(IUNIT)
        CCNAME=' '
        READ(IUNIT) IWIT,JWIT,KWIT,XFINT,YFINT,DXINT,DYINT,ZLEV
        IF     (IWIT.NE.IWI .OR.
     +          JWIT.NE.JWI .OR.
     +          KWIT.NE.KWI .OR.
     +          ABS(XFINT-XFIN).GT.0.01 .OR.
     +          ABS(YFINT-YFIN).GT.0.01 .OR.
     +          ABS(DXINT-DXIN).GT.0.01 .OR.
     +          ABS(DYINT-DYIN).GT.0.01     ) THEN
          WRITE(6,9000) IUNIT,IWI,JWI,KWI,XFIN,YFIN,DXIN,DYIN
          CALL ZHFLSH(6)
          STOP
        ENDIF
        IUNIT=74
        CALL ZHOPEN(IUNIT, 'UNFORMATTED', 'OLD', 0)
        READ(IUNIT)
        CCNAME=' '
        READ(IUNIT) IWIT,JWIT,KWIT,XFINT,YFINT,DXINT,DYINT,ZLEV
        IF     (IWIT.NE.IWI .OR.
     +          JWIT.NE.JWI .OR.
     +          KWIT.NE.KWI .OR.
     +          ABS(XFINT-XFIN).GT.0.01 .OR.
     +          ABS(YFINT-YFIN).GT.0.01 .OR.
     +          ABS(DXINT-DXIN).GT.0.01 .OR.
     +          ABS(DYINT-DYIN).GT.0.01     ) THEN
          WRITE(6,9000) IUNIT,IWI,JWI,KWI,XFIN,YFIN,DXIN,DYIN
          CALL ZHFLSH(6)
          STOP
        ENDIF
      ELSEIF (ICTYPE.EQ.2) THEN
        IUNIT=71
        CALL ZHOPEN(IUNIT, 'UNFORMATTED', 'OLD', 0)
        READ(IUNIT)
        CCNAME=' '
        READ(IUNIT) IWIT,JWIT,KWIT,XFINT,YFINT,DXINT,DYINT,ZLEV
        IF     (IWIT.NE.IWI .OR.
     +          JWIT.NE.JWI .OR.
     +          KWIT.NE.KWI .OR.
     +          ABS(XFINT-XFIN).GT.0.01 .OR.
     +          ABS(YFINT-YFIN).GT.0.01 .OR.
     +          ABS(DXINT-DXIN).GT.0.01 .OR.
     +          ABS(DYINT-DYIN).GT.0.01     ) THEN
          WRITE(6,9000) IUNIT,IWI,JWI,KWI,XFIN,YFIN,DXIN,DYIN
          CALL ZHFLSH(6)
          STOP
        ENDIF
        IUNIT=73
        CALL ZHOPEN(IUNIT, 'UNFORMATTED', 'OLD', 0)
        READ(IUNIT)
        CCNAME=' '
        READ(IUNIT) IWIT,JWIT,KWIT,XFINT,YFINT,DXINT,DYINT,ZLEV
        IF     (IWIT.NE.IWI .OR.
     +          JWIT.NE.JWI .OR.
     +          KWIT.NE.KWI .OR.
     +          ABS(XFINT-XFIN).GT.0.01 .OR.
     +          ABS(YFINT-YFIN).GT.0.01 .OR.
     +          ABS(DXINT-DXIN).GT.0.01 .OR.
     +          ABS(DYINT-DYIN).GT.0.01     ) THEN
          WRITE(6,9000) IUNIT,IWI,JWI,KWI,XFIN,YFIN,DXIN,DYIN
          CALL ZHFLSH(6)
          STOP
        ENDIF
      ELSE
        IUNIT=72
        CALL ZHOPEN(IUNIT, 'UNFORMATTED', 'OLD', 0)
        READ(IUNIT)
        CCNAME=' '
        READ(IUNIT) IWIT,JWIT,KWIT,XFINT,YFINT,DXINT,DYINT,ZLEV
        IF     (IWIT.NE.IWI .OR.
     +          JWIT.NE.JWI .OR.
     +          KWIT.NE.KWI .OR.
     +          ABS(XFINT-XFIN).GT.0.01 .OR.
     +          ABS(YFINT-YFIN).GT.0.01 .OR.
     +          ABS(DXINT-DXIN).GT.0.01 .OR.
     +          ABS(DYINT-DYIN).GT.0.01     ) THEN
          WRITE(6,9000) IUNIT,IWI,JWI,KWI,XFIN,YFIN,DXIN,DYIN
          CALL ZHFLSH(6)
          STOP
        ENDIF
        IUNIT=73
        CALL ZHOPEN(IUNIT, 'UNFORMATTED', 'OLD', 0)
        READ(IUNIT)
        CCNAME=' '
        READ(IUNIT) IWIT,JWIT,KWIT,XFINT,YFINT,DXINT,DYINT,ZLEV
        IF     (IWIT.NE.IWI .OR.
     +          JWIT.NE.JWI .OR.
     +          KWIT.NE.KWI .OR.
     +          ABS(XFINT-XFIN).GT.0.01 .OR.
     +          ABS(YFINT-YFIN).GT.0.01 .OR.
     +          ABS(DXINT-DXIN).GT.0.01 .OR.
     +          ABS(DYINT-DYIN).GT.0.01     ) THEN
          WRITE(6,9000) IUNIT,IWI,JWI,KWI,XFIN,YFIN,DXIN,DYIN
          CALL ZHFLSH(6)
          STOP
        ENDIF
      ENDIF
C
C     DEFINE THE GRID COORDINATES.
C
      IF     (IWI*DXIN.GE.359.9) THEN
        IF     (ABS(IWI * DXIN - 360.0) .GT. 0.01) THEN
          WRITE(6,9050)
          CALL ZHFLSH(6)
          STOP
        ENDIF
        IWIX = IWI + 1
        IBD(1) = 3
        IBD(2) = 3
        IBD(3) = 2
        IBD(4) = 2
      ELSE
        IWIX = IWI
        IBD(1) = 2
        IBD(2) = 2
        IBD(3) = 2
        IBD(4) = 2
      ENDIF
C
C     CONVERT HYCOM LON,LAT TO CLIMATOLOGY ARRAY COORDS.
C
      XLIN  = XFIN + (IWIX-1)*DXIN
      XAMIN = 2*IWI
      XAMAX = 0
      DO J= 1,JDM
        DO I= 1,IDM
          XOV = PLON(I,J)
          IF     (XOV.LT.XFIN) THEN
            XOV = XOV + 360.0
          ELSEIF (XOV.GE.XLIN) THEN
            XOV = XOV - 360.0
          ENDIF
C
          XAF(I,J) = 3.0 + (XOV - XFIN)/DXIN
C
          IF     (MOD(J,100).EQ.1 .OR. J.EQ.JDM) THEN
            IF     (MOD(I,10).EQ.1 .OR. I.EQ.IDM) THEN
              WRITE(6,'("I,J,LONV,XAF =",2I5,2F10.3)') I,J,XOV,XAF(I,J)
            ENDIF
          ENDIF
          XAMIN  = MIN( XAMIN, XAF(I,J) )
          XAMAX  = MAX( XAMAX, XAF(I,J) )
        ENDDO
      ENDDO
C
      YAMIN = 2*JWI
      YAMAX = 0
      DO I= 1,IDM
        DO J= 1,JDM
          YOU = PLAT(I,J)
C
          YAF(I,J) = 3.0 + (YOU - YFIN)/DYIN
C
          IF     (MOD(I,100).EQ.1 .OR. I.EQ.IDM) THEN
            IF     (MOD(J,10).EQ.1 .OR. J.EQ.JDM) THEN
              WRITE(6,'("I,J,LATU,YAF =",2I5,2F10.3)') I,J,YOU,YAF(I,J)
            ENDIF
          ENDIF
          YAMIN  = MIN( YAMIN, YAF(I,J) )
          YAMAX  = MAX( YAMAX, YAF(I,J) )
        ENDDO
      ENDDO
C
      WRITE(6,6200) XAMIN,XAMAX,YAMIN,YAMAX
      CALL ZHFLSH(6)
C
C     CHECK THAT THE INTERPOLATION IS 'SAFE',
C
      IF     (INT(XAMIN).LT.3 .OR. INT(XAMAX).GT.IWI+2 .OR.
     +        INT(YAMIN).LT.3 .OR. INT(YAMAX).GT.JWI+2     ) THEN
        WRITE(6,9150)
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
C     PROCESS ALL THE CLIM RECORDS.
C
c      DO J= 1,JDM
c        DO I= 1,IDM
c          RT(I,J) =  ZERO
c        ENDDO
c      ENDDO
C
      DO 810 KREC= 1,KWI
C
C       READ THE INPUT CLIMS.
C
        IF     (ICTYPE.EQ.1) THEN
          READ(71) NITSEAIN
                   NITSEAIN = NITSEAIN * 6.625 * 12.01
          READ(72) PHOSEAIN
                   PHOSEAIN = PHOSEAIN * 30.97 * 12.01 
          READ(73) SILSEAIN
                   SILSEAIN = SILSEAIN * 6.625 * 12.01
          READ(74) OXYSEAIN
                   OXYSEAIN = OXYSEAIN * 32.0
        ENDIF

C AS17072012: Missing values are negative - set to zero
        DO  J= 1,JWI
          DO  I= 1,IWI
            if (NITSEAIN(I,J).lt.0.0) NITSEAIN(I,J)=0.0
            if (PHOSEAIN(I,J).lt.0.0) PHOSEAIN(I,J)=0.0
            if (SILSEAIN(I,J).lt.0.0) SILSEAIN(I,J)=0.0
            if (OXYSEAIN(I,J).lt.0.0) OXYSEAIN(I,J)=0.0
         end do
        end do
C
C       COPY INTO THE (LARGER) INTERPOLATION ARRAYS.
C
        DO 310 J= 1,JWI
          DO 311 I= 1,IWI
            NITSEAI(I+2,J+2) = dble(NITSEAIN(I,J))
            PHOSEAI(I+2,J+2) = dble(PHOSEAIN(I,J))
            SILSEAI(I+2,J+2) = dble(SILSEAIN(I,J))
            OXYSEAI(I+2,J+2) = dble(OXYSEAIN(I,J))
  311     CONTINUE
  310   CONTINUE
C
C       FILL IN THE PADDING AREA AS NECESSARY.
C
        IF     (INT(XAMAX).GE.IWI+1) THEN
          IF     (IWIX.GT.IWI) THEN
            DO 320 J= 3,JWI+2
              NITSEAI(IWI+3,J) = NITSEAI(3,J)
              NITSEAI(IWI+4,J) = NITSEAI(4,J)
              PHOSEAI(IWI+3,J) = PHOSEAI(3,J)
              PHOSEAI(IWI+4,J) = PHOSEAI(4,J)
              SILSEAI(IWI+3,J) = SILSEAI(3,J)
              SILSEAI(IWI+4,J) = SILSEAI(4,J)
              OXYSEAI(IWI+3,J) = OXYSEAI(3,J)
              OXYSEAI(IWI+4,J) = OXYSEAI(4,J)
  320       CONTINUE
          ELSE
            DO 325 J= 3,JWI+2
              NITSEAI(IWI+3,J) = 2.0*NITSEAI(IWI+2,J) -     
     +                               NITSEAI(IWI+1,J)
              NITSEAI(IWI+4,J) = 3.0*NITSEAI(IWI+2,J) - 
     +                           2.0*NITSEAI(IWI+1,J)
              PHOSEAI(IWI+3,J) = 2.0*PHOSEAI(IWI+2,J) -     
     +                               PHOSEAI(IWI+1,J)
              PHOSEAI(IWI+4,J) = 3.0*PHOSEAI(IWI+2,J) - 
     +                           2.0*PHOSEAI(IWI+1,J)
              SILSEAI(IWI+3,J) = 2.0*SILSEAI(IWI+2,J) -     
     +                               SILSEAI(IWI+1,J)
              SILSEAI(IWI+4,J) = 3.0*SILSEAI(IWI+2,J) - 
     +                           2.0*SILSEAI(IWI+1,J)
              OXYSEAI(IWI+3,J) = 2.0*OXYSEAI(IWI+2,J) -     
     +                               OXYSEAI(IWI+1,J)
              OXYSEAI(IWI+4,J) = 3.0*OXYSEAI(IWI+2,J) - 
     +                           2.0*OXYSEAI(IWI+1,J)
  325       CONTINUE
          ENDIF
        ENDIF
        IF     (INT(XAMIN).LE.3) THEN
          IF     (IWIX.GT.IWI) THEN
            DO 330 J= 3,JWI+2
              NITSEAI(1,J) = NITSEAI(IWI+1,J)
              NITSEAI(2,J) = NITSEAI(IWI+2,J)
              PHOSEAI(1,J) = PHOSEAI(IWI+1,J)
              PHOSEAI(2,J) = PHOSEAI(IWI+2,J)
              SILSEAI(1,J) = SILSEAI(IWI+1,J)
              SILSEAI(2,J) = SILSEAI(IWI+2,J)
              OXYSEAI(1,J) = OXYSEAI(IWI+1,J)
              OXYSEAI(2,J) = OXYSEAI(IWI+2,J)
  330       CONTINUE
          ELSE
            DO 335 J= 3,JWI+2
              NITSEAI(1,J) = 3.0*NITSEAI(3,J) - 2.0*NITSEAI(4,J)
              NITSEAI(2,J) = 2.0*NITSEAI(3,J) -     NITSEAI(4,J)
              PHOSEAI(1,J) = 3.0*PHOSEAI(3,J) - 2.0*PHOSEAI(4,J)
              PHOSEAI(2,J) = 2.0*PHOSEAI(3,J) -     PHOSEAI(4,J)
              SILSEAI(1,J) = 3.0*SILSEAI(3,J) - 2.0*SILSEAI(4,J)
              SILSEAI(2,J) = 2.0*SILSEAI(3,J) -     SILSEAI(4,J)
              OXYSEAI(1,J) = 3.0*OXYSEAI(3,J) - 2.0*OXYSEAI(4,J)
              OXYSEAI(2,J) = 2.0*OXYSEAI(3,J) -     OXYSEAI(4,J)
  335       CONTINUE
          ENDIF
        ENDIF
        IF     (INT(YAMAX).GE.JWI+1) THEN
          DO 340 I= 1,IWI+4
            NITSEAI(I,JWI+3) = 2.0*NITSEAI(I,JWI+2) -     
     +                             NITSEAI(I,JWI+1)
            NITSEAI(I,JWI+4) = 3.0*NITSEAI(I,JWI+2) - 
     +                         2.0*NITSEAI(I,JWI+1)
            PHOSEAI(I,JWI+3) = 2.0*PHOSEAI(I,JWI+2) -     
     +                             PHOSEAI(I,JWI+1)
            PHOSEAI(I,JWI+4) = 3.0*PHOSEAI(I,JWI+2) - 
     +                         2.0*PHOSEAI(I,JWI+1)
            SILSEAI(I,JWI+3) = 2.0*SILSEAI(I,JWI+2) -     
     +                             SILSEAI(I,JWI+1)
            SILSEAI(I,JWI+4) = 3.0*SILSEAI(I,JWI+2) - 
     +                         2.0*SILSEAI(I,JWI+1)
            OXYSEAI(I,JWI+3) = 2.0*OXYSEAI(I,JWI+2) -     
     +                             OXYSEAI(I,JWI+1)
            OXYSEAI(I,JWI+4) = 3.0*OXYSEAI(I,JWI+2) - 
     +                         2.0*OXYSEAI(I,JWI+1)
  340     CONTINUE
        ENDIF
        IF     (INT(YAMIN).LE.3) THEN
          DO 350 I= 1,IWI+4
            NITSEAI(I,1) = 3.0*NITSEAI(I,3) - 2.0*NITSEAI(I,4)
            NITSEAI(I,2) = 2.0*NITSEAI(I,3) -     NITSEAI(I,4)
            PHOSEAI(I,1) = 3.0*PHOSEAI(I,3) - 2.0*PHOSEAI(I,4)
            PHOSEAI(I,2) = 2.0*PHOSEAI(I,3) -     PHOSEAI(I,4)
            SILSEAI(I,1) = 3.0*SILSEAI(I,3) - 2.0*SILSEAI(I,4)
            SILSEAI(I,2) = 2.0*SILSEAI(I,3) -     SILSEAI(I,4)
            OXYSEAI(I,1) = 3.0*OXYSEAI(I,3) - 2.0*OXYSEAI(I,4)
            OXYSEAI(I,2) = 2.0*OXYSEAI(I,3) -     OXYSEAI(I,4)
  350     CONTINUE
        ENDIF
C
C       INTERPOLATE FROM NATIVE TO MODEL CLIM GRIDS.
C       ALSO INFORCE A STABLE DENSITY PROFILE.
C
        IF     (IFRZ.EQ.0) THEN
C         IGNORE ICE.
          TFRZ = -HUGE(TFRZ)
        ELSEIF (IFRZ.EQ.1) THEN
C         ASSUME ICE FORMS (I.E. MIN SST) AT -1.8 DEGC.
          TFRZ = -1.8
        ENDIF
        IF     (ICTYPE.EQ.1) THEN
          IF     (INTERP.EQ.0) THEN
            CALL LINEAR(NITM,XAF,YAF,IDM,IDM,JDM,
     +                  NITSEAI,IWI+4,IWI+4,JWI+4)
            CALL LINEAR(PHOM,XAF,YAF,IDM,IDM,JDM,
     +                  PHOSEAI,IWI+4,IWI+4,JWI+4)
            CALL LINEAR(SILM,XAF,YAF,IDM,IDM,JDM,
     +                  SILSEAI,IWI+4,IWI+4,JWI+4)
            CALL LINEAR(OXYM,XAF,YAF,IDM,IDM,JDM,
     +                  OXYSEAI,IWI+4,IWI+4,JWI+4)
          ELSE
            CALL CUBSPL(NITM,XAF,YAF,IDM,IDM,JDM,
     +                  NITSEAI,IWI+4,IWI+4,JWI+4, 
     +                  IBD, FXI,FYI,WQSEA3,WK)
            CALL CUBSPL(PHOM,XAF,YAF,IDM,IDM,JDM,
     +                  PHOSEAI,IWI+4,IWI+4,JWI+4, 
     +                  IBD, FXI,FYI,WQSEA3,WK)
            CALL CUBSPL(SILM,XAF,YAF,IDM,IDM,JDM,
     +                  SILSEAI,IWI+4,IWI+4,JWI+4, 
     +                  IBD, FXI,FYI,WQSEA3,WK)
            CALL CUBSPL(OXYM,XAF,YAF,IDM,IDM,JDM,
     +                  OXYSEAI,IWI+4,IWI+4,JWI+4, 
     +                  IBD, FXI,FYI,WQSEA3,WK)
          ENDIF
        ENDIF
C
C       WRITE OUT STATISTICS.
C
        CALL MINMAX(NITM,IDM,JDM, XMIN,XMAX)
        CALL AVERMS(NITM,IDM,JDM, XAVE,XRMS)
        WRITE(6,8100) 'NITSEA', XMIN,XMAX,XAVE,XRMS
C
        CALL MINMAX(PHOM,IDM,JDM, XMIN,XMAX)
        CALL AVERMS(PHOM,IDM,JDM, XAVE,XRMS)
        WRITE(6,8100) 'PHOSEA', XMIN,XMAX,XAVE,XRMS
C
        CALL MINMAX(SILM,IDM,JDM, XMIN,XMAX)
        CALL AVERMS(SILM,IDM,JDM, XAVE,XRMS)
        WRITE(6,8100) 'SILSEA', XMIN,XMAX,XAVE,XRMS

        CALL MINMAX(OXYM,IDM,JDM, XMIN,XMAX)
        CALL AVERMS(OXYM,IDM,JDM, XAVE,XRMS)
        WRITE(6,8100) 'OXYSEA', XMIN,XMAX,XAVE,XRMS
C
C       DIAGNOSTIC PRINTOUT.
C
        IF     (MIN(ITEST,JTEST).GT.0) THEN
          WRITE(6,*)
          WRITE(6,'(A,2I5,I3,A,F8.2,A,3F6.2)')
     +      'I,J,K =',ITEST,JTEST,KREC,
     +      '   ZLEV =',ZLEV(KREC),
     +      '   NIT,PHO,SIL,OXY =',NITM(ITEST,JTEST),
     +                             PHOM(ITEST,JTEST),
     +                             SILM(ITEST,JTEST),
     +                             OXYM(ITEST,JTEST)
          WRITE(6,*)
          CALL ZHFLSH(6)
        ENDIF
C
C       WRITE OUT HYCOM CLIMS.
C
        CALL ZAIOWR(NITM,MSK,.FALSE., XMIN,XMAX, 10, .FALSE.)
        WRITE(10,4102) '              nitrate',ZLEV(KREC),XMIN,XMAX
C
        CALL ZAIOWR(PHOM,MSK,.FALSE., XMIN,XMAX, 11, .FALSE.)
        WRITE(11,4102) '            phosphate',ZLEV(KREC),XMIN,XMAX
C
        CALL ZAIOWR(SILM,MSK,.FALSE., XMIN,XMAX, 12, .FALSE.)
        WRITE(12,4102) '             silicate',ZLEV(KREC),XMIN,XMAX
C
        CALL ZAIOWR(OXYM,MSK,.FALSE., XMIN,XMAX, 13, .FALSE.)
        WRITE(13,4102) '               oxygen',ZLEV(KREC),XMIN,XMAX
C
        WRITE(6,6300) KREC,ZLEV(KREC)
        CALL ZHFLSH(6)
  810 CONTINUE
C
      CALL ZAIOCL(10)
      CLOSE( UNIT=10)
      CALL ZAIOCL(11)
      CLOSE( UNIT=11)
      CALL ZAIOCL(12)
      CLOSE( UNIT=12)
      CALL ZAIOCL(13)
      CLOSE( UNIT=13)
C
C     SUMMARY.
C
      WRITE(6,6400) KWI
      CALL ZHFLSH(6)
      STOP
C
 4101 FORMAT(A79)
 4102 FORMAT(A,': depth,range = ',F7.1,1P2E16.7)
 6000 FORMAT(1X,A,2X,A40 //)
 6200 FORMAT(/ 1X,'MIN,MAX I COORDS = ',F8.2,',',F8.2 
     +       / 1X,'MIN,MAX J COORDS = ',F8.2,',',F8.2 /)
 6300 FORMAT(10X,'WRITING CLIM RECORD',I3,'     ZLEV =',F7.1 /)
 6400 FORMAT(I5,' LEVEL CLIMATOLOGY COMPLETED.')
 8100 FORMAT(1X,A,': MIN=',F13.8,' MAX=',F13.8,
     +             ' AVE=',F13.8,' RMS=',F13.8)
 9000 FORMAT(// 20X,'*****  ERROR ON UNIT -',I3,
     +   ' INPUT DOES NOT AGREE WITH PARAMETERS  *****' //
     +   1X,'IWI,JWI,KWI = ',I5,  I10,  I4 /
     +   1X,'XFIN,YFIN   = ',F8.2,F10.2    /
     +   1X,'DXIN,DYIN   = ',F9.3, F9.3    //)
 9050 FORMAT(// 20X,'**********  ERROR - ',
     +   'INPUT IS GLOBAL AND IWI*DXIN IS NOT 360 DEGREES  ****' //)
 9150 FORMAT(// 20X,'**********  ERROR - ',
     +   'THE OUTPUT GRID IS NOT INSIDE THE INPUT GRID  **********' //)
C     END OF PROGRAM WNDINT.
      END



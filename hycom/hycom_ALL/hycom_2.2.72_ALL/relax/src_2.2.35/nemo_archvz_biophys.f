      PROGRAM NEMO_ARCHVZ
      USE MOD_ZA  ! HYCOM array I/O interface
!      use mod_xc  ! HYCOM communication interface

      IMPLICIT NONE
C
C     DIAGNOSTIC/DEBUGGING VARIABLES.
C
*     LOGICAL, PARAMETER :: LDEBUG = .FALSE.
      LOGICAL, PARAMETER :: LDEBUG = .TRUE.
C
      INTEGER        ITEST,JTEST
      COMMON/DEBUGI/ ITEST,JTEST
      SAVE  /DEBUGI/
C
C     BLKDAT VARIABLES.
C
      CHARACTER*79 CTITLE
      CHARACTER*40 SIGFMT
      INTEGER      IVERSN,IEXPT,YRFLAG,KDM,LEVTOP,MONTH,SIGVER,
     +             NHYBRD,NSIGMA,THFLAG
      INTEGER      JDW
      REAL*4       DP00S,DP00,DP00X,DP00F,DS00,DS00X,DS00F,ISOTOP,
     +             DP00I,
     +             SIGMA(9999),THBASE,THKMIN,BLK
      LOGICAL      VSIGMA,ISOPYC
C
C     I/O VARIABLES.
C
      CHARACTER PREAMBL(5)*79,CLINE*80
      REAL*4    HMINA,HMINB,HMAXA,HMAXB,ZLEVK
C
C     INPUT CLIM ARRAYS.
C
      REAL*4               :: ZLEV(999)
      REAL*4,  ALLOCATABLE :: TZ(:,:,:),SZ(:,:,:),RZ(:,:,:)
C
C     OUTPUT ARRAYS
C
      REAL*4,  ALLOCATABLE :: PMIX(:,:),PM(:,:),TM(:,:),SM(:,:),RM(:,:)
C
C     OTHER ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: DEPTH(:,:),PLAT(:,:),
     +                        RKM1(:,:),PKM1(:,:),PKM2(:,:),
     +                        PCM0(:,:),PCM1(:,:),PCM2(:,:),
     +                        PZBOT(:,:),PZTOP(:,:),
     +                        WORK(:,:),
     +                        SIG3D(:,:,:),SIG3D_TMP(:)
      REAL*4               :: DP0K(  9999),DS0K(  9999),
     +                        DPCK(0:9999),DSCK(0:9999)

      REAL*4 :: dummy
      integer :: iosnstep,nrec,irec,tlevel
Cmostafa      integer, allocatable :: coord(:), tlevel1(:)
      integer       :: coord(999), tlevel1(999)

      character(len=5) :: char5
Cmostafa      character(len=8), allocatable :: cfld(:)
      character(len=8)    :: cfld(999)

      REAL*4 :: rday,dens

      REAL*4,  ALLOCATABLE :: TZ1(:,:),SZ1(:,:),RZ1(:,:)
      REAL*4,  ALLOCATABLE :: TH1(:,:), TH(:,:,:)
      REAL*4,  ALLOCATABLE :: VV1(:,:), VV(:,:,:)
      REAL*4,  ALLOCATABLE :: UV1(:,:), UV(:,:,:)
      REAL*4,  ALLOCATABLE :: ZL(:)
      REAL*4,  ALLOCATABLE :: UM(:,:), VM(:,:)
      REAL*4,  ALLOCATABLE :: PU(:,:), PV(:,:)

      REAL*4,  ALLOCATABLE :: MTG1(:,:)
      REAL*4,  ALLOCATABLE :: SRFH(:,:)
      REAL*4,  ALLOCATABLE :: SUFLX(:,:)
      REAL*4,  ALLOCATABLE :: SAFLX(:,:)
      REAL*4,  ALLOCATABLE :: BLDP(:,:)
      REAL*4,  ALLOCATABLE :: MIXDP(:,:)
      REAL*4,  ALLOCATABLE :: UBRTP(:,:)
      REAL*4,  ALLOCATABLE :: VBRTP(:,:)
C
C**********
C*
C 1)  FROM A Z-LEVEL CLIMATOLOGY ON THE HYCOM REGION GRID,
C      CREATE A 'ISOPYCNAL' CLIMATOLOGY SUITABLE FOR INPUT TO HYCOM.
C     THIS VERSION USES LINEAR INTERPOLATION BETWEEN Z-LEVELS, AND
C      ISOPYCNALS HALF WAY BETWEEN TARGET DENSITIES.
C
C      ONLY FOR USE WITH HYCOM 2.2.58 OR LATER.
C
C 2)  PARAMETERS:
C
C     MODEL GRID SPECIFICATION (W.R.T. PRESSURE GRID):
C
C        IDM    = 1ST DIMENSION OF MAJOR (HYCOM) MODEL ARRAYS
C        JDM    = 2ND DIMENSION OF MAJOR (HYCOM) MODEL ARRAYS
C        KZ     = NUMBER OF Z-LEVELS IN INPUT CLIMATOLOGY
C
C 3)  INPUT:
C        ON UNIT 51:  BATHYMETRY FILE
C        ON UNIT 51A: BATHYMETRY FILE
C        ON UNIT 52:  SPACIALLY VARYING ISOPYCNAL TARGET DENSITY FILE
C        ON UNIT 52A: SPACIALLY VARYING ISOPYCNAL TARGET DENSITY FILE
C        ON UNIT 72:  TEMP. Z-LEVEL CLIM FILE
C        ON UNIT 73:  SALN. Z-LEVEL CLIM FILE
C        ON UNIT 99:  A SUBSET OF blkdat.input FROM TARGET SIMULATION
c --- 'flnm_i' = name of original  archive file
c --- 'flnm_o' = name of target    archive file
C        'knemo' = number   of nemo levels
C        'flnm_z' = name of nemo cell interface depths (text file), or "NONE"
c --- 'flag_t' =  temperature flage, or "NONE" no conversion
c --- 'flag_s' =  salinity flage, or "NONE" no conversion
c --- 'flag_u' =  u-velocity flage, or "NONE" no conversion
c --- 'flag_v' =  v-velocity flage, or "NONE" no conversion

C     OUTPUT:
C        ON UNIT 10:  TEMP. LAYER CLIM FILE
C        ON UNIT 10A: TEMP. LAYER CLIM FILE
C        ON UNIT 11:  SALN. LAYER CLIM FILE
C        ON UNIT 11A: SALN. LAYER CLIM FILE
C        ON UNIT 12:  INTF. LAYER CLIM FILE
C        ON UNIT 12A: INTF. LAYER CLIM FILE
C        ON UNIT 21:  DUMMY HYCOM ARCHIVE FILE
C        ON UNIT 21A: DUMMY HYCOM ARCHIVE FILE
C
C 4)  ALAN J. WALLCRAFT,  NAVAL RESEARCH LABORATORY, MARCH 2001
C                                                 AND JUNE 2009.
C     MOSTAFA BAKHODAY-PASKYABI, NERSC, Bergen, Norway, September 2017
C    
C
C**********
C
      REAL*4     ZERO,ONE,RADIAN
      PARAMETER (ZERO=0.0, ONE=1.0, RADIAN=57.2957795)
C
      INTEGER I,J,K,KK,KZ,KZTOP,L
C --- MOSTAFA: BEGIN
      real*4        :: MODEL_DAY
      INTEGER   KNEMO,NRECL,indx
      CHARACTER*240 flnm_z,flnm_o
      CHARACTER*40 flag_t,flag_s,flag_u,flag_v,flag_th
      REAL*4    onem,spval
      PARAMETER (spval=2.0**100)
      REAL*4,    parameter   :: hspval=0.5*2.0**100  ! half spval
      CHARACTER*40 flag_no3,flag_po4,flag_si
      REAL*4,  ALLOCATABLE :: NO3(:,:,:),PO4(:,:,:),SI(:,:,:)
      REAL*4,  ALLOCATABLE :: NO3M(:,:),PO4M(:,:),SIM(:,:)
      integer  dummy_index
      CHARACTER*80 dummy_name

C --- MOSTAFA: END
      REAL*4  DP0KF,DPMS,DS0KF,DSMS,SZK,RZK,TZK,TIME,THK,THIKMN,
     +        PZMID,RZLOC,DMIN,QDEP,
     +        PAVE,XAVE,XMAX,XMIN,TMIN,SMIN,PMIN,ZJ,ZZ,Q,
     +        PINTEG,SIGMAA,SIGMAB,ZBOT,ZTOP
C
      REAL*4  SIG_V,SOFSIG_V,TOFSIG_V
      CHARACTER*512 archvfile , archvfileo
      CHARACTER*80 fldname
      real   :: lrdens
c     INTEGER, EXTERNAL  :: IARGC
      INTEGER, INTRINSIC  :: IARGC

      INTEGER IOS,ISTEP
C
      CALL XCSPMD
C
      onem  = 9806.0   ! g/thref

      if (IARGC() == 1)  THEN
         call getarg(1,archvfile)
         i = index(archvfile,".",back=.true.)

         ! Strip .a or .v ending
         if (i > 0 ) then
            archvfile = archvfile(1:i-1)
         end if
         print '(a)',"Archive file:"//trim(archvfile)
      else 
         read (*,'(a)') archvfile
         write (6,'(2a)') ' input file: ',trim(archvfile)
         i = index(archvfile,".",back=.true.)
         ! Strip .a or .v ending
         if (i > 0 ) then
             archvfile = archvfile(1:i-1)
         end if
         print '(a)',"Archive file:"//trim(archvfile)

         read (*,'(a)') archvfileo
         write (6,'(2a)') ' input file: ',trim(archvfileo)
         i = index(archvfileo,".",back=.true.)
         ! Strip .a or .v ending
         if (i > 0 ) then
            archvfileo = archvfileo(1:i-1)
         end if
         print '(a)',"Archive file:"//trim(archvfileo)

         read (*,'(i6)') KNEMO
         write (6,'(a,i6)') 'MERCATOR number of layers: ',KNEMO

         read (*,'(a)') flnm_z
         write (6,'(2a)') 'MERCATOR depth file: ',trim(flnm_z)
         write(6,*)

         read (*,'(a)') flag_t
         write (6,'(2a)') 'Temperature conversion: ',flag_t

         read (*,'(a)') flag_s
         write (6,'(2a)') 'Salinity conversion: ',flag_s

         read (*,'(a)') flag_th
         write (6,'(2a)') 'Thickness conversion: ',flag_th

         read (*,'(a)') flag_u
         write (6,'(2a)') 'U-velocity conversion: ',flag_u

         read (*,'(a)') flag_v
         write (6,'(2a)') 'V-velocity conversion: ',flag_v

         read (*,'(a)') flag_no3
         write (6,'(2a)') 'BIO: NO3: ',flag_no3

         read (*,'(a)') flag_po4
         write (6,'(2a)') 'BIO: PO4: ',flag_po4

         read (*,'(a)') flag_si
         write (6,'(2a)') 'BIO: SI: ',flag_si


         call flush(6)

         ! exit(1)
      end if

      CALL INIT_SPEC(ARCHVFILE,NREC,KZ,CFLD,COORD,TLEVEL1,MODEL_DAY)
      IF (KZ.NE.KNEMO) THEN
        WRITE(6,*) KZ,KNEMO
        WRITE(6,*) 'ERROR-wrong KZ (archv.[b] and KNEMO).'
        WRITE(6,*)
        CALL ZHFLSH(6)
        STOP
      ENDIF
C      ALLOCATE( CFLD(NREC) )
C      ALLOCATE( TLEVEL1(NREC) )
C      ALLOCATE( COORD(NREC) )
C      CALL FIELDS_SPEC(ARCHVFILE,CFLD,COORD,TLEVEL1,NREC)
C
      ALLOCATE( TZ(KZ+1,IDM,JDM) ) !
      ALLOCATE( SZ(KZ+1,IDM,JDM) )
      ALLOCATE( RZ(KZ+1,IDM,JDM) )
      ALLOCATE(      ZL(KZ+1) )    ! THIS IS MY ADDITION TO THE CODE @CAGLAR
      ALLOCATE(      UM(IDM,JDM) ) !
      ALLOCATE(      VM(IDM,JDM) ) !
      ALLOCATE(      PU(IDM,JDM) ) !
      ALLOCATE(      PV(IDM,JDM) ) !
      ALLOCATE(    PMIX(IDM,JDM) )
      ALLOCATE(      PM(IDM,JDM) )
      ALLOCATE(      TM(IDM,JDM) )
      ALLOCATE(      SM(IDM,JDM) )
      ALLOCATE(      RM(IDM,JDM) )
      ALLOCATE(   DEPTH(IDM,JDM) )
      ALLOCATE(    PLAT(IDM,JDM) )
      ALLOCATE(    RKM1(IDM,JDM) )
      ALLOCATE(    PKM1(IDM,JDM) )
      ALLOCATE(    PKM2(IDM,JDM) )
      ALLOCATE(    PCM0(IDM,JDM) )
      ALLOCATE(    PCM1(IDM,JDM) )
      ALLOCATE(    PCM2(IDM,JDM) )
      ALLOCATE(   PZBOT(IDM,JDM) )
      ALLOCATE(   PZTOP(IDM,JDM) )
      ALLOCATE(    WORK(IDM,JDM) )
      ALLOCATE(     MSK(IDM,JDM) )
!Mostafa TODO: need to becoem more clever in distingishing between phy &
!bio variables
      if     (flag_no3.ne."NONE") then
          ALLOCATE(   NO3M(IDM,JDM) )
          ALLOCATE(   NO3(KZ+1,IDM,JDM) )
      endif
      if     (flag_po4.ne."NONE") then
          ALLOCATE(    PO4M(IDM,JDM) )
          ALLOCATE(    PO4(KZ+1,IDM,JDM) )
      endif
      if     (flag_si.ne."NONE") then
          ALLOCATE(     SIM(IDM,JDM) )
          ALLOCATE(     SI(KZ+1,IDM,JDM) )
      endif



C
C     BLKDAT INPUT.
C
      CALL ZHOPEN(99, 'FORMATTED', 'OLD', 0)
C
C --- 'ctitle' = climatology title
C
      WRITE(6,*)
      READ(99,'(A79)') CTITLE
      WRITE(6,'(A79)') CTITLE
C
C --- 'month'  = month (1 to 12)
C --- 'sigver' = version of the equation of state
C --- 'levtop' = top level of input clim. to use (optional, default 1)
C --- 'iversn' = hycom version number x10
C --- 'iexpt'  = experiment number x10
C --- 'yrflag' = days in year flag (0=360,1=366,2=366Jan1)
C --- 'mapflg' = map flag (0=mercator,1=rotated,2=uniform,3=beta-plane)
C --- 'idm   ' = longitudinal array size
C --- 'jdm   ' = latitudinal  array size
C --- 'jdw   ' = width of zonal average (optional, default 0)
C --- 'itest ' = grid point where detailed diagnostics are desired
C --- 'jtest ' = grid point where detailed diagnostics are desired
C --- 'kdm   ' = longitudinal array size
C --- 'nhybrd' = number of hybrid levels (0=all isopycnal)
C --- 'nsigma' = number of sigma  levels (nhybrd-nsigma z-levels)
C
      WRITE(6,*)
      CALL BLKINI(SIGVER, 'sigver')
      CALL BLKINI2(I,J,  'levtop','iversn')
      IF     (J.EQ.1) THEN
        LEVTOP = I
        CALL BLKINI(IVERSN,'iversn')
      ELSE
        IVERSN = I
        LEVTOP = 1
        write(6,'(a6," =",i6)') 'levtop',LEVTOP
      ENDIF
      WRITE(6,*)
      CALL BLKINI(IEXPT, 'iexpt ')
      CALL BLKINI(YRFLAG,'yrflag')
      WRITE(6,*)
      CALL BLKINI(I,     'idm   ')
      CALL BLKINI(J,     'jdm   ')
C
      IF     (I.NE.IDM) THEN
        WRITE(6,*)
        WRITE(6,*) 'ERROR - wrong IDM'
        WRITE(6,*)
        CALL ZHFLSH(6)
        STOP
      ELSEIF (J.NE.JDM) THEN
        WRITE(6,*)
        WRITE(6,*) 'ERROR - wrong JDM'
        WRITE(6,*)
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      CALL BLKINI2(I,J, 'jdw   ','itest ')
      IF     (J.EQ.1) THEN !jdw
        JDW = I
        CALL BLKINI(ITEST, 'itest ')
      ELSE !itest
        ITEST = I
        JDW   = 0
      ENDIF
      CALL BLKINI(JTEST, 'jtest ')
      CALL BLKINI(KDM,   'kdm   ')
      CALL BLKINI(NHYBRD,'nhybrd')
      CALL BLKINI(NSIGMA,'nsigma')
C
      ISOPYC = NHYBRD.EQ.0
      write(*,*)ISOPYC,'ISOPYC'
      write(*,*)NHYBRD,'NHYBRD'
C
      IF     (IVERSN.LE.20) THEN
C
C ---   'dp00s'  = sigma   spacing minimum thickness (m)
C ---   'dp00'   = z-level spacing minimum thickness (m)
C ---   'dp00x'  = z-level spacing maximum thickness (m)
C ---   'dp00f'  = z-level spacing stretching factor (1.0=const.spacing)
C
        CALL BLKINR(DP00S, 'dp00s ','(a6," =",f10.4," m")')
        CALL BLKINR(DP00,  'dp00  ','(a6," =",f10.4," m")')
        CALL BLKINR(DP00X, 'dp00x ','(a6," =",f10.4," m")')
        CALL BLKINR(DP00F, 'dp00f ','(a6," =",f10.4," ")')
        DS00   = MAX(DP00S,0.01)
        DS00X  = DS00
        DS00F  = 1.0
        IF     (ISOPYC) THEN
          ISOTOP = -1.0  !turn off isotop
        ELSE
          ISOTOP = 0.01  !first layer always fixed
        ENDIF
      ELSE
C
C ---   'isotop' = shallowest depth for isopycnal layers (m), optional
C ---   'dp00'   = deep    z-level spacing minimum thickness (m)
C ---   'dp00x'  = deep    z-level spacing maximum thickness (m)
C ---   'dp00f'  = deep    z-level spacing stretching factor (1.0=const.z)
C ---   'ds00'   = shallow z-level spacing minimum thickness (m)
C ---   'ds00x'  = shallow z-level spacing maximum thickness (m)
C ---   'ds00f'  = shallow z-level spacing stretching factor (1.0=const.z)
C ---   'dp00i'  = deep iso-pycnal spacing minimum thickness (m)
C
C ---   the above describe a system that is isopycnal or:
C ---       z in    deep water, based on dp00,dp00x,dp00f
C ---       z in shallow water, based on ds00,ds00x,ds00f and nsigma
C ---       sigma between them, based on ds00,ds00x,ds00f and nsigma
C
C ---   away from the surface, the minimum layer thickness is dp00i.
C
C ---   for z-only set nsigma=0 (and ds00,ds00x,ds00f=dp00,dp00x,dp00f)
C ---   for sigma-z (shallow-deep) use a very small ds00
C ---    (pure sigma-z also has ds00f=dp00f and ds00x=dp00x*ds00/dp00)
C ---   for z-sigma (shallow-deep) use a very large dp00 (not recommended)
C ---   for sigma-only set nsigma=kdm, dp00 large, and ds00 small
C
C --- version 2.2.58 definition of deep and shallow z-levels.
C --- terrain following starts at depth sum(dp0k(k),k=1,nsigma) and
C --- ends at depth sum(ds0k(k),k=1,nsigma), and the depth of the k-th
C --- layer interface varies linearly with total depth between these
C --- two reference depths.
C
C --- previous to 2.2.58, it was layer thickness (not layer interface
C --- depth) that varied linearly with total depth.  These two approachs
C --- are identical for "pure sigma-z", but differ if ds00f/=dp00f.
C
        CALL BLKINR2(BLK,J,'isotop','(a6," =",f10.4," m")',
     &                     'dp00  ','(a6," =",f10.4," m")')
        IF     (J.EQ.1) THEN
          ISOTOP = BLK
          CALL BLKINR2(BLK,J,'dp0k  ','(a6," =",f10.4," m")',
     &                       'dp00  ','(a6," =",f10.4," m")')
          IF     (J.EQ.1) THEN
            DP0K(1) = BLK
            DP00    = -1.0 !no d[sp]00* input
          ELSE
            DP00    = BLK
          ENDIF
        ELSE
          DP00   = BLK
          IF     (ISOPYC) THEN
            ISOTOP = -1.0  !turn off isotop
          ELSE
            ISOTOP = 0.01  !first layer always fixed
          ENDIF
          write(6,'(a6," =",f10.4," m")') 'isotop',ISOTOP
        ENDIF
        IF     (DP00.GE.0.0) then
          CALL BLKINR(DP00X, 'dp00x ','(a6," =",f10.4," m")')
          CALL BLKINR(DP00F, 'dp00f ','(a6," =",f10.4," ")')
          CALL BLKINR(DS00,  'ds00  ','(a6," =",f10.4," m")')
          CALL BLKINR(DS00X, 'ds00x ','(a6," =",f10.4," m")')
          CALL BLKINR(DS00F, 'ds00f ','(a6," =",f10.4," ")')
        ELSE
          DP00X = DP0K(1)
          DP00F = 0.0
          DO K=2,KDM
            CALL BLKINR(DP0K(K),'dp0k  ','(a6," =",f10.4," m")')
            DP00X = MAX( DP00X, DP0K(K) )
            DP00F = DP00F + DP0K(K)/DP0K(K-1)
          ENDDO
C         dp00f: skip equal sized near bottom levels
          DO K=KDM,2,-1
            IF     (DP0K(K).NE.DP0K(K-1) .OR. K.EQ.2) THEN
              EXIT
            ENDIF
            DP00F = DP00F - 1.0
          ENDDO
          DP00F = DP00F / REAL(K-1)
          DO K=1,NSIGMA
            CALL BLKINR(DS0K(K),'ds0k  ','(a6," =",f10.4," m")')
          ENDDO
        ENDIF
      ENDIF
      CALL BLKINR(DP00I, 'dp00i ','(a6," =",f10.4," m")')
C
      IF (NSIGMA.LE.1) THEN
        NSIGMA=1
        IF     (DP00.GE.0.0) then
          DS00  =DP00
          DS00X =DP00X
          DS00F =DP00F
        ELSE
          DS0K(1)=DP0K(1)
        ENDIF
      ENDIF
C
C --- 'thflag' = reference pressure flag (0=Sigma-0, 2=Sigma-2)
      WRITE(6,*)
      CALL BLKINI(THFLAG,'thflag')
      IF     (THFLAG.EQ.0) THEN
        WRITE(SIGFMT,'(A,I2,A)')
     &    '(a6," =",f10.4," sigma-0 (sigver=',SIGVER,')")'
        IF     (MOD(SIGVER,2).NE.1 .OR. SIGVER.GE.40) THEN
          WRITE(6,*)
          WRITE(6,*) 'ERROR - thflag=0 not compatible with sigver =',
     &               SIGVER
          WRITE(6,*)
          CALL ZHFLSH(6)
          STOP
        ENDIF
      ELSEIF (THFLAG.EQ.2) THEN
        WRITE(SIGFMT,'(A,I2,A)')
     &    '(a6," =",f10.4," sigma-2 (sigver=',SIGVER,')")'
        IF     (MOD(SIGVER,2).NE.0 .OR. SIGVER.GE.40) THEN
          WRITE(6,*)
          WRITE(6,*) 'ERROR - thflag=2 not compatible with sigver =',
     &               SIGVER
          WRITE(6,*)
          CALL ZHFLSH(6)
          STOP
        ENDIF
      ELSEIF (THFLAG.EQ.4) THEN
        WRITE(SIGFMT,'(A,I2,A)')
     &    '(a6," =",f10.4," sigma-4 (sigver=',SIGVER,')")'
        IF     (SIGVER.LT.40) THEN
          WRITE(6,*)
          WRITE(6,*) 'ERROR - thflag=4 not compatible with sigver =',
     &               SIGVER
          WRITE(6,*)
          CALL ZHFLSH(6)
          STOP
        ENDIF
      ELSE
        WRITE(6,*)
        WRITE(6,*) 'ERROR - thflag must be 0 or 2 or 4'
        WRITE(6,*)
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
C --- 'thbase' = reference density (sigma units)
      WRITE(6,*)
      CALL BLKINR(THBASE,'thbase',SIGFMT)
C
C --- 'vsigma' = spacially varying isopycnal layer target densities (0=F,1=T)
      WRITE(6,*)
      CALL BLKINL(VSIGMA,'vsigma')
C --- 'sigma ' = isopycnal layer target densities (sigma units)
      DO K=1,KDM
        CALL BLKINR(SIGMA(K),'sigma ',SIGFMT)
        IF     (K.GT.1) THEN
          IF     (SIGMA(K).LT.SIGMA(K-1)) THEN
            WRITE(6,*)
            WRITE(6,*) 'ERROR - sigma(k) must be > sigma(k-1)'
            WRITE(6,*)
            CALL ZHFLSH(6)
            STOP
          ENDIF
        ENDIF
      ENDDO
C
C --- 'thkmin' = minimum mixed-layer thickness (m)
      WRITE(6,*)
      CALL BLKINR(THKMIN,'thkmin','(a6," =",f10.4," m")')
      WRITE(6,*)
      CLOSE(UNIT=99)
C
C     CALCULATE DP0K AND DS0K
C
      IF     (DP00.LT.0.0) then
C
C       ALREADY INPUT
C
        DPMS=0.0
        DPCK(0)=0.0
C
        DO K=1,KDM
          DPMS=DPMS+DP0K(K)
          DPCK(K)=DPMS
          IF     (THFLAG.EQ.0) THEN
            WRITE(6,6000) K,DP0K(K),DPMS,SIGMA(K),' Sig-th'
          ELSEIF (THFLAG.EQ.2) THEN
            WRITE(6,6000) K,DP0K(K),DPMS,SIGMA(K),' Sigma2'
          ENDIF
        ENDDO
        WRITE(6,*)
        CALL ZHFLSH(6)
        DSMS=0.0
        DSCK(0)=0.0
        DO K=1,NSIGMA
          DSMS=DSMS+DS0K(K)
          DSCK(K)=DSMS
        ENDDO
        DO K= NSIGMA+1,KDM
          DS0K(K)=0.0
          DSCK(K)=DSMS
        ENDDO
      ELSE
C
C       MINIMUM (DEEP) LAYER THICKNESSES.
C
        IF     (ISOPYC) THEN
          DP0K(1)=THKMIN
        ELSE
          DP0K(1)=DP00
        ENDIF
        DPMS=DP0K(1)
        DPCK(0)=0.0
        DPCK(1)=DPMS
        IF     (THFLAG.EQ.0) THEN
          WRITE(6,6000) 1,DP0K(1),DPMS,SIGMA(1),' Sig-th'
        ELSEIF (THFLAG.EQ.2) THEN
          WRITE(6,6000) 1,DP0K(1),DPMS,SIGMA(1),' Sigma2'
        ENDIF
 6000   FORMAT(       'k =',I3,
     +         '   thkns =',F6.1,' m',
     +         '   depth =',F8.1,' m',
     +         '   density =',F7.3,A)
C
        DP0KF=ONE
        DO K=2,KDM
          DP0KF=DP0KF*DP00F
          IF     (K.LE.NHYBRD) THEN
            DP0K(K)=MIN(DP00*DP0KF,DP00X)
          ELSE
            DP0K(K)=ZERO
          ENDIF
          DPMS=DPMS+DP0K(K)
          DPCK(K)=DPMS
          IF     (THFLAG.EQ.0) THEN
            WRITE(6,6000) K,DP0K(K),DPMS,SIGMA(K),' Sig-th'
          ELSEIF (THFLAG.EQ.2) THEN
            WRITE(6,6000) K,DP0K(K),DPMS,SIGMA(K),' Sigma2'
          ENDIF
        ENDDO
        WRITE(6,*)
        CALL ZHFLSH(6)
C
C       MINIMUM (SHALLOW) LAYER THICKNESSES.
C
        IF     (ISOPYC) THEN
          DS0K(1)=THKMIN
        ELSE
          DS0K(1)=DS00
        ENDIF
        DSMS=DS0K(1)
        DSCK(0)=0.0
        DSCK(1)=DSMS
        DS0KF=ONE
        DO K=2,NSIGMA
          DS0KF=DS0KF*DS00F
          DS0K(K)=MIN(DS00*DS0KF,DS00X)
          DSMS=DSMS+DS0K(K)
          DSCK(K)=DSMS
        ENDDO
        DO K= NSIGMA+1,KDM
          DS0K(K)=0.0
          DSCK(K)=DSMS
        ENDDO
      ENDIF !DP0K,DS0K
C
C     TOPOGRAPHY INPUT.
C
      CALL ZAIOST
C
      CALL ZHOPNC(51,"regional.depth.b",'FORMATTED', 'OLD', 0)
      READ (51,'(A79)') PREAMBL
      READ (51,'(A)')   CLINE
      CLOSE(UNIT=51)
      WRITE(6,'(/(1X,A79))') PREAMBL,CLINE
C
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
C
      CALL ZAIOPF('regional.depth.a','OLD', 51)
      CALL ZAIORD(DEPTH,MSK,.FALSE., HMINA,HMAXA, 51)
      CALL ZAIOCL(51)
C
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b topography files not consistent:',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
C     INITIALIZE LAND MASK.
C
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (ABS(DEPTH(I,J)).LT.2.0**90.0) THEN
            MSK(  I,J) = 1
          ELSE
            MSK(  I,J) = 0
            DEPTH(I,J) = ZERO
          ENDIF
        ENDDO
      ENDDO
C
C     TARGET DENSITIES.
C

      ALLOCATE( SIG3D(IDM,JDM,KDM) )
      ALLOCATE( SIG3D_TMP(KDM) )

C
      IF     (.NOT.VSIGMA) THEN
        DO K= 1,KDM
          SIG3D(:,:,K) = SIGMA(K)
        ENDDO
      ELSE
        CALL ZAIOPN('OLD', 52)
        DO K= 1,KDM
          CALL ZAIORD(SIG3D(1,1,K),MSK,.FALSE., HMINA,HMAXA, 52)
          IF     (HMINA.GT.SIGMA(K)+0.005 .OR.
     &            HMAXA.LT.SIGMA(K)-0.005     ) THEN
            WRITE(6,'(/ a,i3,a /)')
     &        'ERROR - VARIABLE TARGET DENSITY FOR LAYER',
     &        K,' IS NOT CONSISTENT WITH SIGMA(K)'
            WRITE(6,*) 'SIGMA(K)    = ',SIGMA(K)
            WRITE(6,*) 'HMINA,HMAXA = ',HMINA,HMAXA
            STOP
          ENDIF
        ENDDO
        CALL ZAIOCL(52)
      ENDIF
C
C     CHECK ITEST,JTEST.
C
      IF     (MIN(ITEST,JTEST).GT.0) THEN
        IF     (ITEST.GT.IDM) THEN
          WRITE(6,'(/ a /)') 'error - itest > idm'
          CALL ZHFLSH(6)
          STOP
        ELSEIF (JTEST.GT.JDM) THEN
          WRITE(6,'(/ a /)') 'error - jtest > jdm'
          CALL ZHFLSH(6)
          STOP
        ELSEIF(MSK(ITEST,JTEST).EQ.0) THEN
          WRITE(6,'(/ a /)') 'error - itest,jtest is a land point'
          CALL ZHFLSH(6)
          STOP
        ENDIF
      ENDIF
C
C     LATITUDE GRID INPUT.
C
      CALL ZHOPNC(31, 'regional.grid.b', 'FORMATTED', 'OLD', 0)
      CALL ZAIOPF('regional.grid.a', 'OLD', 31)
C
      READ(31,*) ! skip idm
      READ(31,*) ! skip jdm
      READ(31,*) ! skip mapflg
      READ(31,*) ! skip plon
      CALL ZAIOSK(31)
      READ(31,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
      CALL ZAIORD(PLAT,MSK,.FALSE., HMINA,HMAXA, 31)
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (plat):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF

      CLOSE(UNIT=72)
C
C     Z-LEVEL CLIMATOLOGY INPUT.
C

C
C THIS SECTION IS A FORK OF HYCOM_IO ROUTINES USED IN HYC2PROJ PROGRAM.
C I HAVE REMOVED UNNECESSARY PARTS FROM THE ORIGINAL CODES, MERGED IT
C WITH KNUT'S (KAL) ADDITIONS TO THIS SCRIPT. THE DEFAULT SCRIPT THAT WE USE
C HERE REQUIRES SEPARATE FILES OF TEMP, SAL ETC... BUT WE HAVE A SINGLE
C OUTPUT FROM NEMO GLOBAL DOMAIN WITH MANY VARIABLES. FOR FUTURE PURPOSES
C AS WELL, I DECIDED TO MODIFY THE CODE ACCORDINGLY.
C


      CALL READ_DEPTH(KZ,ZL,flnm_z)
C
c     DOES NOT WORK
C     SIMILAR PRINCIPAL BELOW, SO LEFT IT MALFUNCTIONING
C     MAYBE IN THE FUTURE, THIS CAN BE FIXED
C
C     /////////////// 2D - SURFACE FIELDS ///////////
      ALLOCATE(      MTG1(IDM,JDM))
      ALLOCATE(      SRFH(IDM,JDM))
      ALLOCATE(      SUFLX(IDM,JDM))
      ALLOCATE(      SAFLX(IDM,JDM))
      ALLOCATE(      BLDP(IDM,JDM))
      ALLOCATE(      MIXDP(IDM,JDM))
      ALLOCATE(      UBRTP(IDM,JDM))
      ALLOCATE(      VBRTP(IDM,JDM))
C /////////////////// NOW 3D /////////////////////////////

      ALLOCATE(      TZ1(IDM,JDM) )
      ALLOCATE(      SZ1(IDM,JDM) )
      ALLOCATE(      TH1(IDM,JDM) ) ! thickness
      ALLOCATE(      UV1(IDM,JDM) )
      ALLOCATE(      VV1(IDM,JDM) )

      ALLOCATE(      TH(KZ+1,IDM,JDM) ) ! knut assigned kz+1 to other
      ALLOCATE(      UV(KZ+1,IDM,JDM) )
      ALLOCATE(      VV(KZ+1,IDM,JDM) )


      CALL READ_ARCHIVE(KZ,IDM,JDM,NREC,KDM,COORD,
     +     tlevel1,ARCHVFILE,cfld,MTG1,SRFH,SUFLX,SAFLX,
     +     BLDP,MIXDP,UBRTP,VBRTP,TZ,SZ,UV,VV,RZ,ZL,
     +     flag_t,flag_s,flag_th,flag_u,flag_v,ISOPYC,
     +     SIGVER,LEVTOP,SIG3D,DEPTH,TZ1,SZ1,UV1,VV1,TH,TH1,
     +     flag_no3,flag_po4,flag_si,NO3,PO4,SI)
C
C    UV1 and VV1 are on P-cell and are converted into the u/v cells after
C    subroutine GRIRD_REMAPPING
C


C
C     DIAGNOSTIC PRINTOUT.
C
      IF     (MIN(ITEST,JTEST).GT.0) THEN
        WRITE(6,*)
        DO K= 1,KZ
          WRITE(6,'(A,2I5,I3,A,F9.4,A,3F7.3,A,2F7.3)')
     +     'I,J,K =',ITEST,JTEST,K,
     +     '  ZLEV =',ZL(K),
     +     '   R,T,S =',RZ(K,ITEST,JTEST),
     +                  TZ(K,ITEST,JTEST),
     +                  SZ(K,ITEST,JTEST),
     +     '  U,V =',UV(K,ITEST,JTEST),
     +                VV(K,ITEST,JTEST)
        ENDDO
        CALL ZHFLSH(6)
      ENDIF

      PREAMBL(1) = CTITLE
      IF     (NSIGMA.EQ.1) THEN
        WRITE(PREAMBL(2),4000) IEXPT/10,MOD(IEXPT,10),
     +                         NHYBRD,0,
     +                         DP00,DP00,DP00X,DP00F
      ELSE
        WRITE(PREAMBL(2),4000) IEXPT/10,MOD(IEXPT,10),
     +                         NHYBRD,NSIGMA,
     +                         DS00,DP00,DP00X,DP00F
      ENDIF
      MONTH =0
      TIME  =0


c ---   output is in "*.[AB]"
      open (unit=21,file=trim(archvfileo)//".b",form='formatted',
     &          status='new',action='write')
      call zaiopf(trim(archvfileo)//".a",'new', 21)

c ---   bio output is in "*.[AB]"
c ---   if you are using 'archm' for nesting use the following
cc      dummy_index=index(archvfileo,'archm')
c ---
      if (flag_no3.ne."NONE" .or.flag_po4.ne."NONE"
     &    .or.flag_si.ne."NONE" ) then
         dummy_index=index(archvfileo,'archv')
         dummy_name=archvfileo(dummy_index:dummy_index+4)//"_fabm"
     +      //archvfileo(dummy_index+5:dummy_index+17)
         open (unit=211,file=trim(archvfileo(1:dummy_index-1))//
     &     trim(dummy_name)//".b",form='formatted',
     &          status='new',action='write')
         call zaiopf(trim(archvfileo(1:dummy_index-1))//
     &     trim(dummy_name)//".a",'new', 211)
         WRITE(211,4200) MONTH,PREAMBL(1),PREAMBL(2),
     &               IVERSN,IEXPT,YRFLAG,IDM,JDM
      endif
      TIME=MODEL_DAY

      YRFLAG = MAX(0,MIN(3,YRFLAG))
      WRITE(21,4200) MONTH,PREAMBL(1),PREAMBL(2),
     &               IVERSN,IEXPT,YRFLAG,IDM,JDM

      IF     (YRFLAG.EQ.0) THEN  ! 360 days, starting Jan 16
        TIME = (MONTH-1)*30.0
      ELSEIF (YRFLAG.EQ.1) THEN  ! 366 days, starting Jan 16
        TIME = (MONTH-1)*30.5
      ELSEIF (YRFLAG.EQ.2) THEN  ! 366 days, starting Jan 1
        TIME = (MONTH-1)*30.5 + 15.0
      ENDIF


C
C THIS IS IMPORTANT TO EXPLAIN
C ORIGINAL CODE CALCULATES A MIXED LAYER DEPTH HERE BASED
C ON T=0.5C DIFFERENCE FROM THE SURFACE AND CALCULATES 
C DENSITY FOR THE MLD AND TRIES TO LOCATE IT WITHIN OUR 
C Z-LEV DENSITIES. MOST PROBABLY THE RESULTING MLD WILL 
C NOT BE USED AS UPPER LAYERS ARE NOT ISOPYCNAL
C 

      DO J= 1,JDM
        DO I= 1,IDM
          IF (DEPTH(I,J).GT.ZERO .AND. DEPTH(I,J) .LT. 2.0**99.0) THEN
            PCM0(I,J) = ZERO
            PCM1(I,J) = ZERO
            PCM2(I,J) = ZERO
            PKM2(I,J) = ZERO
            PKM1(I,J) = ZERO
            TM(  I,J) = TZ(1,I,J)-0.25
            SM(  I,J) = SZ(1,I,J)
            SIGMAA    = SIG_V(TZ(1,I,J)-0.5,SM(I,J),SIGVER)
            CALL FIND_DENSITY(RZLOC,SIGMAA,RZ(1,I,J),KZ+1,1)
            IF     (RZLOC.EQ.0.0) THEN
              PMIX(I,J) = ZL(2)
            ELSEIF (RZLOC.EQ.KZ+1) THEN
              PMIX(I,J) = DEPTH(I,J)
            ELSE
              L = RZLOC
              Q = RZLOC - L
              PMIX(I,J) = (1.0-Q)*ZL(L) + Q*ZL(L+1)
            ENDIF
            PMIX(I,J) = MAX( THKMIN, MIN( PMIX(I,J), DEPTH(I,J) ) )
              if (ldebug.and.i.eq.itest .and. j.eq.jtest) then
                WRITE(6,'(A,2F10.3)')
     +            'PMIX,RZLOC =',PMIX(I,J),RZLOC
                call zhflsh(6)
              endif !debug
          ELSE
            PCM0(I,J) = ZERO
            PCM1(I,J) = ZERO
            PCM2(I,J) = ZERO
            PKM2(I,J) = ZERO
            PKM1(I,J) = ZERO
            PMIX(I,J) = THKMIN
            TM(  I,J) = ZERO
            SM(  I,J) = ZERO
            RM(  I,J) = ZERO
C
            PM(  I,J) = ZERO
C
            DEPTH(I,J) = ZERO  ! should already be zero
          ENDIF
        ENDDO
      ENDDO
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (DEPTH(I,J).GT.ZERO) THEN
            WORK(I,J) = ZERO
          ENDIF
        ENDDO
      ENDDO

!      CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)

      CALL ZAIOWR(MTG1,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'montg1  ',MONTH,TIME,0,ZERO,XMIN,XMAX
      CALL ZAIOWR(SRFH,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'srfhgt  ',MONTH,TIME,0,ZERO,XMIN,XMAX
      CALL ZAIOWR(SUFLX,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'surflx  ',MONTH,TIME,0,ZERO,XMIN,XMAX
      CALL ZAIOWR(SAFLX,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'salflx  ',MONTH,TIME,0,ZERO,XMIN,XMAX
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (DEPTH(I,J).GT.ZERO) THEN
            WORK(I,J) = PMIX(I,J)*9806.0
          ENDIF
        ENDDO
      ENDDO
      CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'bl_dpth ',MONTH,TIME,0,ZERO,XMIN,XMAX
      CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'mix_dpth',MONTH,TIME,0,ZERO,XMIN,XMAX
      DO J= 1,JDM
        DO I= 1,IDM
          IF(DEPTH(I,J).GT.ZERO .AND.TZ(1,I,J) .LT. 1000.0 ) THEN
            TM(  I,J) = TZ(1,I,J)-0.25
            SM(  I,J) = SZ(1,I,J)
            WORK(I,J) = SIG_V(TM(I,J),SM(I,J),SIGVER) - THBASE
          ELSE
            TM(  I,J) = ZERO
            SM(  I,J) = ZERO
            WORK(I,J) = ZERO
          ENDIF
        ENDDO
      ENDDO
      CALL ZAIOWR(TM,  MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'tmix    ',MONTH,TIME,0,  ZERO,XMIN,XMAX
      CALL ZAIOWR(SM,  MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'smix    ',MONTH,TIME,0,  ZERO,XMIN,XMAX
      CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'thmix   ',MONTH,TIME,0,THBASE,XMIN,XMAX
      DO J= 1,JDM
        DO I= 1,IDM
          WORK(I,J) = ZERO
        ENDDO
      ENDDO
      CALL ZAIOWR(WORK,MSK,.FALSE., XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'umix    ',MONTH,TIME,0,ZERO,XMIN,XMAX
      CALL ZAIOWR(WORK,MSK,.FALSE., XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'vmix    ',MONTH,TIME,0,ZERO,XMIN,XMAX
      CALL ZAIOWR(UBRTP,MSK,.FALSE., XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'u_btrop ',MONTH,TIME,0,ZERO,XMIN,XMAX
      CALL ZAIOWR(VBRTP,MSK,.FALSE., XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'v_btrop ',MONTH,TIME,0,ZERO,XMIN,XMAX
      
C
C     NOW THE WATER COLUMN
C
C     PROCESS ALL THE LAYERS.
C
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (DEPTH(I,J).GT.ZERO) THEN
             RKM1(I,J) = ZERO
             PCM0(I,J) = ZERO
             PCM1(I,J) = ZERO
             PCM2(I,J) = ZERO
             PKM1(I,J) = ZERO
             PKM2(I,J) = ZERO
             PM(I,J) = ZERO
            PZTOP(I,J) = 1.0  !surface
            PZBOT(I,J) = 1.0  !surface
          ENDIF
        ENDDO
      ENDDO


      DO 810 K= 1,KDM
        DO J= 1,JDM
          DO I= 1,IDM


            IF     ((DEPTH(I,J)).GT.ZERO) THEN

            SIG3D_TMP=SIG3D(I,J,:)
            call GRIRD_REMAPPING(i,j,K,KDM,KZ,NSIGMA,SIG3D_TMP,RM(I,J),
     +          PM(I,J),RZ(:,I,J),DEPTH(I,J),DSCK,DPCK,ISOTOP,
     +          ISOPYC,PKM1(I,J),PZBOT(I,J),PZTOP(I,J),RKM1(I,J),
     +          ZL,PMIX(I,J),PCM0(I,J),PCM1(I,J),PCM2(I,J),
     +          PKM2(I,J),DP0K,DP00I,DS0K,flag_t,flag_s,
     +          flag_u,flag_v,itest,jtest,SM(I,J),TM(I,J),
     +          PU(I,J),PV(I,J),UV(:,I,J),VV(:,I,J),
     +          ldebug,sigver,TZ(:,I,J),SZ(:,I,J),
     +          flag_no3,flag_po4,flag_si,NO3(:,I,J),PO4(:,I,J),
     +          SI(:,I,J),NO3M(I,J),PO4M(I,J),SIM(I,J))

            ENDIF  !DEPTH>0

          ENDDO  !J=1,JDM
        ENDDO  !I=1,IDM
C ---   CONVERT U AND V FROM P-CELL TO U- AND V-CELLS
        CALL P2UV(IDM,JDM,PU,PV,UM,VM)


C
C   Write bio [ab] files
C
        if     (flag_no3.ne."NONE") then
           NO3M(:,2)=NO3M(:,3)
           CALL ZAIOWR(NO3M,  MSK,.TRUE.,  XMIN,XMAX, 211, .FALSE.)
           WRITE(211,4201) 'ECO_no3    ',MONTH,TIME,K,SIGMA(K),XMIN,XMAX
        endif
        if     (flag_po4.ne."NONE") then
           PO4M(:,2)=PO4M(:,3)
           CALL ZAIOWR(PO4M,  MSK,.TRUE.,  XMIN,XMAX, 211, .FALSE.)
           WRITE(211,4201) 'ECO_pho    ',MONTH,TIME,K,SIGMA(K),XMIN,XMAX
        endif
        if     (flag_si.ne."NONE") then
           SIM(:,2)=SIM(:,3)
           CALL ZAIOWR(SIM,  MSK,.TRUE.,  XMIN,XMAX, 211, .FALSE.)
           WRITE(211,4201) 'ECO_sil    ',MONTH,TIME,K,SIGMA(K),XMIN,XMAX
        endif

C
C       WRITE OUT DUMMY ARCHIVE.
C
        DO J= 1,JDM
          DO I= 1,IDM
            WORK(I,J) = ZERO
          ENDDO
        ENDDO
        if     (flag_u.ne."NONE") then
            CALL ZAIOWR(UM,MSK,.FALSE., XMIN,XMAX, 21, .FALSE.)
            WRITE(21,4201) 'u-vel.  ',MONTH,TIME,K,SIGMA(K),XMIN,XMAX
        endif
        if     (flag_v.ne."NONE") then
            CALL ZAIOWR(VM,MSK,.FALSE., XMIN,XMAX, 21, .FALSE.)
            WRITE(21,4201) 'v-vel.  ',MONTH,TIME,K,SIGMA(K),XMIN,XMAX
        endif
        if     (flag_th.ne."NONE") then

            DO J= 1,JDM
                DO I= 1,IDM
                    WORK(I,J) = (PM(I,J) - PKM1(I,J))*9806.0
                ENDDO
            ENDDO
            CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
            WRITE(21,4201) 'thknss  ',MONTH,TIME,K,SIGMA(K),XMIN,XMAX
        endif
        if     (flag_t.ne."NONE") then

            CALL ZAIOWR(TM,  MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
            WRITE(21,4201) 'temp    ',MONTH,TIME,K,SIGMA(K),XMIN,XMAX
        endif

        if     (flag_s.ne."NONE") then

            CALL ZAIOWR(SM,  MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
            WRITE(21,4201) 'salin   ',MONTH,TIME,K,SIGMA(K),XMIN,XMAX
        endif
        DO J= 1,JDM
          DO I= 1,IDM
            WORK(I,J) = RM(I,J) - THBASE
          ENDDO
        ENDDO
        CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
        WRITE(21,4201) 'density ',MONTH,TIME,K,SIGMA(K),XMIN,XMAX

C
C       PREPARE FOR NEXT K-LOOP
C
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (DEPTH(I,J).GT.ZERO) THEN
               PCM2(I,J) =  PCM1(I,J)
               PCM1(I,J) =  PCM0(I,J)
               PKM2(I,J) =  PKM1(I,J)
               PKM1(I,J) =  PM(  I,J)
               RKM1(I,J) =  RM(  I,J)
              PZTOP(I,J) = PZBOT(I,J)
              PZBOT(I,J) = 0.0
            ENDIF
          ENDDO
        ENDDO
C
 810  CONTINUE  !K=1,KDM
C

      CALL ZAIOCL(21)

      CALL ZAIOCL(211)
      stop
C
C
C


C
 4000 FORMAT('Expt ',I2.2,'.',I1.1,
     +       '  nhybrd=',I2,
     +        ' nsigma=',I2,
     +          ' ds00=',F5.2,
     +          ' dp00=',F5.2,
     +         ' dp00x=',F6.1,
     +         ' dp00f=',F5.3)
 4101 FORMAT(A79)
 4102 FORMAT(A,': month,layer,dens,range = ',I2.2,I4.2,F7.3,1P2E16.7)
 4200 FORMAT(
     + 'Dummy HYCOM archive from climatology for month ',I2.2,'.' /
     + A80/A80/
     + '1234567890123456789012345678901234567890',
     + '1234567890123456789012345678901234567890'/
     & I5,4X,'''iversn'' = hycom version number x10'/
     & I5,4X,'''iexpt '' = experiment number x10'/
     & I5,4x,'''yrflag'' = days in year flag'/
     & I5,4x,'''idm   '' = longitudinal array size'/
     & I5,4x,'''jdm   '' = latitudinal  array size'/
     & 'field       time step  model day',
     & '  k  dens        min              max')
 4201 FORMAT(a8,' =',i11,f11.2,i3,f7.3,1p2e16.7)
 5000 FORMAT(A40)
 5500 FORMAT(6E13.6)
 6200 FORMAT(/ 1X,'MIN,MAX I COORDS = ',F8.2,',',F8.2 
     +       / 1X,'MIN,MAX J COORDS = ',F8.2,',',F8.2 /)
 6300 FORMAT(10X,'WRITING CLIM LAYER',I5,'    SIGMA =',F7.3 /)
 8100 FORMAT(1X,A,': min=',F9.2,' ave=',F9.2,' max=',F9.2,
     +   '   (k,sigma =',i3,F7.2,')')
 8200 FORMAT(' k =',I3.2,' j = ',I4.4,' to ',I4.4,
     +       ' sig =',F7.3,' lat =',F6.1,' inf =',F8.2)
C     END OF PROGRAM WNDINT.
      END














      SUBROUTINE LAYSTAT(PM, THICK, IDM,JDM, PMIN,PAVE,PMAX)
      IMPLICIT NONE
C
      INTEGER IDM,JDM
      REAL*4  PM(IDM,JDM),THICK(IDM,JDM), PMIN,PAVE,PMAX
C
C --- CALCULATE STATISTICS FOR PM.
C --- ONLY WHERE LAYER THICKNESS IS AT LEAST 10 CM.
C --- AVERAGE DOES NOT ALLOW FOR VARIATION IN GRID CELL SIZE.
C
      REAL*4     TENCM
      PARAMETER (TENCM=0.1)
C
      INTEGER I,J,IJSUM
      REAL*4  PPMIN,PPMAX
      REAL*8  PPSUM
C
      PPMIN =  1.E10
      PPMAX = -1.E10
      IJSUM =  0
      PPSUM =  0.0D0
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (THICK(I,J).GT.TENCM) THEN
            PPMIN = MIN( PPMIN, PM(I,J) )
            PPMAX = MAX( PPMAX, PM(I,J) )
            IJSUM = IJSUM + 1
            PPSUM = PPSUM + PM(I,J)
          ENDIF
        ENDDO
      ENDDO
      IF     (IJSUM.NE.0) THEN
        PMIN = PPMIN
        PMAX = PPMAX
        PAVE = PPSUM / IJSUM
      ELSE
        PMIN = 99.9
        PMAX = 99.9
        PAVE = 99.9
      ENDIF
      RETURN
C     END OF LAYSTAT.
      END
      SUBROUTINE FIND_DENSITY(DENLOC,DENTARG,RZ,KZ,MINZ)
      IMPLICIT NONE
C
      INTEGER KZ,MINZ
      REAL*4  DENLOC,DENTARG,RZ(KZ)
C
C     FIND EXACT LOCATION IN LAYER SPACE OF DENTARG
C     SEARCHING WITHIN RZ(MINZ:KZ)
C
C     ASSUME RZ IS MONOTONICALLY NON-DECREASING.
C
C     RETURN 0.0 IF DENTARG < RZ(MINZ)
C     RETURN KZ  IF DENTARG > RZ(KZ)
C
      INTEGER K
C

      IF     (DENTARG.LT.RZ(MINZ)) THEN
        DENLOC = 0.0
      ELSEIF (DENTARG.EQ.RZ(MINZ)) THEN
        DENLOC = MINZ
      ELSEIF (DENTARG.GE.RZ(KZ)) THEN
        DENLOC = KZ
      ELSE
        DO K= MINZ+1,KZ
          IF     (DENTARG.LE.RZ(K)) THEN !dentarg>rz(k-1)
            DENLOC = K-1 + (DENTARG-RZ(K-1))/(RZ(K)-RZ(K-1))
            EXIT
          ENDIF
        ENDDO !k
      ENDIF
      RETURN
      END
      SUBROUTINE FIND_DEPTH(ZLOC,ZTARG,Z,KZ,MINZ)
      IMPLICIT NONE
C
      INTEGER KZ,MINZ
      REAL*4  ZLOC,ZTARG,Z(KZ)
C
C     FIND EXACT LOCATION IN LAYER SPACE OF ZTARG
C     SEARCHING WITHIN Z(MINZ:KZ)
C
C     RETURN 0.0 IF ZTARG < Z(MINZ)
C     RETURN KZ  IF ZTARG > Z(KZ)
C
      INTEGER K
C
      IF     (ZTARG.LT.Z(MINZ)) THEN
        ZLOC = 0.0
      ELSEIF (ZTARG.EQ.Z(MINZ)) THEN
        ZLOC = MINZ
      ELSEIF (ZTARG.GE.Z(KZ)) THEN
        ZLOC = KZ
      ELSE
        DO K= MINZ+1,KZ
          IF     (ZTARG.LE.Z(K)) THEN !ztarg>z(k-1)
            ZLOC = K-1 + (ZTARG-Z(K-1))/(Z(K)-Z(K-1))
            EXIT
          ENDIF
        ENDDO !k
      ENDIF
      RETURN
      END
      REAL*4 FUNCTION SIG_V(TT,SS,SIGVER)
      IMPLICIT NONE
      INTEGER SIGVER
      REAL*4  TT,SS
C
C     SIGVER WRAPPER FOR SIG
C
      REAL*8 SS8,TT8
      REAL*8 SIG_1,SIG_2,SIG_3,SIG_4,SIG_5,SIG_6,SIG_7,SIG_8,
     &       SIG_46,SIG_48
C
      TT8 = TT
      SS8 = SS
      IF     (SIGVER.GT.40) THEN
        IF     (SIGVER.EQ.46) THEN
          SIG_V = SIG_46(TT8,SS8)
        ELSEIF (SIGVER.EQ.48) THEN
          SIG_V = SIG_48(TT8,SS8)
        ENDIF
      ELSEIF (MOD(SIGVER,2).EQ.1) THEN
        IF     (SIGVER.EQ.1) THEN
          SIG_V = SIG_1(TT8,SS8)
        ELSEIF (SIGVER.EQ.3) THEN
          SIG_V = SIG_3(TT8,SS8)
        ELSEIF (SIGVER.EQ.5) THEN
          SIG_V = SIG_5(TT8,SS8)
        ELSEIF (SIGVER.EQ.7) THEN
          SIG_V = SIG_7(TT8,SS8)
        ENDIF
      ELSE
        IF     (SIGVER.EQ.2) THEN
          SIG_V = SIG_2(TT8,SS8)
        ELSEIF (SIGVER.EQ.4) THEN
          SIG_V = SIG_4(TT8,SS8)
        ELSEIF (SIGVER.EQ.6) THEN
          SIG_V = SIG_6(TT8,SS8)
        ELSEIF (SIGVER.EQ.8) THEN
          SIG_V = SIG_8(TT8,SS8)
        ENDIF
      ENDIF
      RETURN
      END
      REAL*4 FUNCTION SOFSIG_V(RR,TT,SIGVER)
      IMPLICIT NONE
      INTEGER SIGVER
      REAL*4  RR,TT
C
C     SIGVER WRAPPER FOR SOFSIG
C
      REAL*8 RR8,TT8
      REAL*8 SOFSIG_1,SOFSIG_2,SOFSIG_3,SOFSIG_4,
     &       SOFSIG_5,SOFSIG_6,SOFSIG_7,SOFSIG_8,
     &       SOFSIG_46,SOFSIG_48
C
      RR8 = RR
      TT8 = TT
      IF     (SIGVER.GT.40) THEN
        IF     (SIGVER.EQ.46) THEN
          SOFSIG_V = SOFSIG_46(RR8,TT8)
        ELSEIF (SIGVER.EQ.48) THEN
          SOFSIG_V = SOFSIG_48(RR8,TT8)
        ENDIF
      ELSEIF (MOD(SIGVER,2).EQ.1) THEN
        IF     (SIGVER.EQ.1) THEN
          SOFSIG_V = SOFSIG_1(RR8,TT8)
        ELSEIF (SIGVER.EQ.3) THEN
          SOFSIG_V = SOFSIG_3(RR8,TT8)
        ELSEIF (SIGVER.EQ.5) THEN
          SOFSIG_V = SOFSIG_5(RR8,TT8)
        ELSEIF (SIGVER.EQ.7) THEN
          SOFSIG_V = SOFSIG_7(RR8,TT8)
        ENDIF
      ELSE
        IF     (SIGVER.EQ.2) THEN
          SOFSIG_V = SOFSIG_2(RR8,TT8)
        ELSEIF (SIGVER.EQ.4) THEN
          SOFSIG_V = SOFSIG_4(RR8,TT8)
        ELSEIF (SIGVER.EQ.6) THEN
          SOFSIG_V = SOFSIG_6(RR8,TT8)
        ELSEIF (SIGVER.EQ.8) THEN
          SOFSIG_V = SOFSIG_8(RR8,TT8)
        ENDIF
      ENDIF
      RETURN
      END
      REAL*4 FUNCTION TOFSIG_V(RR,SS,SIGVER)
      IMPLICIT NONE
      INTEGER SIGVER
      REAL*4  RR,SS
C
C     SIGVER WRAPPER FOR TOFSIG
C
      REAL*8 RR8,SS8
      REAL*8 TOFSIG_1,TOFSIG_2,TOFSIG_3,TOFSIG_4,
     &       TOFSIG_5,TOFSIG_6,TOFSIG_7,TOFSIG_8,
     &       TOFSIG_46,TOFSIG_48
C
      RR8 = RR
      SS8 = SS
      IF     (SIGVER.GT.40) THEN
        IF     (SIGVER.EQ.46) THEN
          TOFSIG_V = TOFSIG_46(RR8,SS8)
        ELSEIF (SIGVER.EQ.48) THEN
          TOFSIG_V = TOFSIG_48(RR8,SS8)
        ENDIF
      ELSEIF (MOD(SIGVER,2).EQ.1) THEN
        IF     (SIGVER.EQ.1) THEN
          TOFSIG_V = TOFSIG_1(RR8,SS8)
        ELSEIF (SIGVER.EQ.3) THEN
          TOFSIG_V = TOFSIG_3(RR8,SS8)
        ELSEIF (SIGVER.EQ.5) THEN
          TOFSIG_V = TOFSIG_5(RR8,SS8)
        ELSEIF (SIGVER.EQ.7) THEN
          TOFSIG_V = TOFSIG_7(RR8,SS8)
        ENDIF
      ELSE
        IF     (SIGVER.EQ.2) THEN
          TOFSIG_V = TOFSIG_2(RR8,SS8)
        ELSEIF (SIGVER.EQ.4) THEN
          TOFSIG_V = TOFSIG_4(RR8,SS8)
        ELSEIF (SIGVER.EQ.6) THEN
          TOFSIG_V = TOFSIG_6(RR8,SS8)
        ELSEIF (SIGVER.EQ.8) THEN
          TOFSIG_V = TOFSIG_8(RR8,SS8)
        ENDIF
      ENDIF
      RETURN
      END
      REAL*8 FUNCTION SIG_1(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_7term.h'
      SIG_1 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_1(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      INCLUDE '../../include/stmt_fns_SIGMA0_7term.h'
      SOFSIG_1 = SOFSIG(RR8,TT8)
      END
      REAL*8 FUNCTION TOFSIG_1(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_7term.h'
      TOFSIG_1 = TOFSIG(RR8,SS8)
      END
      REAL*8 FUNCTION SIG_3(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_9term.h'
      SIG_3 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_3(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      INCLUDE '../../include/stmt_fns_SIGMA0_9term.h'
      SOFSIG_3 = SOFSIG(RR8,TT8)
      END
      REAL*8 FUNCTION TOFSIG_3(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_9term.h'
      TOFSIG_3 = TOFSIG(RR8,SS8)
      END
      REAL*8 FUNCTION SIG_5(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_17term.h'
      SIG_5 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_5(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      REAL*8, PARAMETER :: TOL=1.D-6
      INTEGER NN
      REAL*8  SN,SO
      REAL*8  SOFSIG_7
      INCLUDE '../../include/stmt_fns_SIGMA0_17term.h'
C     sofsig via Newton iteration from a 12-term 1st guess
      SN = SOFSIG_7(RR8,TT8)  !non-negative
      DO NN= 1,10
        SO = SN
        SN = SO - (SIG(TT8,SO)-RR8)/DSIGDS(TT8,SO)
        IF     (NN.EQ.10 .OR. ABS(SN-SO).LT.TOL) THEN
          EXIT
        ENDIF
      ENDDO !nn
      SOFSIG_5 = SN
      END
      REAL*8 FUNCTION TOFSIG_5(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      REAL*8, PARAMETER :: TOL=1.D-6
      INTEGER NN
      REAL*8  TN,TO
      REAL*8  TOFSIG_7
      INCLUDE '../../include/stmt_fns_SIGMA0_17term.h'
C     sofsig via Newton iteration from a 12-term 1st guess
      TN = TOFSIG_7(RR8,SS8)  !non-negative
      DO NN= 1,10
        TO = TN
        TN = TO - (SIG(TO,SS8)-RR8)/DSIGDT(TO,SS8)
        IF     (NN.EQ.10 .OR. ABS(TN-TO).LT.TOL) THEN
          EXIT  
        ENDIF 
      ENDDO !nn
      TOFSIG_5 = TN
      END
      REAL*8 FUNCTION SIG_7(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_12term.h'
      SIG_7 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_7(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      INCLUDE '../../include/stmt_fns_SIGMA0_12term.h'
      SOFSIG_7 = SOFSIG(RR8,TT8)
      END
      REAL*8 FUNCTION TOFSIG_7(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_12term.h'
      TOFSIG_7 = TOFSIG(RR8,SS8)
      END
      REAL*8 FUNCTION SIG_2(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_7term.h'
      SIG_2 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_2(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      INCLUDE '../../include/stmt_fns_SIGMA2_7term.h'
      SOFSIG_2 = SOFSIG(RR8,TT8)
      END
      REAL*8 FUNCTION TOFSIG_2(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_7term.h'
      TOFSIG_2 = TOFSIG(RR8,SS8)
      END
      REAL*8 FUNCTION SIG_4(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_9term.h'
      SIG_4 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_4(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      INCLUDE '../../include/stmt_fns_SIGMA2_9term.h'
      SOFSIG_4 = SOFSIG(RR8,TT8)
      END
      REAL*8 FUNCTION TOFSIG_4(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_9term.h'
      TOFSIG_4 = TOFSIG(RR8,SS8)
      END
      REAL*8 FUNCTION SIG_6(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_17term.h'
      SIG_6 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_6(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      REAL*8, PARAMETER :: TOL=1.D-6
      INTEGER NN
      REAL*8  SN,SO
      REAL*8  SOFSIG_8
      INCLUDE '../../include/stmt_fns_SIGMA2_17term.h'
C     sofsig via Newton iteration from a 12-term 1st guess
      SN = SOFSIG_8(RR8,TT8)  !non-negative
      DO NN= 1,10
        SO = SN
        SN = SO - (SIG(TT8,SO)-RR8)/DSIGDS(TT8,SO)
        IF     (NN.EQ.10 .OR. ABS(SN-SO).LT.TOL) THEN
          EXIT
        ENDIF
      ENDDO !nn
      SOFSIG_6 = SN
      END
      REAL*8 FUNCTION TOFSIG_6(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      REAL*8, PARAMETER :: TOL=1.D-6
      INTEGER NN
      REAL*8  TN,TO
      REAL*8  TOFSIG_8
      INCLUDE '../../include/stmt_fns_SIGMA2_17term.h'
C     sofsig via Newton iteration from a 12-term 1st guess
      TN = TOFSIG_8(RR8,SS8)  !non-negative
      DO NN= 1,10
        TO = TN
        TN = TO - (SIG(TO,SS8)-RR8)/DSIGDT(TO,SS8)
        IF     (NN.EQ.10 .OR. ABS(TN-TO).LT.TOL) THEN
          EXIT  
        ENDIF 
      ENDDO !nn
      TOFSIG_6 = TN
      END
      REAL*8 FUNCTION SIG_8(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_12term.h'
      SIG_8 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_8(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      INCLUDE '../../include/stmt_fns_SIGMA2_12term.h'
      SOFSIG_8 = SOFSIG(RR8,TT8)
      END
      REAL*8 FUNCTION TOFSIG_8(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_12term.h'
      TOFSIG_8 = TOFSIG(RR8,SS8)
      END
      REAL*8 FUNCTION SIG_46(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA4_17term.h'
      SIG_46 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_46(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      REAL*8, PARAMETER :: TOL=1.D-6
      INTEGER NN
      REAL*8  SN,SO
      REAL*8  SOFSIG_48
      INCLUDE '../../include/stmt_fns_SIGMA4_17term.h'
C     sofsig via Newton iteration from a 12-term 1st guess
      SN = SOFSIG_48(RR8,TT8)  !non-negative
      DO NN= 1,10
        SO = SN
        SN = SO - (SIG(TT8,SO)-RR8)/DSIGDS(TT8,SO)
        IF     (NN.EQ.10 .OR. ABS(SN-SO).LT.TOL) THEN
          EXIT
        ENDIF
      ENDDO !nn
      SOFSIG_46 = SN
      END
      REAL*8 FUNCTION TOFSIG_46(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      REAL*8, PARAMETER :: TOL=1.D-6
      INTEGER NN
      REAL*8  TN,TO
      REAL*8  TOFSIG_48
      INCLUDE '../../include/stmt_fns_SIGMA4_17term.h'
C     sofsig via Newton iteration from a 12-term 1st guess
      TN = TOFSIG_48(RR8,SS8)  !non-negative
      DO NN= 1,10
        TO = TN
        TN = TO - (SIG(TO,SS8)-RR8)/DSIGDT(TO,SS8)
        IF     (NN.EQ.10 .OR. ABS(TN-TO).LT.TOL) THEN
          EXIT  
        ENDIF 
      ENDDO !nn
      TOFSIG_46 = TN
      END
      REAL*8 FUNCTION SIG_48(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA4_12term.h'
      SIG_48 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_48(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      INCLUDE '../../include/stmt_fns_SIGMA4_12term.h'
      SOFSIG_48 = SOFSIG(RR8,TT8)
      END
      REAL*8 FUNCTION TOFSIG_48(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA4_12term.h'
      TOFSIG_48 = TOFSIG(RR8,SS8)
      END

      subroutine FieldArchive3D(field,idm,jdm,kz,
     &  cfld,fldname,coord,tlevel,tlevel1,nrec,filebase)
      implicit none
      integer :: k,nrec
      integer,          intent(in)  :: idm,jdm,kz,tlevel
      character(len=*), intent(in)  :: fldname,cfld(nrec)
      character(len=*), intent(in)  :: filebase
      integer,          intent(in)  :: coord(nrec),tlevel1(nrec)
      real,             intent(out) :: field(kz,idm,jdm)
      real    field2D(idm,jdm)

      do k=1,kz
          call FieldArchive(field2D,idm,jdm,
     &    cfld,fldname,k,coord,1,tlevel1,nrec,filebase)
          field(k,:,:) = field2D
      end do
      end subroutine


      recursive subroutine FieldArchive(field,idm,jdm,
     & cfld,fldname,localcoord,coord,tlevel,tlevel1,nrec,filebase)
      implicit none
      integer   nrec
      integer,          intent(in)  :: idm,jdm,coord(nrec),tlevel1(nrec)
      integer         , intent(in)  :: localcoord,tlevel
      real,             intent(out) :: field(idm,jdm)
      character(len=*), intent(in)  :: fldname,cfld(nrec)
      character(len=*), intent(in)  :: filebase
      real*4 :: A(idm,jdm), AMN, AMX, spval,undef
      integer :: indx
      call indexFromH(cfld,fldname,coord,localcoord,tlevel1,tlevel,
     & indx,nrec)
        if (indx/=-1) then
            spval=1e30 !! CAREFUL HERE, ORIGINALLY IT WAS
                       !! spval = undef
                       !! COULDN'T COMPILE LIKE THAT
            call READRAW(A,AMN,AMX,IDM,JDM,.false.,spval,
     &                  trim(filebase),indx)
            field=A
        else
            print '(A)', 'Could not get field "'
     &            //fldname//'" at coordinate:'
            print '(i3.3)',localcoord
            field(:,:)=spval
        end if
      end subroutine 
      
      subroutine indexFromH(cfld,fldname,coord,localcoord,
     &           tlevel1,tlevel,indx,nrec)
      implicit none
      integer :: irec,nrec
      character(len=*), intent(in)  :: cfld(nrec),fldname
      integer         , intent(in)  :: coord(nrec),tlevel1(nrec)
      integer         , intent(in)  :: localcoord,tlevel
      integer         , intent(out) :: indx
      indx=-1
      do irec=1,nrec
C --- BEGIN: MOSTAFA
c      if (trim(fldname)==trim(cfld(irec)) .and.
c     &   tlevel1(irec)==tlevel) then
c         indx=irec
c      end if
C --- END: MOSTAFA
         if (trim(fldname)==trim(cfld(irec)) .and.
     &   localcoord==coord(irec) .and. tlevel1(irec)==tlevel ) then
            indx=irec
         end if
      end do
      end subroutine

C
C
C
      SUBROUTINE READRAW(A,AMN,AMX,IDM,JDM,LSPVAL,SPVAL,CFILE1,K)
      IMPLICIT NONE
!
      REAL*4     SPVALH
      PARAMETER (SPVALH=2.0**99)
c!
      REAL*4,           INTENT(OUT) :: A(IDM,JDM)
      REAL*4,           INTENT(OUT) :: AMN,AMX
      INTEGER,          INTENT(IN)  :: IDM,JDM
      LOGICAL,          INTENT(IN)  :: LSPVAL
      REAL*4,           INTENT(INOUT)  :: SPVAL
      INTEGER,          INTENT(IN)  :: K
      CHARACTER(len=*), INTENT(IN)  :: CFILE1
!
      REAL*4 :: PADA(4096)
!     MOST OF WORK IS DONE HERE.
!

      INTEGER      I,J,IOS,NRECL
      INTEGER NPAD,GET_NPAD
!
      IF(.NOT.LSPVAL) THEN
        SPVAL = SPVALH
      ENDIF
!
!!! Calculate the number of elements padded!!!!!!!!!!!!!!!!!!!!!!!!
      NPAD=GET_NPAD(IDM,JDM)
      INQUIRE( IOLENGTH=NRECL) A,PADA(1:NPAD)
!
      OPEN(UNIT=73, FILE=CFILE1, FORM='UNFORMATTED', STATUS='old',
     &     ACCESS='DIRECT', RECL=NRECL,IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILE1(1:LEN_TRIM(CFILE1))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
!
      READ(73,REC=K,IOSTAT=IOS) A
      close(73)
!
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read record ',K,
     &             ' from '//CFILE1(1:LEN_TRIM(CFILE1))
        CALL EXIT(4)
      ENDIF
!
      AMN =  SPVALH
      AMX = -SPVALH
      DO J= 1,JDM
      DO I=1,IDM
         IF     (A(I,J).LE.SPVALH) THEN
            AMN = MIN( AMN, A(I,J) )
            AMX = MAX( AMX, A(I,J) )
         ELSEIF (LSPVAL) THEN
            A(I,J) = SPVAL
         ENDIF
      END DO
      END DO
!                 
      RETURN
      END SUBROUTINE

      INTEGER FUNCTION GET_NPAD(IDM,JDM)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IDM,JDM
         GET_NPAD = 4096 - MOD(IDM*JDM,4096)
         GET_NPAD = mod(GET_NPAD,4096)
      END FUNCTION
C 
C 
C
       subroutine GRIRD_REMAPPING(i,j,K,KDM,KZ,NSIGMA,SIG3D,RM,PM,RZ,
     +  DEPTH,DSCK,DPCK,ISOTOP,ISOPYC,PKM1,PZBOT,PZTOP,RKM1,ZL,PMIX,
     +  PCM0,PCM1,PCM2,PKM2,DP0K,DP00I,DS0K,flag_t,flag_s,
     +          flag_u,flag_v,itest,jtest,SM,TM,UM,VM,UV,VV,
     +          ldebug,sigver,TZ,SZ,
     +  flag_no3,flag_po4,flag_si,NO3,PO4,SI,NO3M,PO4M,SIM)

        IMPLICIT NONE
         
        integer :: K,KDM,kz,i,j,NSIGMA,itest,jtest,sigver
        real*4, intent(in)    :: RZ(kz+1),SIG3D(kdm),PMIX,
     +         depth, DSCK(0:kdm+1),   DPCK(0:kdm+1),
     +         ISOTOP,PKM1,PZTOP,RKM1,ZL(kz+1),
     +         DP0K(kdm+1),DP00I,DS0K(kdm+1),
     +         UV(kz+1),VV(kz+1),TZ(kz+1),SZ(kz+1),
     +         NO3(kz+1),PO4(kz+1),SI(kz+1)

        real*4, intent(inout) ::  PM,RM
        real*4,  intent(inout) :: PZBOT,PCM0,PkM2,PCM1,PCM2,
     +         UM,VM,SM,TM,NO3M,PO4M,SIM
       logical    ::  ISOPYC,ldebug
        CHARACTER*40    :: flag_t,flag_s,flag_u,flag_v
       CHARACTER*40    :: flag_no3,flag_po4,flag_si

       real*4 ::qdep,Q,sigmaa,RZLOC,DPMS,PZMID,THIKMN,THK
       integer :: L,kztop
       REAL*4     ZERO,ONE,SOFSIG_V,dmin
       PARAMETER (ZERO=0.0, ONE=1.0)
C
C               BEGIN SUBROUTINE
C
                IF     (K.GT.1) THEN
                RM =     SIG3D(K)
              ELSE
                RM = MIN(SIG3D(1),RZ(1) )
              ENDIF
               if (i.eq.600 .and. j.eq.2) then
                      WRITE(6,'(A,F5.2,A,I2)')
     +                  'RM: ilk loop =',RM,'  K=',K
               endif !debug

C
C             FIND RM AND PM (I).
C



              QDEP = MAX( 0.0, MIN( 1.0,
     +                    (DEPTH  - DSCK(NSIGMA)) /
     +                    (DPCK(NSIGMA) - DSCK(NSIGMA))  ) )
              DMIN = (1.0-QDEP)*DSCK(K-1) + QDEP*DPCK(K-1)
                if (i.eq.600 .and. j.eq.2) then
                  WRITE(6,'(A,I3,F10.3,F8.3)')
     +              'ISOTOP - K,DMIN =',K,DMIN,ISOTOP
                endif !debug
              IF     (ISOPYC .AND. K.EQ.1) THEN
C
C               UPPER LAYER IS THE MIXED LAYER.
C
                PM = MIN( DEPTH, PMIX )
              ELSEIF (K.EQ.KDM) THEN
C
C               LOWEST LAYER.
C
                PM = DEPTH
                if (i.eq.600 .and. j.eq.2) then
                  WRITE(6,'(A,I3,A,F8.3,A,F8.3)')
     +              'LOWEST LAYER =',K,'  PM=',PM,
     +              '  DEPTH=',DEPTH
                endif !debug
              ELSEIF (DMIN.LE.ISOTOP) THEN
C
C               FIXED GRID NEAR THE SURFACE.
C

                DMIN    = (1.0-QDEP)*DSCK(K) + QDEP*DPCK(K)
                PM = MIN( DEPTH, MAX( PKM1, DMIN ) )
                  if (i.eq.600 .and. j.eq.2) then
                    WRITE(6,'(A,I3,F10.3,A,F7.3,A,F7.3)')
     +                'FIXED - K,PM =',K,PM,' PKM1=',PKM1
                    WRITE(6,'(A,F8.3)')' PZBOT=',PZBOT,
     +              ' PZTOP=',PZTOP
                  endif !debug
              ELSEIF (RM.LT.RKM1) THEN
C
C               ZERO THICKNESS LAYER REQUIRED FOR STABILITY.
C
                PM = PKM1
                  if (i.eq.600 .and. j.eq.2) then
                    WRITE(6,'(A,I3,F10.3,A,F8.3)')
     +            'ZERO THICKNESS - K,PM =',K,PM,' PKM1=',PKM1
                    WRITE(6,'(A,I3,F10.3,A,F8.3)')
     +           'ZERO THICKNESS - K,RM =',K,RM,' RKM1=',RKM1
                    WRITE(6,'(A,F8.3,A,F8.3)')' PZBOT=',PZBOT,
     +              ' PZTOP=',PZTOP
                  endif !debug

              ELSE  !K.LT.KDM
C
C               INTERFACE IS HALF WAY BETWEEN TARGET DENSITIES.
C

                SIGMAA = 0.5*(SIG3D(K) + SIG3D(K+1))

                KZTOP  = PZTOP

                CALL FIND_DENSITY(PZBOT,
     +                            SIGMAA,  RZ,KZ+1,KZTOP)

                IF     (PZBOT.EQ.0.0) THEN
                  PM = PKM1
                ELSEIF (PZBOT.EQ.KZ+1) THEN
                  PM = DEPTH
                ELSE
                  L = PZBOT
                  Q = PZBOT - L
                  PM = MIN( DEPTH,
     +                        MAX( PKM1,
     +                         (1.0-Q)*ZL(L) + Q*ZL(L+1) ) )

                ENDIF
                  if (i.eq.600 .and. j.eq.2) then
                    WRITE(6,'(A,F10.3,A,F7.3)')
     +                'PM=',PM,' PZBOT=',PZBOT
                    WRITE(6,'(A,I3,A,F7.3,A,F7.3)')'KZTOP=',KZTOP,
     +                ' SIGMAA=',SIGMAA,' RZ',RZ(1)
                  endif !debug
              ENDIF


C
C             MODIFY RM AND PM (IF NECESSARY).  (II)
C


              KZTOP  = PZTOP
              CALL FIND_DEPTH(PZBOT,
     +                        PM, ZL,KZ+1,KZTOP)
              
                if (i.eq.itest .and. j.eq.jtest) then
                  WRITE(6,'(A,2F10.3)')
     +              'PM,PZBOT =',PM,PZBOT
                endif !debug
C
              SIGMAA = SIG3D(K)
              KZTOP  = PZTOP
              CALL FIND_DENSITY(RZLOC,
     +                          SIGMAA,  RZ,KZ+1,KZTOP)
                if (i.eq.itest .and. j.eq.jtest) then
                  WRITE(6,'(A,2F10.3)')
     +              'RZ,RZLOC =',SIGMAA,RZLOC
                endif !debug
              IF     (RZLOC.LT.PZTOP .OR.
     +                RZLOC.GT.PZBOT     ) THEN
                PZMID = 0.5*(PZTOP + PZBOT)
              ELSE
                PZMID = RZLOC
              ENDIF
                if (i.eq.itest .and. j.eq.jtest) then
                  WRITE(6,'(A,2F10.3)')
     +              'XX,PZMID =',SIGMAA,PZMID
                endif !debug
C
              IF     (DMIN.LE.ISOTOP) THEN
                DPMS      = PM
                PCM0 = DPMS
              ELSE
                IF     (DP0K(K).LE.DP00I) THEN
                  Q = DP0K(K)
                ELSE
                  Q  = MAX( DP00I,
     &                      DP0K(K) * DP0K(K)/
     &                                MAX( DP0K( K),
     &                                     PKM1-PCM1 ) )
                ENDIF
                PCM0 = PCM1 +
     +                        MIN(Q, (1.0-QDEP)*DS0K(K)+QDEP*DP0K(K))
                DPMS      = PKM1 +
     +                        MIN(Q, (1.0-QDEP)*DS0K(K)+QDEP*DP0K(K))
                  if (i.eq.itest .and. j.eq.jtest) then
                    WRITE(6,'(A,3F10.3)')
     +                'PKM1,PCM1,S =',PKM1,PCM1,
     +                                DP0K(K)/
     +                                MAX( DP0K( K),
     +                                     PKM1-PCM1 )
                    WRITE(6,'(A,3F10.3)')
     +                'DP0K,Q,DPMS =',DP0K(K),Q,DPMS
                  endif !debug
              ENDIF
              IF     (.NOT. ISOPYC .AND.
     +                PM.EQ.DEPTH) THEN
                IF     (PM.GT.PKM1) THEN
C
C                 MAKE LOWEST NON-ZERO LAYER EXACTLY ISOPYNCAL,
C                 UNLESS THE LAYER ABOVE IS NOT ISOPYCNAL.
C
                  IF     (K.GT.1) THEN
                    IF     (RKM1.NE.SIG3D(K-1)) THEN
                      THK    = PKM1-PKM2
                      IF     (DP0K(K-1).LE.DP00I) THEN
                        Q = DP0K(K-1)
                      ELSE
                        Q  = MAX( DP00I,
     &                            DP0K(K-1) * DP0K(K-1)/
     &                                      MAX( DP0K(K-1),
     &                                           PKM2-PCM2 ) )
                      ENDIF
                      THIKMN = MIN(Q, (1.0-QDEP)*DS0K(K)+QDEP*DP0K(K))
                    ELSE  ! LAYER ABOVE IS DEFINATELY ISOPYCNAL
                      THK    = 1.0
                      THIKMN = 0.0
                    ENDIF
                  ELSE
                    THK    = 1.0
                    THIKMN = 1.0
                  ENDIF
                  IF     (        THIKMN .EQ.ZERO .OR.
     +                    ABS(THK-THIKMN).GT.0.01     ) THEN
                    RM = MAX(RKM1,SIG3D(K))
                  ELSE
C
C                   LOWEST NON-ZERO LAYER IS NOT-ISOPYCNAL.
C
                    L = PZMID
                    Q = PZMID - L
                    RM = (1.0-Q)*RZ(L) + Q*RZ(L+1)
                      if (i.eq.itest .and. j.eq.jtest) then
                        WRITE(6,'(A,I4,F8.4)')
     +                    'LRM: L,Q =',L,Q
                      endif !debug
                  ENDIF
                ELSE
                  RM =     RKM1
                ENDIF
              ELSEIF (.NOT. ISOPYC .AND.
     +                PM.LE.MIN(DEPTH,DPMS)) THEN
C
C               HYBRID (MINIMUM THICKNESS) LAYER.
C
                PM = MIN(DPMS,DEPTH)
                KZTOP  = PZTOP
                CALL FIND_DEPTH(PZBOT,
     +                          PM, ZL,KZ+1,KZTOP)
                  if (i.eq.itest .and. j.eq.jtest) then
                    WRITE(6,'(A,2F10.3)')
     +                'PM,PZBOT =',PM,PZBOT
                  endif !debug
                IF     (RZLOC.LT.PZTOP .OR.
     +                  RZLOC.GT.PZBOT     ) THEN
                  PZMID = 0.5*(PZTOP + PZBOT)
                ELSE
                  PZMID = RZLOC
                ENDIF
                  if (i.eq.itest .and. j.eq.jtest) then
                    WRITE(6,'(A,2F10.3)')
     +                'PM,PZMID =',PM,PZMID
                  endif !debug
                  if (i.eq.itest .and. j.eq.jtest) then
                    WRITE(6,'(A,I3,2F10.3)')
     +                'HYBRID - K,PKM1,PM =',K,PKM1,PM
                  endif !debug
                IF     (PM.GT.PKM1) THEN
                  L = PZMID
                  Q = PZMID - L
                  RM = (1.0-Q)*RZ(L) + Q*RZ(L+1)
                    if (i.eq.itest.and. j.eq.jtest) then
                      WRITE(6,'(A,I4,2F8.4)')
     +                  'HRM: L,Q =',L,Q,RM
                    endif !debug
C               ELSE
C                 RM UNCHANGED
                ENDIF
              ENDIF
C
C             FIND T & S.
C
              IF     (PM.GT.PKM1) THEN
                L = PZMID
                Q = PZMID - L
                TM = (1.0-Q)*TZ(L) + Q*TZ(L+1) !
                !SM = (1.0-Q)*SZ(L) + Q*SZ(L+1)
                if(flag_u.ne."NONE")
     &                 UM=(1.0-Q)*UV(L)+Q*UV(L+1) ! from Z to isop.
                if(flag_v.ne."NONE")
     &                 VM=(1.0-Q)*VV(L)+Q*VV(L+1) !
                if(flag_no3.ne."NONE")
     &                 NO3M=(1.0-Q)*NO3(L)+Q*NO3(L+1) ! from Z to isop.
                if(flag_po4.ne."NONE")
     &                 PO4M=(1.0-Q)*PO4(L)+Q*PO4(L+1) ! from Z to isop.
                if(flag_si.ne."NONE")
     &                 SIM=(1.0-Q)*SI(L)+Q*SI(L+1) ! from Z to isop.
                  if (ldebug.and.i.eq.itest.and. j.eq.jtest) then
                    WRITE(6,'(A,I4,3F8.4)')
     +                ' TM: L,Q =',L,Q,TM,RM
                  endif !debug
C             ELSE
C               TM UNCHANGED
              ENDIF
              !L = PZMID
              !Q = PZMID - L
              !SM = MIN((1.0-Q)*SZ(L) + Q*SZ(L+1),SOFSIG_V(RM,TM,SIGVER))  !! MOSTAFA
              SM=SOFSIG_V(RM,TM,SIGVER)


        end
        
C
C      READ 2D and 3D fields from ARCHIVE FILE
C


        SUBROUTINE READ_ARCHIVE(KZ,IDM,JDM,NREC,KDM,COORD,
     +     tlevel1,ARCHVFILE,cfld,MTG1,SRFH,SUFLX,SAFLX,
     +     BLDP,MIXDP,UBRTP,VBRTP,TZ,SZ,UV,VV,RZ,ZL,
     +     flag_t,flag_s,flag_th,flag_u,flag_v,ISOPYC,
     +     SIGVER,LEVTOP,SIG3D,DEPTH,TZ1,SZ1,UV1,VV1,TH,TH1,
     +     flag_no3,flag_po4,flag_si,NO3,PO4,SI)
        IMPLICIT NONE
C
        INTEGER :: nrec,KZ,IDM,JDM,KDM,
     +          SIGVER,LEVTOP
        INTEGER  :: coord(nrec),tlevel1(nrec)
        REAL*4, intent(Out):: MTG1(IDM,JDM),SRFH(IDM,JDM),
     +              SAFLX(IDM,JDM),BLDP(IDM,JDM),MIXDP(IDM,JDM),
     +              SUFLX(IDM,JDM),UBRTP(IDM,JDM),VBRTP(IDM,JDM)
        REAL*4, intent(in) ::  ZL(KZ+1),SIG3D(IDM,JDM,KDM+1),
     +         DEPTH(IDM,JDM)
        REAL*4, intent(out) ::TZ1(IDM,JDM),SZ1(IDM,JDM),
     +         UV1(IDM,JDM),VV1(IDM,JDM),UV(KZ+1,IDM,JDM),
     +         VV(KZ+1,IDM,JDM),RZ(KZ+1,IDM,JDM),
     +         TH(KZ+1,IDM,JDM),SZ(KZ+1,IDM,JDM),
     +         TH1(IDM,JDM),TZ(KZ+1,IDM,JDM),
     +         NO3(KZ+1,IDM,JDM),PO4(KZ+1,IDM,JDM),SI(KZ+1,IDM,JDM)


        REAL*4          :: SIG_V,TOFSIG_V,spval,spval_
        REAL*4          :: PU(IDM,JDM),PV(IDM,JDM)
        LOGICAL         :: ISOPYC

        CHARACTER*40, INTENT(IN)::flag_t,flag_s,flag_u,
     +                      flag_v,flag_th,
     +                      flag_no3,flag_po4,flag_si
        CHARACTER*256, INTENT(IN) :: ARCHVFILE
        CHARACTER(len=8), INTENT(IN) :: cfld(NREC)
        REAL*4          :: NO1(IDM,JDM),PO1(IDM,JDM),SI1(IDM,JDM)
C
        INTEGER   ::  K,I,J

C
C        2D FIELDS
C
         spval_=2.00**99
         spval=2.00**99


         call FieldArchive(MTG1,IDM,JDM,cfld,'montg1    ',
     &     1,coord,1,tlevel1,nrec,trim(archvfile)//".a")

         call FieldArchive(SRFH,IDM,JDM,cfld,'srfhgt    ',
     &     0,coord,1,tlevel1,nrec,trim(archvfile)//".a")

         call FieldArchive(SUFLX,IDM,JDM,cfld,'surflx    ',
     &     0,coord,1,tlevel1,nrec,trim(archvfile)//".a")

         call FieldArchive(SAFLX,IDM,JDM,cfld,'salflx    ',
     &     0,coord,1,tlevel1,nrec,trim(archvfile)//".a")

         call FieldArchive(BLDP,IDM,JDM,cfld,'bl_dpth    ',
     &     0,coord,1,tlevel1,nrec,trim(archvfile)//".a")

         call FieldArchive(MIXDP,IDM,JDM,cfld,'mix_dpth    ',
     &     0,coord,1,tlevel1,nrec,trim(archvfile)//".a")

         call FieldArchive(UBRTP,IDM,JDM,cfld,'u_btrop    ',
     &     0,coord,1,tlevel1,nrec,trim(archvfile)//".a")

         call FieldArchive(VBRTP,IDM,JDM,cfld,'v_btrop    ',
     &     0,coord,1,tlevel1,nrec,trim(archvfile)//".a")


C
C        3D FIELDS
C


      DO K=1,KZ ! LOOP THROUGH LEVELS
C
C
C
        if     (flag_u.ne."NONE") then
         call FieldArchive(UV1,IDM,JDM,cfld,'u-vel.    ',
     &     K,coord,1,tlevel1,nrec,trim(archvfile)//".a")
        endif
        if     (flag_v.ne."NONE") then
         call FieldArchive(VV1,IDM,JDM,cfld,'v-vel.    ',
     &     K,coord,1,tlevel1,nrec,trim(archvfile)//".a")
        endif
        if     (flag_th.ne."NONE") then
         call FieldArchive(TH1,IDM,JDM,cfld,'thknss    ',
     &     K,coord,1,tlevel1,nrec,trim(archvfile)//".a")
        endif
        if     (flag_t.ne."NONE") then
         call FieldArchive(TZ1,IDM,JDM,cfld,'temp    ',
     &     K,coord,1,tlevel1,nrec,trim(archvfile)//".a")
        endif
        if     (flag_s.ne."NONE") then
         call FieldArchive(SZ1,IDM,JDM,cfld,'salin    ',
     &     K,coord,1,tlevel1,nrec,trim(archvfile)//".a")
        endif
        if     (flag_no3.ne."NONE") then
         call FieldArchive(NO1,IDM,JDM,cfld,'ECO_no3    ',
     &     K,coord,1,tlevel1,nrec,trim(archvfile)//".a")
        endif
        if     (flag_po4.ne."NONE") then
         call FieldArchive(PO1,IDM,JDM,cfld,'ECO_pho    ',
     &     K,coord,1,tlevel1,nrec,trim(archvfile)//".a")
        endif 
        if     (flag_si.ne."NONE") then
         call FieldArchive(SI1,IDM,JDM,cfld,'ECO_sil    ',
     &     K,coord,1,tlevel1,nrec,trim(archvfile)//".a")
        endif


        ! CAUTION:  CONVERT VELOCITY COMPONENTS FROM P-CELL TO U- AND V-CELL, RESPECTIVELY HERE
        ! TODO:  
        CALL UV2P(IDM,JDM,UV1,UV1,PU,PV)


         ! CAUTION: for some reason temperature at level-75 is FILL_value
         !          but salinity is not. check the code
         !          READING NOT DONE PROPERLY ??
         !
         ! THE ORIGINAL "RELAXI_ARCHV.F" CODE DOES CHECKS ON HMINA,HMINB AND
         ! ZLEVS AMONG DIFFERENT INPUT FILES. I SKIPPED THOSE HERE BECAUSE
         ! WE ARE READING SINGLE FILE, SO PRESUMABLY ZLEVS WILL BE CONSISTENT.     
         ! CODE CHECKS IF HMINB AND HMAXB ARE EQUAL, MEANING THAT THE LEVEL IS
         ! ALL CONSTANT VALUES, AND ASSIGNS TO EVERY POINT HMINB. LATER ON IT
         ! ALSO CHECKS MINS AND MAXS FROM .A FILE AND COMPARES THEM WITH .B
         ! ONES. IF INCONSISTENT, STOPS THE PROGRAM. MAYBE THESE SHOULD BE
         ! ADDED LATER.

C ---    U AND V ARE NOW ON P-CELL
         DO J= 1,JDM
            DO I= 1,IDM               ! assign 3D
                 if     (flag_t.ne."NONE")
     &                  TH(K,I,J) = MIN(SPVAL_,TH1(I,J)) ! layer thickness
                 if     (flag_s.ne."NONE")
     &                  SZ(K,I,J) = MIN(SPVAL_,SZ1(I,J)) ! salinity
                 if     (flag_u.ne."NONE")
     &                  UV(K,I,J) = MIN(SPVAL_,PU(I,J)) ! u-velocity

                 if     (flag_v.ne."NONE")
     &                  VV(K,I,J) = MIN(SPVAL_,PV(I,J)) ! v-velocity

                 if     (flag_no3.ne."NONE")
     &                  NO3(K,I,J) = MIN(SPVAL_,NO1(I,J)) ! u-velocity
                 if     (flag_po4.ne."NONE")
     &                  PO4(K,I,J) = MIN(SPVAL_,PO1(I,J)) ! u-velocity
                 if     (flag_si.ne."NONE")
     &                  SI(K,I,J) = MIN(SPVAL_,SI1(I,J)) ! v-velocity
            ENDDO
         ENDDO

         ! THESE LOOPS HERE CALCULATE DENSITY FROM T AND S
         IF     (K.EQ.1) THEN
            DO J= 1,JDM
               DO I= 1,IDM
                    TZ(K,I,J) = TZ1(I,J) ! temperature is retrieved here
                                         ! in the original code
                    RZ(K,I,J) = SIG_V(TZ(K,I,J),SZ(K,I,J),SIGVER)
               END DO
            END DO

            ELSEIF (.NOT.ISOPYC) THEN !
                                      ! BLKDAT I USE SELECTS THIS CONDITION
                                      !
            ! RZ MUST BE MONOTONICALLY NON-DECREASING (NEAR THE BOTTOM).
               DO J= 1,JDM
                  DO I= 1,IDM
                       TZ(K,I,J) = TZ1(I,J)
                       RZ(K,I,J) = SIG_V(TZ(K,I,J),SZ(K,I,J),SIGVER)   
                     IF (RZ(K,I,J).LT.RZ(K-1,I,J) .AND.
     &                  ZL(MIN(K+3,KZ+1)).GE.DEPTH(I,J)) THEN
                        RZ(K,I,J) = RZ(K-1,I,J)
                        TZ(K,I,J) = TZ(K-1,I,J)
                     ENDIF
                  END DO
               END DO

               ELSE
               ! LIMIT MAXIMUM DENSITY TO SIGMA(KDM)
               DO J= 1,JDM
                  DO I= 1,IDM
                       TZ(K,I,J) = MIN(SPVAL_,TZ1(I,J))

!                        IF     (MAX(TZ(K,I,J),
!     &                        SZ(K,I,J) ).GT.2.0**90) THEN
!                            RZ(K,I,J) = SPVAL
!                        else
                            RZ(K,I,J)=SIG_V(TZ(K,I,J),SZ(K,I,J),SIGVER)
!                        endif

                       IF (RZ(K,I,J).LT.RZ(K-1,I,J)) THEN
                          RZ(K,I,J) = RZ(K-1,I,J)
                          TZ(K,I,J) = TZ(K-1,I,J)
                       END IF
                       IF (RZ(K,I,J).GT.SIG3D(I,J,KDM)) THEN
                         IF (RZ(MAX(K-1,1),I,J).EQ.SIG3D(I,J,KDM)) THEN
                            RZ(K,I,J) = SIG3D(I,J,KDM)
                            TZ(K,I,J) = TZ(K-1,I,J)
                            SZ(K,I,J) = SZ(K-1,I,J)
                            if     (flag_u.ne."NONE")
     &                                  UV(K,I,J) = UV(K-1,I,J) ! SINCE SZ SHOWED UP HERE IN THE
                            if     (flag_v.ne."NONE")
     &                                  VV(K,I,J) = VV(K-1,I,J) ! ORIGINAL CODE, I INCLUDED UV AND
                                                    ! VV AS WELL. DOUBLE CHECK IF ITS OK
                            if     (flag_no3.ne."NONE")
     &                                  NO3(K,I,J) = NO3(K-1,I,J)
                            if     (flag_po4.ne."NONE")
     &                                  PO4(K,I,J) = PO4(K-1,I,J)
                            if     (flag_si.ne."NONE")
     &                                  SI(K,I,J) = SI(K-1,I,J)

                            ELSE
                            RZ(K,I,J) = SIG3D(I,J,KDM)
                            TZ(K,I,J) = TOFSIG_V(RZ(K,I,J),
     &                                  SZ(K,I,J),SIGVER)
                         ENDIF
                       END IF
                  END DO
               END DO
           END IF
 
      END DO ! END OF LEVEL LOOP


      IF     (LEVTOP.GT.1) THEN
        DO J= 1,JDM
          DO I= 1,IDM
            DO K= 1,LEVTOP-1
              RZ(K,I,J) = RZ(LEVTOP,I,J)
              TZ(K,I,J) = TZ(LEVTOP,I,J)
            ENDDO
          ENDDO
        ENDDO
      ENDIF


      DO J= 1,JDM
        DO I= 1,IDM
          RZ(KZ+1,I,J) = RZ(KZ,I,J) + 0.001
          TZ(KZ+1,I,J) = TZ(KZ,I,J)
        ENDDO
      ENDDO



        END SUBROUTINE READ_ARCHIVE
        
        
c
c --- create u and v grids for a given p-grid matrix.
c
c --- for compatibility:
c ---   idm,jdm are input from regional.grid.b,
c ---   AUTHOR: MOSTAFA BAKHODAY-PASKYABI
c
c --- define the 4 staggered grids.
c
        SUBROUTINE P2UV(IDM,JDM,PU,PV,UU,VV)
        IMPLICIT NONE
        REAL*4 , INTENT(INOUT) :: PU(IDM,JDM),PV(IDM,JDM)
        INTEGER,  INTENT(IN):: IDM,JDM
        REAL*4 , INTENT(OUT) :: UU(IDM,JDM),VV(IDM,JDM)
        REAL*4,    parameter   :: hspval=0.5*2.0**100  ! half spval
        INTEGER   :: I,J
C
            do j= 1,jdm
                do i= 1,idm
c
                    if (PU(i,j)  .gt. hspval) then
                        PU(i,j) = 0.0
                    else
                        if     (j.ne.1) then
                            UU(i,j) = 0.5*(PU(i,j-1) + PU(i,j))
                        else
                            UU(i,j) = 2.0*PU(i,2) - PU(i,3)
                        endif
                    ENDIF
c
                    if (PV(i,j)  .gt. hspval) then
                        PV(i,j) = 0.0
                    else
                        if     (i.ne.1) then
                            VV(i,j) = 0.5*(PV(i-1,j) + PV(i,j))
                        else
                            VV(i,j) = 2.0*PV(2,j) - PV(3,j)
                        endif
                    ENDIF
                enddo
            enddo

        END SUBROUTINE P2UV


c
c --- create A P-GRID from specific u and v grids.
c
c --- for compatibility:
c ---   idm,jdm are input from regional.grid.b,
c ---   AUTHOR: MOSTAFA BAKHODAY-PASKYABI
c
c --- define the 4 staggered grids.
c
        SUBROUTINE UV2P(IDM,JDM,UU,VV,PU,PV)
        IMPLICIT NONE
        REAL*4 , INTENT(IN) :: UU(IDM,JDM),VV(IDM,JDM)
        INTEGER,  INTENT(IN):: IDM,JDM
        REAL*4 , INTENT(OUT) :: PU(IDM,JDM),PV(IDM,JDM)
        REAL*4,    parameter   :: hspval=0.5*2.0**100  ! half spval
        INTEGER   :: I,J
C
        do j= 1,jdm
            do i= 1,idm
c
                if (UU(i,j)  .gt. hspval) then
                    PU(i,j) = 0.0
                else
                    if     (j.ne.jdm) then
                        PU(i,j) = 0.5*(UU(i,j) + UU(i,j+1))
                    else
                        PU(i,j) = 2.0*UU(i,j-1) - UU(i,j-2)
                    endif
                endif
c
                if (VV(i,j)  .gt. hspval) then
                    PV(i,j) = 0.0
                else
                    if     (i.ne.idm) then
                        PV(i,j) = 0.5*(VV(i,j) + VV(i+1,j))
                    else
                        PV(i,j) = 2.0*VV(i-1,j) - VV(i-2,j)
                    endif
                endif
c
            enddo
        enddo



        END SUBROUTINE UV2P
c
c ---
c --- Read depth file
c ---
c
        SUBROUTINE READ_DEPTH(KZ,ZL,flnm_z)
        IMPLICIT NONE
C
        INTEGER :: nrec,KZ,K
        REAL*4 :: ZL(KZ+1)
        CHARACTER*240, INTENT(IN) :: flnm_z

        write(6,*)
        write(6,*) 'Open MERCATOR depth file'
        write(6,*)

        call zhopnc(9, flnm_z, 'FORMATTED', 'OLD', 0)
        do K= 1,KZ
            read(9,'(F10.6)') ZL(K)
            write(6,*) 'LEVEL  ',K,'--',ZL(K)
        enddo
        ZL(KZ+1) = ZL(KZ) + 50.0
c        write(6,*) 'LEVEL  ',KZ+1,'--',ZL(KZ+1)

        close(9)
        write(6,*) 'close ',trim(flnm_z)
        call flush(6)

        END SUBROUTINE READ_DEPTH

c
c ---
c --- extract number of vertical layers from archive file [b].
c ---

        SUBROUTINE INIT_SPEC(ARCHVFILE,NREC,KZ,
     +             CFLD,COORD,TLEVEL1,RDAY)
        IMPLICIT NONE
        INTEGER, INTENT(OUT)  :: NREC,KZ
        CHARACTER(len=8), INTENT(OUT)  :: cfld(999)
        INTEGER, INTENT(OUT) :: TLEVEL1(999),COORD(999)
        CHARACTER*256, INTENT(IN) :: ARCHVFILE
        REAL*4, INTENT(out) :: RDAY
        INTEGER     :: I,IOS,tlevel,K,ISTEP
        REAL*4 :: LRDENS,HMINB,HMAXB
        CHARACTER*256 :: CLINE,fldname
C
        CALL ZHOPNC(72, trim(archvfile)//".b", 'FORMATTED', 'OLD', 0)
        kz=-1
c       Read header section
        do i=1,10
            READ (72,'(A)',end=123,err=123) CLINE
        end do
c       Read fields
        nrec=0;
        tlevel = 1
        do
            READ (72,'(A)',end=123,err=123) CLINE
            I = INDEX(CLINE,'=')
            READ (CLINE(I+1:),*) ISTEP,RDAY,K,LRDENS,HMINB,HMAXB
            call zhflsh(6)
            fldname=trim(cline(1:i-1))
            if (trim(fldname)=='salin') then
                kz=max(kz,k)
            end if
            nrec=nrec+1
            coord(NREC)   = K
            tlevel1(NREC) = tlevel
            cfld(NREC)    = fldname

        end do
123     continue
c        nrec=nrec-1
c        print '(a,i4)',"Number of input layers: ",kz
        CLOSE(72)
        END SUBROUTINE INIT_SPEC
c
c ---
c --- assign all variable to extract field values from  archive file [a].
c ---
c
        SUBROUTINE FIELDS_SPEC(ARCHVFILE,CFLD,COORD,TLEVEL1,
     +             NREC)
        IMPLICIT NONE
        INTEGER, INTENT(IN)  :: NREC
        CHARACTER*256, INTENT(IN) :: ARCHVFILE
        CHARACTER(len=8), INTENT(OUT)  :: cfld(NREC)
C
        INTEGER, INTENT(OUT) :: TLEVEL1(NREC),COORD(NREC)
C        REAL, INTENT(OUT) :: MINVALUE(NREC),MAXVALUE(NREC)
        INTEGER     :: I,IOS,tlevel,K,ISTEP,inrec
        REAL*4 :: LRDENS,HMINB,HMAXB,RDAY
        CHARACTER*256 :: CLINE,fldname
C
        CALL ZHOPNC(72, trim(ARCHVFILE)//".b", 'FORMATTED', 'OLD', 0)
C
C     process the archive header
C
        do i=1,10
            READ (72,'(A)') CLINE
        end do
        inrec=0 
        ios=0
        tlevel = 1
        I=0
        do while(ios==0)
            READ (72,'(A)',iostat=ios) CLINE
            I = INDEX(CLINE,'=')
            READ (CLINE(I+1:),*) ISTEP,RDAY,K,LRDENS,HMINB,HMAXB
            fldname       = trim(CLINE(1:i-1))
            inrec          = inrec+1
            coord(iNREC)   = K
            tlevel1(iNREC) = tlevel
            cfld(iNREC)    = fldname
        end do
        iNREC = iNREC -1
        CLOSE(72)
	     END SUBROUTINE FIELDS_SPEC




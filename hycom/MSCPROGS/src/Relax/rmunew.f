      PROGRAM RLXMSK
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
      REAL       SPVAL
      PARAMETER (SPVAL=2.0**100)
C
C     OUTPUT ARRAYS
C
      REAL,  ALLOCATABLE :: RMU(:,:)
C
C     OTHER ARRAYS.
C
      INTEGER, ALLOCATABLE :: IRMU(:,:)
      REAL,  ALLOCATABLE :: DEPTH(:,:)
      REAL,  ALLOCATABLE :: PLON(:,:)
      REAL,  ALLOCATABLE :: PLAT(:,:)
      integer :: rmaxx_i,rmaxx_j
      real :: celldist,rdist,rmaxx
C
C     OTHER VARIABLES.
C
      REAL    RMUMIN,RMUMAX,HMINA,HMINB,HMAXA,HMAXB
      CHARACTER PREAMBL(5)*79,CLINE*80
C
C     NAMELIST INPUT.
C
      CHARACTER*79   CTITLE
      INTEGER        IF(999),IL(999),JF(999),JL(999)
      REAL           EFOLD(999)
      NAMELIST/MASK/ CTITLE,IF,IL,JF,JL,EFOLD
CKAL -- Special points...
      integer nrpoints, ipoint, nrpointsmax
      parameter  (nrpointsmax = 200 )
      real   pointlon(nrpointsmax), pointlat(nrpointsmax)
      real pointrad(nrpointsmax)
      real :: pntlon,pntlat,prad
      character*3 pointname(nrpointsmax)
      logical pointflag(nrpointsmax)
      logical  pflag
      integer ifst,ilst,jfst,jlst

      integer ios, ios2, ioall
      logical lok, ex, lflag ,recflag
      character*3 char3
      real, external :: spherdist

C
C**********
C*
C 1)  CREATE A HYCOM RELAXATION MASK
C
C 2)  PARAMETERS:
C
C     MODEL GRID SPECIFICATION (W.R.T. PRESSURE GRID):
C
C        IDM    = 1ST DIMENSION OF MAJOR (HYCOM) MODEL ARRAYS
C        JDM    = 2ND DIMENSION OF MAJOR (HYCOM) MODEL ARRAYS
C
C 3)  NAMELIST PARAMETERS
C
C     /MASK/
C        CTITLE      - 2ND LINE OF OUTPUT PREAMBL
C                       (1ST LINE IS ALWAYS 'Relaxation Mask')
C        IF,IL,JF,JL - ARRAY BOX WHERE EFOLD RELAXATION IS APPLIED
C                       = 4*0.0; END OF BOX LIST
C        EFOLD       - RELAXATION E-FOLDING TIME IN DAYS
C                       = 0.0; NO RELAXATION
C
C     THE ALLOWED RANGE FOR IF,IL IS 1 TO IDM-1
C     THE ALLOWED RANGE FOR JF,JL IS 1 TO JDM-1
C
C 4)  INPUT:
C        ON UNIT  5: NAMELIST INPUT
C        ON UNIT 51: BATHYMETRY FILE
C     OUTPUT:
C        ON UNIT 21: RELAXATION MASK
C
C 5)  ALAN J. WALLCRAFT,  NAVAL RESEARCH LABORATORY, JUNE 2000.
C*
C**********
C
      INTEGER I,ISF,ISL,ISEC,J,K
      REAL*4  RVAL,RVMIN,rvmax,rvscale
C
      CHARACTER*1 C(-1:9)
      DATA C / '*', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' /
C
CKAL  --- Efold time and number of boundary cells
      logical,allocatable :: block(:,:)
      integer :: nbnd=7
      integer :: nbnd_pnt=3
      real :: efoldtime=20.0
      logical :: least=.true.
      logical :: lwest=.true.
      logical :: lsouth=.true.
      logical :: lnorth=.true.
CKAL  --- statement function for boundary relaxation
      real ::  rvalxx,bndrmu,jxx
      integer :: nbndxx
      bndrmu(jxx,nbndxx,rvalxx) = 
     &   rvalxx*(nbndxx-jxx)/float(nbndxx-1)
C     bndrmu(jxx,nbndxx,rvalxx) = 1.0
CKAL  --- The above is a linear function for rmu
C
      CALL XCSPMD
      ALLOCATE(   RMU(IDM,JDM) )
      ALLOCATE(  IRMU(IDM,JDM) )
      ALLOCATE( DEPTH(IDM,JDM) )
      ALLOCATE( BLOCK(IDM,JDM) )
C
C     NAMELIST INPUT.
C
CKAL  --- Namelist not used in simplified version
CKAL  DO K= 1,999
CKAL    IF(   K) = 0
CKAL    IL(   K) = 0
CKAL    JF(   K) = 0
CKAL    JL(   K) = 0
CKAL    EFOLD(K) = 0.0
CKAL  ENDDO
CKAL  READ( 5,MASK)
CKAL  WRITE(6,MASK)
CKAL ---- Instead we do this :
      !block=.false.
      inquire(exist=ex,file='rmu.in')
      if (ex) then
         open(567,file='rmu.in')

         ioall = 0.
         ! Read boundary flags (on/off)
         read(567, *,iostat=ios) least  ; ioall = ioall+abs(ios)
         read(567, *,iostat=ios) lwest  ; ioall = ioall+abs(ios)
         read(567, *,iostat=ios) lsouth ; ioall = ioall+abs(ios)
         read(567, *,iostat=ios) lnorth ; ioall = ioall+abs(ios)

         ! Read boundary width
         read(567, *,iostat=ios) nbnd ; ioall = ioall+abs(ios)
         ! Read relaxation time scale
         read(567, *,iostat=ios) rval ; ioall = ioall+abs(ios)
         read(567, *,iostat=ios) nbnd_pnt ; ioall = ioall+abs(ios)

         if (ioall /= 0 ) then
            print *,'Error reading rmu.in'
            stop
         end if

         ! Read Points to use
         ios2=0
         nrpoints=0
         do while (ios2==0) 
            read(567, *,iostat=ios2) char3, pflag, pntlon, pntlat!, prad
            if (ios2==0) then
               lok=.false.
               if (char3=='END') then
                  ios2=-9999
               else
                  nrpoints=nrpoints+1
                  pointname(nrpoints)=char3
                  !pointrad (nrpoints)=prad
                  pointflag(nrpoints)=pflag
                  pointlon (nrpoints)=pntlon
                  pointlat (nrpoints)=pntlat
               end if
            !else
            !   print *,'read error on rmu.in -- make sure '
            !   print *,'you have a working version'
            !   stop 
            end if
         end do
         print *,'extra point read is :',nrpoints,pointname(1),
     &   pointflag(1),pointlon(1),pointlat(1)

         ! Read rectangle overrides
         ios2=0
         do while (ios2==0)
            read(567,*,iostat=ios2) ifst,ilst,jfst,jlst,recflag
            if (ifst<=idm .and. ifst>0 .and.
     &          ilst<=idm .and. ilst>0 .and.
     &          jfst<=jdm .and. jfst>0 .and.
     &          jlst<=jdm .and. jlst>0 .and. ios2==0 ) then
               do i=ifst,ilst
               do j=jfst,jlst
                  block(i,j) = .not.recflag
               end do
               end do
            !else if (ios2==0) then
            !   print *,'Invalid region specified for rectanlge'
            !   print *,'Grid dimension (idm,jdm):',idm,jdm
            !   print *,'ifirst,ilast:',ifst,ilst
            !   print *,'jfirst,jlast:',jfst,jlst
            !   stop
            end if
         end do

         close(567)
      else
         print *,'No rmu.in present'
         stop
      end if

C
C     TOPOGRAPHY INPUT.
C
      CALL ZAIOST
C
      CALL ZHOPNC(51, 'regional.depth.b', 'FORMATTED', 'OLD', 0)
      CALL ZAIOPF('regional.depth.a', 'OLD', 51)
      !CALL ZHOPEN(51, 'FORMATTED', 'OLD', 0)
      !CALL ZAIOPN('OLD', 51)
      READ (51,'(A79)') PREAMBL
      READ (51,'(A)')   CLINE
      CLOSE(UNIT=51)
      WRITE(6,'(/(1X,A79))') PREAMBL,CLINE
C
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
C
      CALL ZAIORD(DEPTH,IRMU,.FALSE., HMINA,HMAXA, 51)
      CALL ZAIOCL(51)
C     Added by KAL  
      WHERE(DEPTH>0.5*2.0**99) DEPTH=0.
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
C     INITIALIZE MASK TO ZERO
C
      DO J= 1,JDM
        DO I= 1,IDM
          RMU(I,J) = 0.0
        ENDDO
      ENDDO
C
C     INITIALIZE SPECIFIED PATCHES.
C
CKAL  RVMIN = 1.E20
CKAL  DO K= 1,999
CKAL    IF     (MIN(IF(K),IL(K),JF(K),JL(K)).GT.0   .AND.
CKAL +                              EFOLD(K).GT.0.0      ) THEN
CKAL      IF     (IF(K).GT.IL(K) .OR.
CKAL +            IF(K).LT.1     .OR.
CKAL +            IL(K).GT.IDM-1     ) THEN
CKAL        WRITE(6,9000) K,IF(K),IL(K),IDM-1
CKAL        CALL ZHFLSH(6)
CKAL        STOP
CKAL      ENDIF
CKAL      IF     (JF(K).GT.JL(K) .OR.
CKAL +            JF(K).LT.1     .OR.
CKAL +            JL(K).GT.JDM-1     ) THEN
CKAL        WRITE(6,9100) K,JF(K),JL(K),JDM-1
CKAL        CALL ZHFLSH(6)
CKAL        STOP
CKAL      ENDIF
C
CKAL      RVAL = 1.0/(EFOLD(K)*86400.0)
CKAL      DO J= JF(K),JL(K)
CKAL        DO I= IF(K),IL(K)
CKAL          IF     (DEPTH(I,J).NE.SPVAL) THEN
CKAL            RMU(I,J) = MAX(RMU(I,J),RVAL)
CKAL          ENDIF
CKAL        ENDDO
CKAL      ENDDO
CKAL      RVMIN = MIN( RVMIN, RVAL )
CKAL    ENDIF
CKAL  ENDDO
C
C
C
C
CKAL  --- A simpler masking method; Applied on boundaries
CKAL  --- where depth 2 points away from border is > 0
CKAL  --- e.g. depth(2,j) > 0 and depth(itdm-1) > 0
      rval=1.0/(efoldtime*86400.0)
      rvmin=1e20
      do i=2,idm-1

         if (lsouth) then
         if (depth(i,2)>0.5.and..not.block(i,1)) then
            do j=2,nbnd
               if (depth(i,j)>0.5) then
                  rmu(i,j)=bndrmu(float(j-1),nbnd,rval)
               end if
               !print *,'east:',i,j,j,nbnd,rmu(i,j)
            end do
         end if
         end if

         if (lnorth) then
         if (depth(i,jdm-1)>0.5.and..not.block(i,jdm)) then
            do j=2,nbnd
               if (depth(i,jdm-j+1)>0.5) then
                  rmu(i,jdm-j+1)=max(bndrmu(float(j-1),nbnd,rval),
     &                           rmu(i,jdm-j+1))
               end if
               !print *,'west:',i,j,j,nbnd, rmu(i,jdm-j+1)
            end do
         end if
         end if

      end do

      do j=2,jdm-1

         if (lwest) then
         if (depth(2,j)>0.5.and..not.block(1,j)) then
            do i=2,nbnd
               if (depth(i,j)>0.5) then
                  rmu(i,j)=max(bndrmu(float(i-1),nbnd,rval),rmu(i,j))
               end if
            !print *,'south:',2,j,i,nbnd,rmu(i,j)
            end do
         end if
         end if

         if (least) then
         if (depth(idm-1,j)>0.5.and..not.block(idm,j)) then
            do i=2,nbnd
               if (depth(idm-i+1,j)>0.5) then
                  rmu(idm-i+1,j)=max(bndrmu(float(i-1),nbnd,rval),
     &                               rmu(idm-i+1,j))
               end if
               !print *,'north:',idm-1,j,i,nbnd, rmu(idm-i+1,j)
            end do
         end if
         end if
      end do
CKAL  --- End of simple mask...
CKAL  --- Start more complex mask: -- Mediterranean inflow
C
      ALLOCATE(PLON(IDM,JDM))
      ALLOCATE(PLAT(IDM,JDM))
      CALL ZHOPNC(31, 'regional.grid.b', 'FORMATTED', 'OLD', 0)
      CALL ZAIOPF('regional.grid.a', 'OLD', 31)
C
      READ(31,*) ! skip idm
      READ(31,*) ! skip jdm
      READ(31,*) ! skip mapflg
C     READ(31,*) ! skip plon
C     CALL ZAIOSK(31)
C
      READ(31,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
      CALL ZAIORD(PLON,IRMU,.FALSE., HMINA,HMAXA, 31)
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
      READ(31,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
      CALL ZAIORD(PLAT,IRMU,.FALSE., HMINA,HMAXA, 31)
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (plat):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF

C --- Find closest point to med (spherdist)
      do ipoint = 1,nrpoints

         if (pointflag(ipoint)) then
         !rmaxx=1e20
         rmaxx_i = -200
         rmaxx_j = -200
         rmaxx=600000 ! 600 km

         print *, 'processing '//pointname(ipoint)
         DO j=1,JDM
         DO i=1,IDM
            rdist=spherdist(pointlon(ipoint),pointlat(ipoint),
     &                      plon(i,j),plat(i,j))
            if (rdist<rmaxx) then
               rmaxx=rdist
               rmaxx_i=i
               rmaxx_j=j
            end if
         end do
         end do

C --- With this point as center, draw out rmu
         do i=1,idm
         do j=1,jdm
            celldist=sqrt(float(i-rmaxx_i)**2+float(j-rmaxx_j)**2)
            celldist=max(celldist,2.)
            if (celldist<nbnd_pnt.and.depth(i,j)>0.5) then
               rmu(i,j)=max(rmu(i,j),bndrmu(celldist-1,nbnd_pnt,rval))
            end if
         end do
         end do
         end if
      end do

C
CKAL  WRITE(6,*)
CKAL  DO K= 1,999
CKAL    IF     (MIN(IF(K),IL(K),JF(K),JL(K)).GT.0   .AND.
CKAL +                              EFOLD(K).GT.0.0      ) THEN
CKAL      RVAL = 1.0/(EFOLD(K)*86400.0)
CKAL      WRITE(6,6200) IF(K),IL(K),JF(K),JL(K),
CKAL +                  EFOLD(K),RVAL,RVAL/RVMIN
CKAL      CALL ZHFLSH(6)
CKAL    ENDIF
CKAL  ENDDO
C
C     WRITE OUT THE MASK
C
      CTITLE = 'Relaxation mask created with program rmunew'
      PREAMBL(1) = 'Relaxation Mask'
      PREAMBL(2) = CTITLE
      PREAMBL(3) = ' '
      PREAMBL(4) = ' '
      WRITE(PREAMBL(5),'(A,2I5)')
     +        'i/jdm =',
     +       IDM,JDM
C
      CALL ZAIOPN('NEW', 21)
      CALL ZAIOWR(RMU,IRMU,.FALSE., RMUMIN,RMUMAX, 21, .FALSE.)
      CALL ZAIOCL(21)
C
      CALL ZHOPEN(21, 'FORMATTED', 'NEW', 0)
      WRITE(21,4101) PREAMBL
      WRITE(21,4102) '     rmu',RMUMIN,RMUMAX
      CLOSE(UNIT=21)
C
      WRITE(6, *)
      WRITE(6, 4101) PREAMBL
      WRITE(6, 4102) ' rmu',RMUMIN,RMUMAX
      WRITE(6, *)
      CALL ZHFLSH(6)
C
C     PRINTOUT THE MASK
C
      rvmin=minval(rmu,rmu>0.)
      rvmax=maxval(rmu,rmu>0.)
      rvscale=(rvmax-rvmin)
      WRITE(6,6000) IDM,JDM,RVMIN
      ISEC = (IDM-1)/100 + 1
      DO K= 1,ISEC
        ISF = (K-1)*100 + 1
        ISL = MIN(IDM, ISF+100-1)
        WRITE(6,6050) ISF,ISL
        DO J= JDM,1,-1
          DO I= ISF,ISL
C           IF     (DEPTH(I,J).EQ.SPVAL) THEN
            IF     (DEPTH(I,J) <0.5) THEN
              IRMU(I,J) = -1
            ELSEIF (RMU(I,J).EQ.0.0) THEN
              IRMU(I,J) =  0
            ELSE
CKAL          IRMU(I,J) =  MIN( 9, NINT( RMU(I,J)/RVMIN ) )
              IRMU(I,J)=MIN(9,NINT(1+8*(RMU(I,J)-RVMIN)/rvscale))
            ENDIF
          ENDDO
          WRITE(6,6100) J,(C(IRMU(I,J)),I=ISF,ISL)
        ENDDO
      ENDDO
      STOP
C
 4101 FORMAT(A79)
 4102 FORMAT(A,': range = ',1P2E16.7)
 6000 FORMAT(1H1 / 30X,'RELAXATION MASK FOR AN',I5,'  BY',I5,' MESH.'
     +           / 31X,'LAND = *,  OCEAN = NINT(RMU /',1PE9.2,')' )
 6050 FORMAT(/ / / 21X,'I =',I5,'  TO',I5,'  :' / /)
 6100 FORMAT(4X,'J =',I5,5X,10(10A1,1X))
 6200 FORMAT('IF,IL,JF,JL =',4I5,'   EFOLD,RVAL =',F6.1,1PE12.3,
     +   '  (',0PF4.1,')')
 9000 FORMAT('ERROR - ILLEGAL IF OR IL FOR K =',I3,
     +       '   MUST HAVE 1<=IF(K)<=IL(K)<=IDM-1' /
     +       'IF(K),IL(K),IDM-1 = ',3I5 /)
 9100 FORMAT('ERROR - ILLEGAL JF OR JL FOR K =',I3,
     +       '   MUST HAVE 1<=JF(K)<=JL(K)<=JDM-1' /
     +       'JF(K),JL(K),JDM-1 = ',3I5 /)
      END

      PROGRAM HYCOM_ARCTIC
      IMPLICIT NONE
C
C  hycom_arctic - Usage:  hycom_arctic   fin.a idm jdm [jdma] type     fout.a
C                         hycom_arctic_g fin.a idm jdm [jdma] type.txt fout.a
C
C                 Outputs a valid arctic patch version of each input field.
C
C  f.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  jdma is the minimum extent for a full arctic, default jdm
C
C  type is the grid type:
C   scalar p=1,q=2,u=3,v=4; vector p=11,u=13,v=14; angle p=21.
C  vector always implies a sign change across the patch mid-point.
C
C  for hycom_arctic_g, type.txt contains a list of the grid types with
C  one type per line and one line per input record
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  September 2001.
C
      REAL*4, ALLOCATABLE :: A(:,:)
      REAL*4              :: PAD(4096)
C
      INTEGER       IOS,L
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      LOGICAL       LTYPE
      INTEGER       IDM,JDM,JDMA,ITYPE,NPAD
      CHARACTER*240 CFILE,CFILET,CFILEO
C
C     READ ARGUMENTS.
C
      CALL GETARG(0,CARG)
      L = LEN_TRIM(CARG)
*     WRITE(6,"(4a)") TRIM(CARG),'"',CARG(L-4:L),'"'
      IF     (CARG(L-1:L).EQ.'_g') THEN
        LTYPE=.TRUE.
      ELSEIF (CARG(L-6:L).EQ.'_arctic') THEN
        LTYPE=.FALSE.
      ELSE
        WRITE(6,'(2a)')
     &    'Usage:  ',
     &    'hycom_arctic or hycom_arctic_g'
        CALL EXIT(1)
      ENDIF
C
      NARG = IARGC()
C
      IF     (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,5000) IDM
        CALL GETARG(3,CARG)
        READ(CARG,5000) JDM
        JDMA   = JDM
        IF     (LTYPE) THEN
          CALL GETARG(4,CFILET)
        ELSE
          CALL GETARG(4,CARG)
          READ(CARG,5000) ITYPE
        ENDIF
        CALL GETARG(5,CFILEO)
      ELSEIF (NARG.EQ.6) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,5000) IDM
        CALL GETARG(3,CARG)
        READ(CARG,5000) JDM
        CALL GETARG(4,CARG)
        READ(CARG,5000) JDMA
        IF     (LTYPE) THEN
          CALL GETARG(5,CFILET)
        ELSE
          CALL GETARG(5,CARG)
          READ(CARG,5000) ITYPE
        ENDIF
        CALL GETARG(6,CFILEO)
      ELSEIF (LTYPE) THEN
        WRITE(6,*) 'Usage: ' //
     +   'hycom_arctic_g fin.a idm jdm [jdma] type.txt fout.a'
        CALL EXIT(1)
      ELSE
        WRITE(6,*) 'Usage: ' //
     +   'hycom_arctic fin.a idm jdm [jdma] type fout.a'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_arctic: could not allocate ',
     +             IDM*JDM,' words for A'
        CALL EXIT(2)
      ENDIF
C
      CALL ARCTIC(A,IDM,JDM,PAD,NPAD, JDMA,ITYPE,LTYPE,
     &            CFILE,CFILET,CFILEO)
      CALL EXIT(0)
 5000 FORMAT(I4)
      END
      SUBROUTINE ARCTIC(A,IDM,JDM, PAD,NPAD,
     +                    JDMA,ITYPE,LTYPE, CFILE,CFILET,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFILET,CFILEO
      LOGICAL       LTYPE
      INTEGER       IDM,JDM,NPAD,JDMA,ITYPE
      REAL*4        A(IDM,JDM),PAD(NPAD)
      REAL*4        PI,AMX,AMN
C
C     MOST OF WORK IS DONE HERE.
C
      CHARACTER*18 CASN
      INTEGER      I,II,J,JJ,K,IOS,NRECL
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
      INQUIRE( IOLENGTH=NRECL) A,PAD
#ifdef CRAY
#ifdef t3e
      IF     (MOD(NRECL,4096).EQ.0) THEN
        WRITE(CASN,8000) NRECL/4096
 8000   FORMAT('-F cachea:',I4.4,':1:0')
        IU8 = 12
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit ',IU8
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
        IU8 = 21                        
        CALL ASNUNIT(IU8,CASN,IOS8)     
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit ',IU8
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
      ENDIF
#else
      CALL ASNUNIT(12,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 12'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
      CALL ASNUNIT(21,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 21'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
#endif
#endif
      OPEN(UNIT=12, FILE=CFILE, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=21, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEO)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      IF     (LTYPE) THEN
        OPEN(UNIT=99, FILE=CFILET, FORM='FORMATTED', STATUS='OLD',
     +           IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'Error: can''t open ',TRIM(CFILET)
          WRITE(6,*) 'ios   = ',ios
          CALL EXIT(3)
        ENDIF
      ENDIF
C
      DO K= 1,HUGE(K)
        READ(12,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILE)
            CALL EXIT(4)
          ELSE  ! end-of-file
            EXIT  !k-loop
          ENDIF
        ENDIF
        IF     (LTYPE) THEN
          READ(99,*,IOSTAT=IOS) ITYPE
          IF     (IOS.NE.0) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILET)
            CALL EXIT(4)
          ENDIF
        ENDIF
        IF     (ITYPE.EQ.1) THEN !p-grid,scalar
          DO J= JDMA,JDM
          JJ = JDMA-1-(J-JDMA)
          DO I= 1,IDM
            II = IDM-MOD(I-1,IDM)
            A(I,J) = A(II,JJ)
          ENDDO !i
          ENDDO !j
        ELSEIF (ITYPE.EQ.11) THEN !p-grid,vector
          DO J= JDMA,JDM
          JJ = JDMA-1-(J-JDMA)
          DO I= 1,IDM
            IF     (A(I,J).NE.SPVAL) THEN
              II = IDM-MOD(I-1,IDM)
              A(I,J) = -A(II,JJ)
            ENDIF
          ENDDO !i
          ENDDO !j
        ELSEIF (ITYPE.EQ.21) THEN !p-grid,pang from regional.grid
          PI = 4.D0*ATAN(1.D0)
          DO J= JDMA,JDM
          JJ = JDMA-1-(J-JDMA)
          DO I= 1,IDM
            II = IDM-MOD(I-1,IDM)
            IF     (A(II,JJ).GE.0.0) THEN
              A(I,J) = PI - A(II,JJ)
            ELSE
              A(I,J) = PI + A(II,JJ)
            ENDIF
          ENDDO !i
          ENDDO !j
        ELSEIF (ITYPE.EQ.2) THEN !q-grid,scalar
          DO J= JDMA,JDM
          JJ = JDMA-(J-JDMA)
          DO I= 1,IDM
            II = MOD(IDM-(I-1),IDM)+1
            A(I,J) = A(II,JJ)
          ENDDO !i
          ENDDO !j
        ELSEIF (ITYPE.EQ.3) THEN !u-grid,scalar
          DO J= JDMA,JDM
          JJ = JDMA-1-(J-JDMA)
          DO I= 1,IDM
            II = MOD(IDM-(I-1),IDM)+1
            A(I,J) = A(II,JJ)
          ENDDO !i
          ENDDO !j
        ELSEIF (ITYPE.EQ.13) THEN !u-grid,vector
          DO J= JDMA,JDM
          JJ = JDMA-1-(J-JDMA)
          DO I= 1,IDM
            IF     (A(I,J).NE.SPVAL) THEN
              II = MOD(IDM-(I-1),IDM)+1
              A(I,J) = -A(II,JJ)
            ENDIF
          ENDDO !i
          ENDDO !j
        ELSEIF (ITYPE.EQ.4) THEN !v-grid,scalar
          DO J= JDMA,JDM
          JJ = JDMA-(J-JDMA)
          DO I= 1,IDM
            II = IDM-MOD(I-1,IDM)
            A(I,J) = A(II,JJ)
          ENDDO !i
          ENDDO !j
        ELSEIF (ITYPE.EQ.14) THEN !v-grid,vector
          J  = JDMA
          JJ = JDMA
          DO I= 1,IDM/2
            IF     (A(I,J).NE.SPVAL) THEN
              II = IDM-MOD(I-1,IDM)
              A(I,J) = -A(II,JJ)
            ENDIF
          ENDDO !i
          DO J= JDMA+1,JDM
          JJ = JDMA-(J-JDMA)
          DO I= 1,IDM
            IF     (A(I,J).NE.SPVAL) THEN
              II = IDM-MOD(I-1,IDM)
              A(I,J) = -A(II,JJ)
            ENDIF
          ENDDO !i
          ENDDO !j
        ELSE
          WRITE(6,*) 'Unknown type (not 1,2,3,4,11,13,14)'
          WRITE(6,*) 'Usage: ' //
     +     'hycom_arctic f.a idm jdm [jdma] type fout.a'
          CALL EXIT(2)
        ENDIF !itype
C
        AMN =  SPVAL
        AMX = -SPVAL
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (A(I,J).NE.SPVAL) THEN
              AMX = MAX( AMX, A(I,J) )
              AMN = MIN( AMN, A(I,J) )
            ENDIF
          ENDDO !i
        ENDDO !j
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif  
        WRITE(21,REC=K,IOSTAT=IOS) A
        WRITE(6,'(a,i2,1p2g16.8)')
     &     'itype,min,max = ',ITYPE,AMN,AMX
      ENDDO !k
      WRITE(6,*)
      WRITE(6,*) K-1,' FIELDS PROCESSED'
      WRITE(6,*)
C     
      CLOSE(11)
      CLOSE(21)
      IF     (LTYPE) THEN
        CLOSE(99)
      ENDIF
      CALL EXIT(0)  ! normal exit
      RETURN
      END

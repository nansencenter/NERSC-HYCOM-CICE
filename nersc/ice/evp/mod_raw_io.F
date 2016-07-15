      module mod_raw_io
      contains


! Modified from Alan Wallcraft's RAW routine by Knut Liseter @ NERSC
! So far only the "I" in "IO" is present
      SUBROUTINE READRAW(A,AMN,AMX,IDM,JDM,LSPVAL,SPVAL,CFILE1,K)
      IMPLICIT NONE
C
      REAL*4     SPVALH
      PARAMETER (SPVALH=1e30)
C
      REAL*4,        INTENT(OUT) :: A(IDM,JDM)
      REAL*4,        INTENT(OUT) :: AMN,AMX
      INTEGER,       INTENT(IN)  :: IDM,JDM
      LOGICAL,       INTENT(IN)  :: LSPVAL
      REAL*4,        INTENT(INOUT)  :: SPVAL
      INTEGER,       INTENT(IN)  :: K
      CHARACTER(len=*), INTENT(IN)  :: CFILE1
C
      REAL*4 :: PADA(4096)
C
C     MOST OF WORK IS DONE HERE.
C

      CHARACTER*18 CASN
      INTEGER      LEN_TRIM
      INTEGER      I,J,IOS,NRECL
      INTEGER NPAD
C
      IF(.NOT.LSPVAL) THEN
        SPVAL = SPVALH
      ENDIF
C
!!! Calculate the number of elements padded!!!!!!!!!!!!!!!!!!!!!!!!
      NPAD=GET_NPAD(IDM,JDM)
C
      INQUIRE( IOLENGTH=NRECL) A,PADA(1:NPAD)
C     
C     
      OPEN(UNIT=11, FILE=CFILE1, FORM='UNFORMATTED', STATUS='old',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILE1(1:LEN_TRIM(CFILE1))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      READ(11,REC=K,IOSTAT=IOS) A
      close(11)
C
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read record ',K,
     &             ' from '//CFILE1(1:LEN_TRIM(CFILE1))
        CALL EXIT(4)
      ENDIF
C
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
C                 
      RETURN
      END SUBROUTINE





      INTEGER FUNCTION GET_NPAD(IDM,JDM)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IDM,JDM
         GET_NPAD = 4096 - MOD(IDM*JDM,4096)
         GET_NPAD = mod(GET_NPAD,4096)
      END FUNCTION
      end module mod_raw_io

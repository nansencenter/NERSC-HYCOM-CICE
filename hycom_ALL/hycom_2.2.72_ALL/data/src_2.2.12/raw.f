       SUBROUTINE RAW(A,IDM,JDM,PAD,NPAD,
     &  LSPVAL,SPVAL, CFILE1,k)
      IMPLICIT NONE
C
      REAL*4     SPVALH
      PARAMETER (SPVALH=2.0**100)
C
      CHARACTER*240 :: CFILE1
      
      LOGICAL      LSPVAL
      INTEGER      IDM,JDM,NPAD
      REAL*4       SPVAL
      REAL*4       A(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C

      CHARACTER*18 CASN
      INTEGER      LEN_TRIM
      INTEGER      I,J,K,IOS,NRECL,MRECL
      REAL*4       AMN,AMX
C
      IF(.NOT.LSPVAL) THEN
        SPVAL = SPVALH
      ENDIF
C
      INQUIRE( IOLENGTH=MRECL) A
      INQUIRE( IOLENGTH=NRECL) A,PAD
      
      
      OPEN(UNIT=11, FILE=CFILE1, FORM='UNFORMATTED', STATUS='unknown',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILE1(1:LEN_TRIM(CFILE1))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF

        READ(11,REC=K,IOSTAT=IOS) A
        AMN =  SPVALH
        AMX = -SPVALH
        DO 210 J= 1,JDM
          DO 212 I=1,IDM

            IF     (A(I,J).LE.SPVALH) THEN
              AMN = MIN( AMN, A(I,J) )
              AMX = MAX( AMX, A(I,J) )
	      
            ELSEIF (LSPVAL) THEN
              A(I,J) = SPVAL
            ENDIF
                  
  212     CONTINUE
  210   CONTINUE
      RETURN
      close(11)
      END

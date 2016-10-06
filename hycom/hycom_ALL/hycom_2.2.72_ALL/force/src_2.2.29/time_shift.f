      PROGRAM TIME_SHIFT
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     WIND/FLUX ARRAYS.
C
      CHARACTER PREAMBL(5)*79,PREAMBL2(5)*79,CNAME*8
C
C     NAMELIST.
C
      REAL*8           FINC,FSTART,WSTART,TSTART,TMAX
      NAMELIST/AFTIME/ FINC,FSTART,WSTART,TSTART,TMAX
      INTEGER          INTERP,ILIMIT,ICOMBI,NYEAR
      NAMELIST/AFFLAG/ INTERP,ILIMIT,ICOMBI,NYEAR
C
C**********
C*
C 1)  FROM ACTUAL-YEAR FLUX FIELD .b FILE
C      SHIFT TO A DIFFERENT START DATE.
C
C 2)  NO 2.
C
C 3)  NAMELIST INPUT:
C
C     /AFTIME/
C        FSTART - TIME OF OUTPUT HEAT FLUX START       (DAYS)
C        WSTART - TIME OF OUTPUT WIND START            (DAYS)
C        TSTART - TIME OF INPUT  HEAT FLUX/WIND START  (DAYS)
C
C     NAMELIST /AFTIME/ IS PATTERNED AFTER /XXTIME/ SO THAT THE
C      MODEL''S STANDARD AWK-BASED RUN SCRIPT CUSTOMIZER CAN ALSO
C      WORK FOR THE WIND/FLUX GENERATION SCRIPT.  IN PARTICULAR, 
C      'WSTART' AND 'FSTART' MUST BE EQUAL SINCE THIS PROGRAM
C      WORKS FOR BOTH WINDS AND FLUXES.
C
C 4)  INPUT:
C        ON UNIT  5:    NAMELIST /AFTIME/
C        ON UNIT 20:    FORMATTED MODEL FIELD FILE TO SHIFT
C     OUTPUT:
C        ON UNIT 10:    FORMATTED MODEL FIELD FILE
C*
C**********
C
      CHARACTER*240 CLINE
C
      INTEGER       IOS,K
      REAL*8        OFFSET,WDAY
C
C     NAMELIST INPUT.
C
      CALL ZHOPEN(6, 'FORMATTED', 'UNKNOWN', 0)
C
      FSTART = 0.0
      WSTART = 0.0
      TSTART = 0.0
      TMAX   = 0.0
      WRITE(6,*) 'READING /AFTIME/'
      CALL ZHFLSH(6)
      READ( 5,AFTIME)
      WRITE(6,AFTIME)
      WRITE(6,*) 
      CALL ZHFLSH(6)
C
      OFFSET = TSTART - WSTART
      WRITE(6,*) 'OFFSET = ',OFFSET
      WRITE(6,*)
C
C     INITIALIZE OUTPUT.
C
      CALL ZHOPEN(10, 'FORMATTED', 'NEW', 0)
      CALL ZHOPEN(20, 'FORMATTED', 'OLD', 0)
      READ( 20,'(A79)') PREAMBL
      WRITE(PREAMBL(4),'(a,f12.5,a,f12.5)') 
     &  'offset from',TSTART,' to',WSTART
      WRITE(10,'(A79)') PREAMBL
      WRITE(6,*)
      WRITE(6, '(A79)') PREAMBL
      WRITE(6,*)
C
      DO K= 1,HUGE(K)
        READ( 20,          '(A)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
        READ( CLINE(27:38),'(F12.5)')        WDAY
        WDAY = WDAY - OFFSET
        WRITE(CLINE(27:38),'(F12.5)')        WDAY
        WRITE(10,          '(A)')       TRIM(CLINE)
        WRITE( 6,          '(A)')       TRIM(CLINE)
        CALL ZHFLSH(6)
      ENDDO
C
      CLOSE( UNIT=10)
C
C     SUMMARY.
C
      WRITE(6,*) K,' RECORDS WRITTEN'
      CALL ZHFLSH(6)
      STOP
      END
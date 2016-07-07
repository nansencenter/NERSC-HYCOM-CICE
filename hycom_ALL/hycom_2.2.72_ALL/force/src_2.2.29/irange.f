      PROGRAM I_RANGE
      IMPLICIT NONE
C
      INTEGER   FLG_SCL
      REAL*4    F4MIN,F4MAX,R4MIN,R4MAX
      REAL*8    FLDRNG,SCLRNG
      REAL*4    SCALE_I,ADD_OFF,SCL_FAC
      INTEGER*2 I2MIN,I2MAX
C
      READ( 5,*)        F4MIN,F4MAX,FLG_SCL
      WRITE(6,*) "IN: ",F4MIN,F4MAX,FLG_SCL
      WRITE(6,*) 
C
          IF     (F4MIN.EQ.F4MAX) THEN
            SCALE_I = 1
            ADD_OFF = F4MIN
          ELSEIF (FLG_SCL.EQ.1) THEN
C
C           FIT TO TOTAL RANGE
C
            FLDRNG = F4MAX - F4MIN
            CALL IRANGE(SCALE_I, FLDRNG,SCLRNG)
            ADD_OFF = 0.5*(F4MIN + F4MAX)
          ELSEIF (FLG_SCL.EQ.2) THEN
C
C           FIT CENTERED ON ZERO
C
            FLDRNG = 2.D0*MAX(ABS(F4MAX),ABS(F4MIN))
            CALL IRANGE(SCALE_I, FLDRNG,SCLRNG)
            ADD_OFF = 0.0
          ELSEIF (FLG_SCL.EQ.3) THEN
C
C           FIT BETWEEN 0.0 AND MAXIMUM
C
            FLDRNG = F4MAX
            CALL IRANGE(SCALE_I, FLDRNG,SCLRNG)
            ADD_OFF = 0.5d0*SCLRNG
          ELSEIF (FLG_SCL.EQ.4) THEN
C
C           FIT BETWEEN MINIMUM AND 0.0
C
            FLDRNG = -F4MIN
            CALL IRANGE(SCALE_I, FLDRNG,SCLRNG)
            ADD_OFF = -0.5d0*SCLRNG + 1.D0/SCALE_I
          ENDIF
        SCL_FAC = 1.0/SCALE_I
C
        WRITE(6,*) "SCALE_I = ",SCALE_I
        WRITE(6,*) "SCL_FAC = ",SCL_FAC
        WRITE(6,*) "ADD_OFF = ",ADD_OFF
        WRITE(6,*) 
        I2MIN = -2**15
        I2MAX =  2**15-1
        R4MIN = I2MIN * SCL_FAC + ADD_OFF
        R4MAX = I2MAX * SCL_FAC + ADD_OFF
        WRITE(6,'(1x,a,e20.11,a,e20.11)')
     &   '  PACKED  MIN=',R4MIN,' MAX=',R4MAX
        WRITE(6,'(1x,a,f20.14,a,f20.14)')
     &   '  PACKED  MIN=',R4MIN,' MAX=',R4MAX
C
        WRITE(6,'(1x,a,f20.14,a,f20.14)')
     &     ' OVERALL MIN=',F4MIN,' MAX=',F4MAX
        I2MIN = MAX( -2**15, MIN( 2**15-1,
     &                     NINT((F4MIN-ADD_OFF)*SCALE_I) ))
        I2MAX = MAX( -2**15, MIN( 2**15-1,
     &                     NINT((F4MAX-ADD_OFF)*SCALE_I) ))
        F4MIN = I2MIN * SCL_FAC + ADD_OFF
        F4MAX = I2MAX * SCL_FAC + ADD_OFF
        WRITE(6,'(1x,a,f20.14,a,f20.14)')
     &     'UNPACKED MIN=',F4MIN,' MAX=',F4MAX
        WRITE(6,'(1x,a,i20,  a,i20 /)')
     &     '  PACKED MIN=',I2MIN,' MAX=',I2MAX
      END
      SUBROUTINE IRANGE(SCALE_I, FLDRNG,SCLRNG)
      IMPLICIT NONE
C
      REAL             SCALE_I
      DOUBLE PRECISION FLDRNG,SCLRNG
C
C     CALCULTE A GOOD SCALE_I FOR FLDRNG
C     SCLRNG IS RETURNED AS THE ACTUAL RANGE WITH SCALE_I
C
      INTEGER*8 I_SCALE
      INTEGER   IP
C
      IF     (FLDRNG.GT.2.D0**13) THEN
        WRITE(6,'(/ a,e16.6 /)')
     &   'error in IRANGE: FLDRNG too big = ',FLDRNG
        STOP
      ENDIF
C
      IF     (FLDRNG.LT.2.D0**-30) THEN
        WRITE(6,'(/ a,e16.6 /)')
     &   'error in IRANGE: FLDRNG too small = ',FLDRNG
        STOP
      ENDIF
C
      I_SCALE = (2.D0**16-1.D0) / FLDRNG
      write(6,*) 'fldrng  = ',fldrng
      write(6,*) 'i_scale = ',i_scale
      DO IP= 1,62
        IF     (I_SCALE.LT.2.D0**IP) THEN
          EXIT
        ENDIF
      ENDDO !ip
      write(6,*) 'ip      = ',ip     
      IF     (IP.GT.2 .AND. I_SCALE.LT.3.D0*2.D0**(IP-2)) THEN
        I_SCALE =      2.D0**(IP-1)
      ELSE
        I_SCALE = 3.D0*2.D0**(IP-2)
      ENDIF
      write(6,*) 'i_scale = ',i_scale
      SCALE_I = I_SCALE
      SCLRNG  = I_SCALE
      SCLRNG  = 2d0**16/SCLRNG  !overestimate by 1/I_SCALE
      write(6,*) 'fldrng  = ',fldrng
      write(6,*) 'sclrng  = ',sclrng 
C
C     END OF IRANGE.
      END

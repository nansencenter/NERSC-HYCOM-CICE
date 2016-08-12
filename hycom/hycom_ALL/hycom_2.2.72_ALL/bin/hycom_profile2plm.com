      PROGRAM HYCOM_PROFILE2PLM
      IMPLICIT NONE
C
C  hycom_profile2pcm - Usage: hycom_profile2pcm archv.txt archp.txt
C
C                 converts a HYCOM text profile file to one
C                 which will plot as the PLM profile
C
C   archv.txt is assumed to be an HYCOM archive text profile file
C   archp.txt will be the output, PLM plot, text profile file
C
C  the input and the output are equivalent HYCOM profiles, but
C  the ouput adds two zero thickness layers at each interface
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  Auguast 2005.
C
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      CHARACTER*240 CFILEA,CFILEC
      CHARACTER*240 CLINE
      REAL          THK,FLAG
      INTEGER       IOS,K,KDM
C
      REAL, ALLOCATABLE :: UVTSR(:,:),P(:)
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.2) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEC)
      ELSE
        WRITE(6,*)
     +    'Usage: hycom_profile2pcm archv.txt archz.txt'
        CALL EXIT(1)
      ENDIF
C
C     OPEN ALL FILES.
C
      OPEN(UNIT=11, FILE=CFILEA, FORM='FORMATTED', STATUS='OLD',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',CFILEA(1:LEN_TRIM(CFILEA))
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=21, FILE=CFILEC, FORM='FORMATTED', STATUS='NEW',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',CFILEC(1:LEN_TRIM(CFILEC))
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(5)
      ENDIF
C
C     COPY PROFILE HEADER TO OUTPUT.
C
      DO K= 1,4
        READ( 11,'(a)') CLINE
        WRITE(21,'(a)') CLINE(1:LEN_TRIM(CLINE))
      ENDDO
      READ( 11,'(a)') CLINE
C
C     READ THE ISOPYCNAL PROFILE.
C
      P(1) =  0.0
      KDM  = -1
      DO K= 1,99
        READ(11,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          IF     (K.NE.KDM+1) THEN
            WRITE(6,*) 'Error: inconsistent input profile'
            CALL EXIT(6)
          ENDIF
          EXIT
        ENDIF
        READ(CLINE,*) KDM,U(K),V(K),T(K),S(K),R(K),THK
        P(K+1) = P(K) + THK
        IF     (THK.EQ.0.0) THEN
          U(K) = U(K-1)
          V(K) = V(K-1)
          T(K) = T(K-1)
          S(K) = S(K-1)
          R(K) = R(K-1)
        ENDIF
      ENDDO
      CLOSE(11)
C
C     OUTPUT, 3 COPIES OF EACH LAYER.
C
      WRITE(21,'(a,a)')
     &  '#   k',
     &  '    utot    vtot    temp    saln    dens    thkns      dpth'
      DO K= 1,KDM
        THK = P(K+1) - P(K)
        WRITE(21,'(i4,1x,2f8.2,3f8.3,f9.3,f10.3)')
     &    3*K-2,U(K),V(K),T(K),S(K),R(K),0.0,P(K)
        WRITE(21,'(i4,1x,2f8.2,3f8.3,f9.3,f10.3)')
     &    3*K-1,U(K),V(K),T(K),S(K),R(K),THK,0.5*(P(K)+P(K+1))
        WRITE(21,'(i4,1x,2f8.2,3f8.3,f9.3,f10.3)')
     &    3*K  ,U(K),V(K),T(K),S(K),R(K),0.0,P(K+1)
      ENDDO
      END

      PROGRAM HYCOM_PROFILE_LOCSIG
      IMPLICIT NONE
C
C  hycom_profile_locsig - Usage:  hycom_profile_locsig archv.txt archs.txt
C
C                 replace density with
C                 locally referenced potential density (sigloc)
C                 in a HYCOM isopycnal text profile file
C
C   archv.txt is assumed to be an HYCOM archive text profile file
C   archs.txt will be the output text profile file, with sigloc
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  December 2006.
C
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      CHARACTER*240 CFILEA,CFILEC,CFORMAT
      CHARACTER*240 CLINE
      REAL          THK,DEPTH,FLAG,alfadt,betads,ROFF
      INTEGER       IOS,K,KDM,KI,KK,KP
C
c------------------------------------------------------------------------
      real    sigloc,dsiglocdt,dsiglocds
      real    c1l,c2l,c3l,c4l,c5l,c6l,c7l
      real    r,s,t,prs
c
c --- sub-coefficients for locally referenced sigma
c --- a fit towards Jackett & McDougall (1995)
      real, parameter, dimension(7) ::
     &  alphap = (/ -0.1364705627213484   , 0.04681812123458564,
     &               0.80700383913187     ,-0.007453530323180844,
     &              -0.002944183249153631 , 0.00003435702568990446,
     &               0.0000348657661057688 /)
     & ,betap  = (/  0.05064226654169138  ,-0.0003571087848996894,
     &              -0.0000876148051892879, 5.252431910751829e-6,
     &               1.579762259448864e-6 ,-3.466867400295792e-8,
     &              -1.687643078774232e-8 /)
     & ,gammap = (/ -5.526396144304812e-6 , 4.885838128243163e-8,
     &               9.96026931578033e-9  ,-7.251389796582352e-10,
     &              -3.987360250058777e-11, 4.006307891935698e-12,
     &               8.26367520608008e-13 /)
c
c --- HYCOM pressure to bar, for locally referenced equations
      real, parameter :: prs2pb=0.1/9806.0  !accurate value
*     real, parameter :: prs2pb=1.e-5       !original value
c
c --- locally referenced sigma, a fit towards Jackett & McDougall (1995)
c --- t: potential temperature; s: psu; prs: pressure
      c1l(prs)=alphap(1)+prs2pb*prs*(betap(1)+prs2pb*prs*gammap(1))
      c2l(prs)=alphap(2)+prs2pb*prs*(betap(2)+prs2pb*prs*gammap(2))
      c3l(prs)=alphap(3)+prs2pb*prs*(betap(3)+prs2pb*prs*gammap(3))
      c4l(prs)=alphap(4)+prs2pb*prs*(betap(4)+prs2pb*prs*gammap(4))
      c5l(prs)=alphap(5)+prs2pb*prs*(betap(5)+prs2pb*prs*gammap(5))
      c6l(prs)=alphap(6)+prs2pb*prs*(betap(6)+prs2pb*prs*gammap(6))
      c7l(prs)=alphap(7)+prs2pb*prs*(betap(7)+prs2pb*prs*gammap(7))
      sigloc(t,s,prs)=c1l(prs)+c3l(prs)*s+
     &       t*(c2l(prs)+c5l(prs)*s+t*(c4l(prs)+c7l(prs)*s+c6l(prs)*t))
      dsiglocdt(t,s,prs)=(c2l(prs)+c5l(prs)*s+
     &       2.0*t*(c4l(prs)+c7l(prs)*s+1.5*c6l(prs)*t))
      dsiglocds(t,s,prs)=(c3l(prs)+t*(c5l(prs)+t*c7l(prs)))
c
c------------------------------------------------------------------------
      REAL, ALLOCATABLE :: SI(:,:),P(:)
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
     +    'Usage:  hycom_profile_locsig archv.txt archs.txt'
        CALL EXIT(1)
      ENDIF
C
C     OPEN ALL FILES.
C
      OPEN(UNIT=11, FILE=CFILEA, FORM='FORMATTED', STATUS='OLD',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILEA)
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=21, FILE=CFILEC, FORM='FORMATTED', STATUS='NEW',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILEC)
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(5)
      ENDIF
C
C     COPY PROFILE HEADER TO OUTPUT.
C
      DO K= 1,99
        READ( 11,'(a)')      CLINE
        IF     (CLINE(1:5).EQ.'#  k ') then
          EXIT
        ENDIF
        WRITE(21,'(a)') TRIM(CLINE)
      ENDDO
C
C     READ THE ISOPYCNAL PROFILE, TO GET KDM.
C
      DO K= 1,99999
        READ(11,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
      ENDDO
      KDM = K-1
C
C     RE-READ THE ISOPYCNAL PROFILE.
C
      ALLOCATE( P(KDM+1), SI(KDM,7) )
C
      REWIND(11)
      DO K= 1,99
        READ( 11,'(a)') CLINE
        IF     (CLINE(1:5).EQ.'#  k ') then
          EXIT
        ENDIF
      ENDDO
      P(1) =  0.0
      DO K= 1,KDM
        READ(11,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'Error: inconsistent input profile'
          CALL EXIT(6)
        ENDIF
        READ(CLINE,*) KI,(SI(K,KK),KK=1,5),THK,DEPTH
        P(K+1) = P(K) + THK
        IF     (THK.EQ.0.0) THEN
          DO KK= 1,5
            SI(K,KK)=SI(K-1,KK)
          ENDDO !kk
        ENDIF
      ENDDO
      CLOSE(11)
C
C     CONVERT TO SIGLOC, KEEP ORIGINAL AND ALT AS TRACERS
C
      K=1
        SI(K,6) = SI(K,5)
        SI(K,5) = 0.0
        SI(K,7) = 0.0
      DO K= 2,KDM
        SI(K,6) = SI(K,5)
        alfadt = 0.5*(dsiglocdt(SI(K-1,3),SI(K-1,4),9806.0*P(K))+
     &                dsiglocdt(SI(K  ,3),SI(K  ,4),9806.0*P(K)) )*
     &               (SI(K-1,3)-SI(K,3))
        betads = 0.5*(dsiglocds(SI(K-1,3),SI(K-1,4),9806.0*P(K))+
     &                dsiglocds(SI(K  ,3),SI(K  ,4),9806.0*P(K)) )*
     &               (SI(K-1,4)-SI(K,4))
        SI(K,5) = SI(K-1,5)-alfadt-betads
        alfadt = dsiglocdt(0.5*(SI(K-1,3)+SI(K,3)),
     &                     0.5*(SI(K-1,4)+SI(K,4)),9806.0*P(K))*
     &                         (SI(K-1,3)-SI(K,3))
        betads = dsiglocds(0.5*(SI(K-1,3)+SI(K,3)),
     &                     0.5*(SI(K-1,4)+SI(K,4)),9806.0*P(K))*
     &                         (SI(K-1,4)-SI(K,4))
        SI(K,7) = SI(K-1,7)-alfadt-betads
      ENDDO
C     make lowest layer potential density the same as before
      ROFF =  SI(KDM,6) - SI(KDM,5)
      DO K= 1,KDM
        SI(K,5) = SI(K,5)+ROFF
      ENDDO
      ROFF =  SI(KDM,6) - SI(KDM,7)
      DO K= 1,KDM
        SI(K,7) = SI(K,7)+ROFF
      ENDDO
C
C     OUTPUT
C
        WRITE(CFORMAT,'(a)')
     &    '(3a)'
        WRITE(21,CFORMAT)
     &      '#  k',
     &      '    utot    vtot  p.temp    saln  p.dens',
     &      '    thkns      dpth  tracer  tracer'
C
          WRITE(CFORMAT,'(a)')
     &      '(i4,2f8.2,3f8.4,f9.3,f10.3,2f8.4)'
C
        DO K= 1,KDM
          THK = P(K+1) - P(K)
          WRITE(21,CFORMAT)
     &      K,(SI(K,KK),KK=1,5),THK,0.5*(P(K)+P(K+1)),
     &        (SI(K,KK),KK=6,7)
        ENDDO !k
      CLOSE(21)
      END

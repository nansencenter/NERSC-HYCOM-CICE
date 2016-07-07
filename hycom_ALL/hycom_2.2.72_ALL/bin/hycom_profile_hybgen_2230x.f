      PROGRAM HYCOM_PROFILE_HYBGEN
      IMPLICIT NONE
C
C  hycom_profile_hybgen - Usage: hycom_profile_hybgen archv.txt blkdat.input archg.txt
C
C                 apply the HYCOM hybrid grid generator to a HYCOM profile.
C                 note that the velocities are unchanged on exit.
C                 based on hybgen version 2.2.30
C
C   archv.txt    is assumed to be an HYCOM archive text profile file
C   blkdat.input is assumed to be an HYCOM parameter file (version 2.2.30)
C   archg.txt    will be the output text profile file
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  August 2003
C
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      CHARACTER*240 CFILEA,CFILEB,CFILEC
      CHARACTER*240 CLINE
      REAL          THK,THKIN,FLAG,ONEM
      REAL          U(99),V(99),T(99),S(99),R(99),DP(99),P(99+1)
      REAL          RIN(99),TIN(99),SIN(99),PIN(99+1)
      REAL*8        SUMI(3),SUMO(3),SUMD(3)
      INTEGER       IOS,K,KDM
      LOGICAL       LDEBUG
C
      LOGICAL       ISOPCM
      INTEGER       KDMBLK,NHYBRD,NSIGMA,HYBMAP,HYBFLG
      REAL          HYBISO,HYBTHK,HYBRLX,QHYBRLX
      REAL          DP00,DP00X,DP00F,DS00,DS00X,DS00F,DP00I,THETA(999)
      REAL          DP0K(999),DS0K(999),DSSK(999),DEPTHS,ISOTOP
      REAL          THKBOT
c
      integer     mnproc,lp
      common/xxx/ mnproc,lp
c
      integer   thflag
      common/th/thflag
c
      mnproc = 1
      lp     = 6
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      CALL GETARG(0,CLINE)
      K = LEN_TRIM(CLINE)
      LDEBUG = CLINE(K-5:K) .EQ. "_debug"
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEB)
        CALL GETARG(3,CFILEC)
      ELSE
        WRITE(6,*)
     +    'Usage: hycom_profile_hybgen archv.txt blkdat.input archg.txt'
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
      OPEN(UNIT=99, FILE=CFILEB, FORM='FORMATTED', STATUS='OLD',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',CFILEB(1:LEN_TRIM(CFILEB))
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
C     INPUT BLKDAT PARAMETERS
C
      DO K= 1,10
        READ(99,*)
      ENDDO
c
c --- 'kdm   ' = number of layers
c --- 'nhybrd' = number of hybrid levels (0=all isopycnal)
c --- 'nsigma' = number of sigma  levels (nhybrd-nsigma z-levels)
c --- 'dp00'   = deep    z-level spacing minimum thickness (m)
c --- 'dp00x'  = deep    z-level spacing maximum thickness (m)
c --- 'dp00f'  = deep    z-level spacing stretching factor (1.0=const.z)
c --- 'ds00'   = shallow z-level spacing minimum thickness (m)
c --- 'ds00x'  = shallow z-level spacing maximum thickness (m)
c --- 'ds00f'  = shallow z-level spacing stretching factor (1.0=const.z)
c --- 'dp00i'  = deep iso-pycnal spacing minimum thickness (m)
c --- 'isotop' = shallowest depth for isopycnal layers (m), <0 from file
c
c --- or, in place of 'dp00','dp00x','dp00f','ds00','ds00x','ds00f' specify:
c --- 'dp0k  ' = layer k deep    z-level spacing minimum thickness (m)
c ---              k=1,kdm; dp0k must be zero for k>nhybrd
c --- 'ds0k  ' = layer k shallow z-level spacing minimum thickness (m)
c ---              k=1,nsigma
c
      call blkini(kdmblk,'kdm   ')
      call blkini(nhybrd,'nhybrd')
      call blkini(nsigma,'nsigma')
      call blkinr2(dp00,k,
     &                   'dp00  ','(a6," =",f10.4," m")',
     &                   'dp0k  ','(a6," =",f10.4," m")' )
      call flush(lp)
      if     (k.eq.1) then !dp00
        call blkinr(dp00x, 'dp00x ','(a6," =",f10.4," m")')
        call blkinr(dp00f, 'dp00f ','(a6," =",f10.4," ")') 
        call blkinr(ds00,  'ds00  ','(a6," =",f10.4," m")')
        call blkinr(ds00x, 'ds00x ','(a6," =",f10.4," m")')
        call blkinr(ds00f, 'ds00f ','(a6," =",f10.4," ")') 
      else !dp0k
        dp0k(1) = dp00
        dp00    = -1.0  !signal that dp00 is not input
        do k=2,kdmblk
          call blkinr(dp0k(k), 'dp0k  ','(a6," =",f10.4," m")')
c
          if      (k.gt.nhybrd .and. dp0k(k).ne.0.0) then
            write(lp,'(/ a,i3 /)')
     &        'error - dp0k must be zero for k>nhybrd'
            call flush(lp)
            stop '(blkdat)'
          endif !k>nhybrd&dp0k(k)!=0
        enddo !k
        do k=1,nsigma
          call blkinr(ds0k(k), 'ds0k  ','(a6," =",f10.4," m")')
        enddo !k
      endif !dp00:dp0k
      call blkinr(dp00i, 'dp00i ','(a6," =",f10.4," m")')
      call blkinr(isotop,'isotop','(a6," =",f10.4," m")')
      call flush(lp)
c
c --- 'thflag' = reference pressure flag (0=Sigma-0, 2=Sigma-2)
c
      READ(99,*) !'saln0 '
      READ(99,*) !'locsig'
      READ(99,*) !'kapflg'
      call blkini(thflag,'thflag')
      READ(99,*) !'thbase'
      READ(99,*) !'vsigma'
c
c --- 'sigma ' = isopycnal layer target densities (sigma units)
c
      do k=1,kdmblk
        call blkinr(THETA(k),'sigma ','(a6," =",f10.4," sigma-0")')
      enddo
      call flush(lp)
c
c --- skip 'iniflg' to 'batrop'.
c
      DO K= 1,21
        READ(99,*)
      ENDDO
c
c --- 'hybrlx' = HYBGEN: inverse relaxation coefficient (time steps)
c                 (1.0 for no relaxation)
c --- 'hybthk' = HYBGEN: rato of layer vs its neighbors should be > hybthk
c                 (0.0 for old behaviour; only active when < dp0k )
c --- 'hybiso' = HYBGEN: Use PCM if layer is within hybiso of target density
c                 (0.0 for no PCM; large to recover pre-2.2.09 behaviour)
c --- 'hybmap' = HYBGEN:  remapper  flag (0=PCM, 1=PLM,    2=PPM,  3=WENO-like)
c --- 'hybflg' = HYBGEN:  generator flag (0=T&S, 1=th&S,   2=th&T)
c
      call blkinr(hybrlx,'hybrlx','(a6," =",f10.4," time steps")')
      call blkinr(hybthk,'hybthk','(a6," =",f10.4," ")')
      call blkinr(hybiso,'hybiso','(a6," =",f10.4," kg/m^3")')
      call blkini(hybmap,'hybmap')
      call blkini(hybflg,'hybflg')
      call flush(lp)
c
      qhybrlx = 1.0/hybrlx
      isopcm  = hybiso.gt.0.0  !use PCM for isopycnal layers?
C
      CLOSE(99)
C
C     INITIALIZE FIXED GRID.
C
      CALL GEOPAR(DP00,DP00X,DP00F,DS00,DS00X,DS00F,DP00I,
     +            NHYBRD,NSIGMA,KDMBLK,
     +            DP0K,DS0K,DSSK)
*     call flush(lp)
C
C     COPY PROFILE HEADER TO OUTPUT.
C
      DO K= 1,99
        READ( 11,'(a)') CLINE
        IF     (CLINE(1:5).EQ.'#  k ') then
          EXIT
        ENDIF
        WRITE(21,'(a)') TRIM(CLINE)
      ENDDO
C
C     READ THE PROFILE.
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
C     TRANSFORM THE PROFILE.
C
      ONEM = 9806.0 !g/thref
      DO K= 1,KDM
        DP(K) = (P(K+1) - P(K))*ONEM
C
        RIN(K) = R(K)
        TIN(K) = T(K)
        SIN(K) = S(K)
        PIN(K) = P(K)
      ENDDO
      PIN(KDM+1) = P(KDM+1)
C
      DEPTHS = P(KDM+1)*ONEM
      ISOTOP = ISOTOP*ONEM
      THKBOT = 10.0
C
      CALL HYBGEN(T,S,R,DP,THETA,KDM,
     +            NHYBRD,ISOPCM,HYBMAP,HYBFLG,HYBISO,HYBTHK, QHYBRLX,
     +            DP00I,DP0K,DS0K,DSSK,DEPTHS,ISOTOP,THKBOT,
     +            LDEBUG)
C
      DO K= 1,KDM
        P(K+1) = P(K) + DP(K)/ONEM
      ENDDO
C
C     DIAGNOSTIC OUTPUT
C
      WRITE(6,*)
      WRITE(6,'(a,a)')
     &  '#  k ',
     &  '   sigma    dens    dens    deld    thkns    thkns'
      DO K= 1,KDM
        THK   = P(  K+1) - P(  K)
        THKIN = PIN(K+1) - PIN(K)
        IF     (ABS(THK-THKIN).LT.0.001) THEN
          WRITE(6,'(i4,1x,4f8.3,2f9.3)')
     &      K,THETA(K),RIN(K),R(K),R(K)-THETA(K),THKIN,THK
        ELSE
          WRITE(6,'(i4,1x,4f8.3,2f9.3,a)')
     &      K,THETA(K),RIN(K),R(K),R(K)-THETA(K),THKIN,THK,' *'
        ENDIF
      ENDDO
C
      SUMI(:)=0.0
      SUMO(:)=0.0
      DO K= 1,KDM
        THK   = P(  K+1) - P(  K)
        THKIN = PIN(K+1) - PIN(K)
        SUMI(1) = SUMI(1) + THKIN*TIN(K)
        SUMI(2) = SUMI(2) + THKIN*SIN(K)
        SUMI(3) = SUMI(3) + THKIN*RIN(K)
        SUMO(1) = SUMO(1) + THK  *T(  K)
        SUMO(2) = SUMO(2) + THK  *S(  K)
        SUMO(3) = SUMO(3) + THK  *R(  K)
      ENDDO
      DO K=1,3
        SUMI(K) = SUMI(K)/PIN(KDM+1)
        SUMO(K) = SUMO(K)/P(  KDM+1)
        SUMD(K) = SUMO(K) - SUMI(K)
      ENDDO
      WRITE(*,'(/A)')       '        LABEL:  T  S  R'
      WRITE(*,'(A,3e14.6)') ' INPUT TOTALS: ',SUMI(:)
      WRITE(*,'(A,3e14.6)') 'OUTPUT TOTALS: ',SUMO(:)
      WRITE(*,'(A,3e14.6)') 'OUT-IN TOTALS: ',SUMD(:)
      WRITE(*,*)
C
C     OUTPUT THE TRANSFORMED PROFILE.
C
      WRITE(21,'(3a)')
     &      '#  k',
     &      '    utot    vtot  p.temp    saln  p.dens',
     &      '    thkns      dpth'
      DO K= 1,KDM
        THK = P(K+1) - P(K)
        WRITE(21,'(i4,2f8.2,3f8.3,f9.3,f10.3)')
     &    K,U(K),V(K),T(K),S(K),R(K),THK,0.5*(P(K)+P(K+1))
      ENDDO
      END

      subroutine blkinr(rvar,cvar,cfmt)
      implicit none
c
      real      rvar
      character cvar*6,cfmt*(*)
c
c     read in one real value
c
      character*6 cvarin
c
      read(99,*) rvar,cvarin
      write(6,cfmt) cvarin,rvar
c
      if     (cvar.ne.cvarin) then
        write(6,*) 
        write(6,*) 'error in blkinr - input ',cvarin,
     +                      ' but should be ',cvar
        write(6,*) 
        stop
      endif
      return
      end
      subroutine blkinr2(rvar,nvar,cvar1,cfmt1,cvar2,cfmt2)
      implicit none
c
      real      rvar
      integer   nvar
      character cvar1*6,cvar2*6,cfmt1*(*),cfmt2*(*)
c
c     read in one real value from stdin,
c     identified as either cvar1 (return nvar=1) or cvar2 (return nvar=2)
c
      character*6 cvarin
c
      read(99,*) rvar,cvarin
c
      if     (cvar1.eq.cvarin) then
        nvar = 1
        write(6,cfmt1) cvarin,rvar
      elseif (cvar2.eq.cvarin) then
        nvar = 2
        write(6,cfmt2) cvarin,rvar
      else
        write(6,*)
        write(6,*) 'error in blkinr2 - input ',cvarin,
     +                      ' but should be ',cvar1,' or ',cvar2
        write(6,*)
        stop
      endif
      return
      end
      subroutine blkini(ivar,cvar)
      implicit none
c
      integer     ivar
      character*6 cvar
c
c     read in one integer value
c
      character*6 cvarin
c
      read(99,*) ivar,cvarin
      write(6,6000) cvarin,ivar
c
      if     (cvar.ne.cvarin) then
        write(6,*) 
        write(6,*) 'error in blkini - input ',cvarin,
     +                      ' but should be ',cvar
        write(6,*) 
        stop
      endif
      return
 6000 format(a6,' =',i6)
      end

      subroutine geopar(dp00,dp00x,dp00f,ds00,ds00x,ds00f,dp00i,
     &                  nhybrd,nsigma,kdm,
     &                  dp0k,ds0k,dssk)
      implicit none
c
      integer nhybrd,nsigma,kdm
      real    dp00,dp00x,dp00f,ds00,ds00x,ds00f,dp00i
      real    dp0k(kdm),ds0k(kdm),dssk(kdm)
c
      integer     mnproc,lp
      common/xxx/ mnproc,lp
c
c --- set up model parameters related to geography
c
      real      dp0kf,dpm,dpms,ds0kf,dsm,dsms,onem,qonem
      integer   k
c
      onem  = 9806.0 !g/thref
      qonem = 1.0/onem
c
c --- calculate dp0k and ds0k?
      if     (dp00.lt.0.0) then
c ---   dp0k and ds0k already input
        dp00 =onem*dp0k(1)
        dp00x=onem*dp0k(kdm-1)
        dp00i=onem*dp00i
        dpms = 0.0
        do k=1,kdm
          dpm     = dp0k(k)
          dpms    = dpms + dpm
          dp0k(k) = dp0k(k)*onem
          if     (mnproc.eq.1) then
          write(lp,135) k,dp0k(k)*qonem,dpm,dpms
          endif
          if     (mnproc.eq.-99) then  ! bugfix that prevents optimization
            write(6,*) 'geopar: dp0k   = ',dp0k(k),k,mnproc
          endif
        enddo !k
        dsms = 0.0
        do k=1,nsigma
          dsm     = ds0k(k)
          dsms    = dsms + dsm
          ds0k(k) = ds0k(k)*onem
          if     (mnproc.eq.1) then
          write(lp,130) k,ds0k(k)*qonem,dsm,dsms
          endif
          if     (mnproc.eq.-99) then  ! bugfix that prevents optimization
            write(6,*) 'geopar: ds0k   = ',ds0k(k),k,mnproc
          endif
        enddo !k
        if     (mnproc.eq.1) then
        write(lp,*)
        endif
      else
c ---   calculate dp0k and ds0k
c
c ---   logorithmic k-dependence of dp0 (deep z's)
        dp00 =onem*dp00
        dp00x=onem*dp00x
        dp00i=onem*dp00i
        if     (nhybrd.eq.0) then
*         dp0k(1)=thkmin*onem
          dp0k(1)=  20.0*onem
        else
          dp0k(1)=dp00
        endif
        dpm  = dp0k(1)*qonem
        dpms = dpm
        write(lp,*)
        write(lp,135) 1,dp0k(1)*qonem,dpm,dpms
 135    format('dp0k(',i2,') =',f7.2,' m',
     &            '    thkns =',f7.2,' m',
     &            '    depth =',f8.2,' m')
c
        dp0kf=1.0
        do k=2,kdm
          dp0kf=dp0kf*dp00f
          if     (k.le.nhybrd) then
            dp0k(k)=min(dp00*dp0kf,dp00x)
          else
            dp0k(k)=0.0
          endif
          dpm  = dp0k(k)*qonem
          dpms = dpms + dpm
          write(lp,135) k,dp0k(k)*qonem,dpm,dpms
          if     (mnproc.eq.-99) then  ! bugfix that prevents optimization
            write(6,*) 'geopar: dp0kf  = ',dp0kf,    mnproc
            write(6,*) 'geopar: dp0k   = ',dp0k(k),k,mnproc
          endif
        enddo
c
c ---   logorithmic k-dependence of ds0 (shallow z-s)
        ds00 =onem*ds00
        ds00x=onem*ds00x
        if     (nhybrd.eq.0) then
*         ds0k(1)=thkmin*onem
          ds0k(1)=  20.0*onem
        else
          ds0k(1)=ds00
        endif
        dsm  = ds0k(1)*qonem
        dsms = dsm
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,130) 1,ds0k(1)*qonem,dsm,dsms
        endif
 130    format('ds0k(',i2,') =',f7.2,' m',
     &            '    thkns =',f7.2,' m',
     &            '    depth =',f8.2,' m')
c
        ds0kf=1.0
        do k=2,nsigma
          ds0kf=ds0kf*ds00f
          ds0k(k)=min(ds00*ds0kf,ds00x)
          dsm  = ds0k(k)*qonem
          dsms = dsms + dsm
          if     (mnproc.eq.1) then
          write(lp,130) k,ds0k(k)*qonem,dsm,dsms
          endif
          if     (mnproc.eq.-99) then  ! bugfix that prevents optimization
            write(6,*) 'geopar: ds0kf  = ',ds0kf,    mnproc
            write(6,*) 'geopar: ds0k   = ',ds0k(k),k,mnproc
          endif
        enddo
        if     (mnproc.eq.1) then
        write(lp,*)
        endif
      endif !input:calculate dp0k,ds0k
c
c --- sigma-depth scale factors
      do k=1,nsigma
        dssk(k)=ds0k(k)/dsms  ! onem * fraction of depths in sigma layer k
      enddo
      do k= nsigma+1,kdm
        ds0k(k)=dp0k(k)
        dssk(k)=0.0           ! these layers are zero in sigma mode
      enddo
c
      return
      end

      real function sig(t,s)
      real t,s
c
      integer   thflag
      common/th/thflag
c
c --- equation of state
c
      real sig0,sig2
c
      if     (thflag.eq.0) then
        sig = sig0(t,s)
      else
        sig = sig2(t,s)
      endif
      end

      real function dsigdt(t,s)
      real t,s
c
      integer   thflag
      common/th/thflag
c
c --- d(sig)/dt
c
      real dsigdt0,dsigdt2
c
      if     (thflag.eq.0) then
        dsigdt = dsigdt0(t,s)
      else
        dsigdt = dsigdt2(t,s)
      endif
      end

      real function dsigds(t,s)
      real t,s
c
      integer   thflag
      common/th/thflag
c
c --- d(sig)/ds
c
      real dsigds0,dsigds2
c
      if     (thflag.eq.0) then
        dsigds = dsigds0(t,s)
      else
        dsigds = dsigds0(t,s)
      endif
      end
c
      real function tofsig(r,s)
      real r,s
c
      integer   thflag
      common/th/thflag
c
c --- temp (deg c) as a function of sigma and salinity (mil)
c
      real tofsig0,tofsig2
c
      if     (thflag.eq.0) then
        tofsig = tofsig0(r,s)
      else
        tofsig = tofsig2(r,s)
      endif
      end

      real function sofsig(r,t)
      real r,t
c
      integer   thflag
      common/th/thflag
c
c --- salinity (mil) as a function of sigma and temperature (deg c)
c
      real sofsig0,sofsig2
c
      if     (thflag.eq.0) then
        sofsig = sofsig0(r,t)
      else
        sofsig = sofsig2(r,t)
      endif
      end

      real function sig0(tx,sx)
      real tx,sx
c
c --- equation of state
c
      include 'hycom_profile_hybgen_SIGMA0.h'
c
      sig0 = sig(tx,sx)
      end

      real function dsigdt0(tx,sx)
      real tx,sx
c
c --- d(sig)/dt
c
      include 'hycom_profile_hybgen_SIGMA0.h'
c
      dsigdt0 = dsigdt(tx,sx)
      end

      real function dsigds0(tx,sx)
      real tx,sx
c
c --- d(sig)/ds
c
      include 'hycom_profile_hybgen_SIGMA0.h'
c
      dsigds0 = dsigds(tx,sx)
      end
c
      real function tofsig0(rx,sx)
      real rx,sx
c
c --- temp (deg c) as a function of sigma and salinity (mil)
c
      include 'hycom_profile_hybgen_SIGMA0.h'
c
      tofsig0 = tofsig(rx,sx)
      end

      real function sofsig0(rx,tx)
      real rx,tx
c
c --- salinity (mil) as a function of sigma and temperature (deg c)
c
      include 'hycom_profile_hybgen_SIGMA0.h'
c
      sofsig0 = sofsig(rx,tx)
      end

      real function sig2(tx,sx)
      real tx,sx
c
c --- equation of state
c
      include 'hycom_profile_hybgen_SIGMA2.h'
c
      sig2 = sig(tx,sx)
      end

      real function dsigdt2(tx,sx)
      real tx,sx
c
c --- d(sig)/dt
c
      include 'hycom_profile_hybgen_SIGMA2.h'
c
      dsigdt2 = dsigdt(tx,sx)
      end

      real function dsigds2(tx,sx)
      real tx,sx
c
c --- d(sig)/ds
c
      include 'hycom_profile_hybgen_SIGMA2.h'
c
      dsigds2 = dsigds(tx,sx)
      end
c
      real function tofsig2(rx,sx)
      real rx,sx
c
c --- temp (deg c) as a function of sigma and salinity (mil)
c
      include 'hycom_profile_hybgen_SIGMA2.h'
c
      tofsig2 = tofsig(rx,sx)
      end

      real function sofsig2(rx,tx)
      real rx,tx
c
c --- salinity (mil) as a function of sigma and temperature (deg c)
c
      include 'hycom_profile_hybgen_SIGMA2.h'
c
      sofsig2 = sofsig(rx,tx)
      end

      subroutine hybgen(temp,saln,th3d,dp,theta,kdm,
     &                  nhybrd,isopcm,hybmap,hybflg,
     &                  hybiso,hybthk, qhybrlx,
     &                  dp00i,dp0k,ds0k,dssk,depths,topiso,
     &                  thkbot, ldebug)
      implicit none
c
      logical isopcm,ldebug
      integer kdm, nhybrd,hybmap,hybflg
      real     temp(1,1,kdm,1),
     &         saln(1,1,kdm,1),
     &         th3d(1,1,kdm,1),
     &           dp(1,1,kdm,1),
     &        theta(1,1,kdm)
      real    hybiso,hybthk,qhybrlx,thkbot
      real    dp00i,dp0k(kdm),ds0k(kdm),dssk(kdm),
     &        depths(1,1),topiso(1,1)
c
      integer     mnproc,lp
      common/xxx/ mnproc,lp
c
      integer klist(1,1)
      real    p(1,1,kdm+1),dpmixl(1,1,1),q2(1,1,1,1),q2l(1,1,1,1)
      real    tracer(1,1,1,1,1),trcflg(1)
      logical mxlkta,thermo
      integer j,kk,n, nstep, i0,j0,itest,jtest
      real    epsil,onem,tenm,tencm,onemm,qonem,thbase
c
      real     sig,dsigdt,dsigds,tofsig,sofsig
      external sig,dsigdt,dsigds,tofsig,sofsig
c
      integer, parameter :: mxtrcr=1
      integer, parameter :: ntracr=0
      logical, parameter :: mxlmy =.false.
c
c --- ---------------------
c --- hybrid grid generator
c --- ---------------------
c
      logical, parameter :: lunmix=.true.     !unmix a too light deepest layer
      logical, parameter :: lconserve=.false. !explicitly conserve each column
c
      double precision asum(  mxtrcr+4,3)
      real             offset(mxtrcr+4)
c
      logical lcm(kdm)             !use PCM for some layers?
      real    s1d(kdm,mxtrcr+4),   !original scalar fields
     &        f1d(kdm,mxtrcr+4),   !final    scalar fields
     &        c1d(kdm,mxtrcr+4,3), !interpolation coefficients
     &        dpi( kdm),           !original layer thicknesses, >= dpthin
     &        dprs(kdm),           !original layer thicknesses
     &        pres(kdm+1),         !original layer interfaces
     &        prsf(kdm+1),         !final    layer interfaces
     &        qhrlx( kdm+1),       !relaxation coefficient, from qhybrlx
     &        dp0ij( kdm),         !minimum layer thickness
     &        dp0cum(kdm+1)        !minimum interface depth
      real    p_hat,p_hat0,p_hat2,p_hat3,hybrlx,
     &        delt,deltm,dels,delsm,q,qtr,qts,thkbop,
     &        zthk,dpthin,dpbot,dpmid,dptop,dpinc
      integer i,k,ka,kp,ktr,l,fixlay,nums1d
      character*12 cinfo
c
      double precision, parameter :: dsmll=1.0d-8
      double precision, parameter ::   zp5=0.5    !for sign function
c
c --- c u s h i o n   function (from Bleck & Benjamin, 1992):
c --- if delp >= qqmx*dp0 >>  dp0, -cushn- returns -delp-
c --- if delp <= qqmn*dp0 << -dp0, -cushn- returns  -dp0-
c
      real       qqmn,qqmx,cusha,cushb
      parameter (qqmn=-4.0, qqmx=2.0)  ! shifted range
*     parameter (qqmn=-2.0, qqmx=4.0)  ! traditional range
*     parameter (qqmn=-4.0, qqmx=6.0)  ! somewhat wider range
      parameter (cusha=qqmn**2*(qqmx-1.0)/(qqmx-qqmn)**2)
      parameter (cushb=1.0/qqmn)
c
      real qq,cushn,delp,dp0
*     include 'stmt_fns.h'
      qq(   delp,dp0)=max(qqmn, min(qqmx, delp/dp0))
      cushn(delp,dp0)=dp0*
     &                (1.0+cusha*(1.0-cushb*qq(delp,dp0))**2)*
     &                max(1.0, delp/(dp0*qqmx))
c
      mxlkta = .false.
      thermo = .true.
      i      = 1
      j      = 1
      kk     = kdm
      n      = 1
      nstep  = 1
      thbase = 0.0
      i0     = 0
      j0     = 0
      if     (ldebug) then
        itest = 1
        jtest = 1
      else
        itest = 0
        jtest = 0
      endif
c
      epsil = 1.0e-11
      onem  = 9806.0 !g/thref
      qonem = 1.0/onem
      tenm  = onem*10.0
      tencm = onem/10.0
      onemm = onem/1000.0
c
      dpthin = 0.001*onemm
      thkbop = thkbot*onem
      hybrlx = 1.0/qhybrlx
c
      if (mxlmy) then
        nums1d = ntracr + 4
      else
        nums1d = ntracr + 2
      endif
c
      if     (.not.isopcm) then
        do k=1,nhybrd
          lcm(k) = .false.  !use same remapper for all layers
        enddo !k
        do k=nhybrd+1,kk
          lcm(k) = .true.   !purely isopycnal layers use PCM
        enddo !k
      endif
c
*     do l=1,isp(j)
*     do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
c
      dp0cum(1)=0.0
      qhrlx( 1)=1.0 !no relaxation in top layer
      dp0ij( 1)=min( dp0k(1), max( ds0k(1), dssk(1)*depths(i,j) ) )
      dp0cum(2)=dp0cum(1)+dp0ij(1)
      qhrlx( 2)=1.0 !no relaxation in top layer
      p(i,j, 2)=p(i,j,1)+dp(i,j,1,n)
      do k=2,kk
c ---   q is dp0k(k) when in surface fixed coordinates
c ---   q is dp00i   when much deeper than surface fixed coordinates
        if     (dp0k(k).le.dp00i) then
          q  =      dp0k(k)
          qts=      0.0     !0 at dp0k
        else
          q  = max( dp00i,
     &              dp0k(k) * dp0k(k)/
     &                        max( dp0k( k),
     &                             p(i,j,k)-dp0cum(k) ) )
          qts= 1.0 - (q-dp00i)/(dp0k(k)-dp00i)  !0 at dp0k, 1 at dp00i
        endif
        qhrlx( k+1)=1.0/(1.0 + qts*(hybrlx-1.0))  !1 at  dp0k, qhybrlx at dp00i
        dp0ij( k)  =min( q,max( ds0k(k), dssk(k)*depths(i,j) ) )
        dp0cum(k+1)=dp0cum(k)+dp0ij(k)
        p(i,j, k+1)=p(i,j,k)+dp(i,j,k,n)
c
        if (i.eq.itest .and. j.eq.jtest) then
          write (lp,'(i6,1x,3f8.3,2f9.3)')
     .    k,dp0ij(k)*qonem,q*qonem,p(i,j,k)*qonem-dp0cum(k)*qonem,
     .                             p(i,j,k)*qonem,dp0cum(k)*qonem
        endif !debug
      enddo !k
c
c --- identify the always-fixed coordinate layers
      fixlay = 1  !layer 1 always fixed
      do k= 2,nhybrd
        if     (dp0cum(k).ge.topiso(i,j)) then
          exit  !layers k to nhybrd can be isopycnal
        endif
        qhrlx(k+1)=1.0  !no relaxation in fixed layers
        fixlay = fixlay+1
      enddo !k
      if (i.eq.itest .and. j.eq.jtest) then
        write(lp,'(a,i3)')
     &        'hybgen, always-fixed coordinate layers: 1 to ',
     &        fixlay
        call flush(lp)
      endif !debug
c
      if (i.eq.itest .and. j.eq.jtest) then
        write (lp,'(a/(i6,1x,2f8.3,2f9.3,f9.3))')
     .  'hybgen:   thkns  minthk     dpth  mindpth   hybrlx',
     .  (k,dp(i,j,k,n)*qonem,   dp0ij(k)*qonem,
     .      p(i,j,k+1)*qonem,dp0cum(k+1)*qonem,
     .      1.0/qhrlx(k+1),
     .   k=1,kk)
      endif !debug
c
c --- identify the deepest layer kp with significant thickness (> dpthin)
c
      kp = 2  !minimum allowed value
      do k=kk,3,-1
        if (p(i,j,k+1)-p(i,j,k).ge.dpthin) then
          kp=k
          exit
        endif
      enddo
c
      k=kp  !at least 2
c
      if (i.eq.itest .and. j.eq.jtest) then
        write(lp,'(a,i3)')
     &        'hybgen, deepest inflated layer:',k
        call flush(lp)
      endif !debug
c
      ka = max(k-2,1)  !k might be 2
      if     (k.gt.fixlay+1 .and.
     &        theta(i,j,k)-epsil.gt.th3d(i,j,k,n) .and.
     &         th3d(i,j,k-1,n)  .gt.th3d(i,j,k,n) .and.
     &         th3d(i,j,ka, n)  .gt.th3d(i,j,k,n)      ) then
c
c ---   water in the deepest inflated layer with significant thickness
c ---   (kp) is too light, and it is lighter than the two layers above.
c ---
c ---   this should only occur when relaxing or nudging layer thickness
c ---   and is a bug (bad interaction with tsadvc) even in those cases
c ---
c ---   entrain the entire layer into the one above
c---    note the double negative in T=T-q*(T-T'), equiv. to T=T+q*(T'-T)
        q = (p(i,j,k+1)-p(i,j,k))/(p(i,j,k+1)-p(i,j,k-1))
        if     (hybflg.eq.0) then  !T&S
          temp(i,j,k-1,n)=temp(i,j,k-1,n)-q*(temp(i,j,k-1,n) -
     &                                       temp(i,j,k,  n)  )
          saln(i,j,k-1,n)=saln(i,j,k-1,n)-q*(saln(i,j,k-1,n) -
     &                                       saln(i,j,k,  n)  )
          th3d(i,j,k-1,n)=sig(temp(i,j,k-1,n),saln(i,j,k-1,n))-thbase
        elseif (hybflg.eq.1) then  !th&S
          th3d(i,j,k-1,n)=th3d(i,j,k-1,n)-q*(th3d(i,j,k-1,n) -
     &                                       th3d(i,j,k,  n)  )
          saln(i,j,k-1,n)=saln(i,j,k-1,n)-q*(saln(i,j,k-1,n) -
     &                                       saln(i,j,k,  n)  )
          temp(i,j,k-1,n)=tofsig(th3d(i,j,k-1,n)+thbase,saln(i,j,k-1,n))
        elseif (hybflg.eq.2) then  !th&T
          th3d(i,j,k-1,n)=th3d(i,j,k-1,n)-q*(th3d(i,j,k-1,n) -
     &                                       th3d(i,j,k,  n)  )
          temp(i,j,k-1,n)=temp(i,j,k-1,n)-q*(temp(i,j,k-1,n) -
     &                                       temp(i,j,k,  n)  )
          saln(i,j,k-1,n)=sofsig(th3d(i,j,k-1,n)+thbase,temp(i,j,k-1,n))
        endif
        do ktr= 1,ntracr
          tracer(i,j,k-1,n,ktr)=tracer(i,j,k-1,n,ktr)-
     &                               q*(tracer(i,j,k-1,n,ktr) -
     &                                  tracer(i,j,k,  n,ktr)  )
        enddo !ktr
        if (mxlmy) then
          q2( i,j,k-1,n)=q2( i,j,k-1,n)-
     &                     q*(q2( i,j,k-1,n)-q2( i,j,k,n))
          q2l(i,j,k-1,n)=q2l(i,j,k-1,n)-
     &                     q*(q2l(i,j,k-1,n)-q2l(i,j,k,n))
        endif
c ---   entrained the entire layer into the one above, so now kp=kp-1
        p(i,j,k) = p(i,j,k+1)
        kp = k-1
        if (i.eq.itest .and. j.eq.jtest) then
          write(lp,'(a,i3,f6.3,5f8.3)')
     &      'hybgen, 11(+):',
     &      k-1,q,temp(i,j,k-1,n),saln(i,j,k-1,n),
     &          th3d(i,j,k-1,n)+thbase,theta(i,j,k-1)+thbase
          call flush(lp)
        endif !debug
      endif
c
      if     (lunmix        .and. !usually .true.
     &        k.gt.fixlay+1 .and.
     &        theta(i,j,k)-epsil.gt.th3d(i,j,k,n) .and.
     &        theta(i,j,k-1)    .lt.th3d(i,j,k,n) .and.
     &        ( th3d(i,j,k,n)- th3d(i,j,k-1,n)).gt.
     &        (theta(i,j,k)  -theta(i,j,k-1)  )*0.001  ) then
c
c ---   water in the deepest inflated layer with significant thickness
c ---   (kp) is too light
c ---
c ---   split layer into 2 sublayers, one near the desired density
c ---   and one exactly matching the T&S properties of layer k-1.
c ---   To prevent "runaway" T or S, the result satisfies either
c ---     abs(T.k - T.k-1) <= abs(T.k-2 - T.k-1) or
c ---     abs(S.k - S.k-1) <= abs(S.k-2 - S.k-1) 
c ---   It is also limited to a 50% change in layer thickness.
c
        if (i.eq.itest .and. j.eq.jtest) then
          write(lp,'(a,i3)')
     &      'hybgen, deepest inflated layer too light   (stable):',k
          call flush(lp)
        endif !debug
c
        delsm=abs(saln(i,j,k-2,n)-saln(i,j,k-1,n))
        dels =abs(saln(i,j,k-1,n)-saln(i,j,k,  n))
        deltm=abs(temp(i,j,k-2,n)-temp(i,j,k-1,n))
        delt =abs(temp(i,j,k-1,n)-temp(i,j,k,  n))
c ---   sanity check on deltm and delsm
        q=min(temp(i,j,k-2,n),temp(i,j,k-1,n),temp(i,j,k,n))
        if     (q.gt. 6.0) then
          deltm=min( deltm,  6.0*(theta(i,j,k)-theta(i,j,k-1)) )
        elseif (q.gt. 0.0) then
          deltm=min( deltm, 10.0*(theta(i,j,k)-theta(i,j,k-1)) )
        else  !(q.le. 0.0)
          deltm=min( deltm, 25.0*(theta(i,j,k)-theta(i,j,k-1)) )
        endif
        delsm=min( delsm, 1.3*(theta(i,j,k)-theta(i,j,k-1)) )
        qts=0.0
        if     (delt.gt.epsil) then
          qts=max(qts, (min(deltm, 2.0*delt)-delt)/delt)  ! qts<=1.0
        endif
        if     (dels.gt.epsil) then
          qts=max(qts, (min(delsm, 2.0*dels)-dels)/dels)  ! qts<=1.0
        endif
        q=(theta(i,j,k)-th3d(i,j,k,  n))/
     &    (theta(i,j,k)-th3d(i,j,k-1,n))
        q=min(q,qts/(1.0+qts))  ! upper sublayer <= 50% of total
        q=qhrlx(k)*q
c ---   qhrlx is relaxation coefficient (inverse baroclinic time steps)
        p_hat=q*(p(i,j,k+1)-p(i,j,k))
        p(i,j,k)=p(i,j,k)+p_hat
        if     (hybflg.eq.0) then  !T&S
          temp(i,j,k,n)=temp(i,j,k,n)+(q/(1.0-q))*(temp(i,j,k,n)  -
     &                                             temp(i,j,k-1,n) )
          saln(i,j,k,n)=saln(i,j,k,n)+(q/(1.0-q))*(saln(i,j,k,n)  -
     &                                             saln(i,j,k-1,n) )
          th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
        elseif (hybflg.eq.1) then  !th&S
          th3d(i,j,k,n)=th3d(i,j,k,n)+(q/(1.0-q))*(th3d(i,j,k,n)  -
     &                                             th3d(i,j,k-1,n) )
          saln(i,j,k,n)=saln(i,j,k,n)+(q/(1.0-q))*(saln(i,j,k,n)  -
     &                                             saln(i,j,k-1,n) )
          temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
        elseif (hybflg.eq.2) then  !th&T
          th3d(i,j,k,n)=th3d(i,j,k,n)+(q/(1.0-q))*(th3d(i,j,k,n)  -
     &                                             th3d(i,j,k-1,n) )
          temp(i,j,k,n)=temp(i,j,k,n)+(q/(1.0-q))*(temp(i,j,k,n)  -
     &                                             temp(i,j,k-1,n) )
          saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
        endif
        if     (ntracr.gt.0 .and. p_hat.ne.0.0) then
          qtr=p_hat/(p(i,j,k)-p(i,j,k-1))  !ok because <1.0 and >0.0
          do ktr= 1,ntracr
            if     (trcflg(ktr).eq.2) then !temperature tracer
              tracer(i,j,k,n,ktr)=tracer(i,j,k,n,ktr)+
     &                           (q/(1.0-q))*(tracer(i,j,k,  n,ktr)-
     &                                        tracer(i,j,k-1,n,ktr))
            else !standard tracer - not split into two sub-layers
              tracer(i,j,k-1,n,ktr)=tracer(i,j,k-1,n,ktr)+
     &                                   qtr*(tracer(i,j,k,  n,ktr)-
     &                                        tracer(i,j,k-1,n,ktr))
cdiag              if (i.eq.itest .and. j.eq.jtest) then
cdiag                write(lp,'(a,i4,i3,5e12.3)')
cdiag     &            'hybgen, 10(+):',
cdiag     &            k,ktr,p_hat,p(i,j,k),p(i,j,k-1),
cdiag     &            qtr,tracer(i,j,k-1,n,ktr)
cdiag                call flush(lp)
cdiag              endif !debug
            endif !trcflg
          enddo !ktr
        endif !tracers
        if (mxlmy .and. p_hat.ne.0.0) then
          qtr=p_hat/(p(i,j,k)-p(i,j,k-1))  !ok because <1.0 and >0.0
          q2( i,j,k-1,n)=q2( i,j,k-1,n)+
     &                     qtr*(q2( i,j,k,n)-q2( i,j,k-1,n))
          q2l(i,j,k-1,n)=q2l(i,j,k-1,n)+
     &                     qtr*(q2l(i,j,k,n)-q2l(i,j,k-1,n))
              if (i.eq.itest .and. j.eq.jtest) then
               write(lp,'(a,i4,i3,6e12.3)')
     &            'hybgen, 10(+):',
     &            k,0,p_hat,p(i,j,k)-p(i,j,k-1),p(i,j,k+1)-p(i,j,k),
     &            qtr,q2(i,j,k-1,n),q2l(i,j,k-1,n)
               call flush(lp)
              endif !debug
        endif
        if (i.eq.itest .and. j.eq.jtest) then
          write(lp,'(a,i3,f6.3,5f8.3)')
     &      'hybgen, 10(+):',
     &      k,q,temp(i,j,k,n),saln(i,j,k,n),
     &          th3d(i,j,k,n)+thbase,theta(i,j,k)+thbase
          call flush(lp)
        endif !debug
        if (i.eq.itest .and. j.eq.jtest) then
          write(lp,'(a,i3,f6.3,5f8.3)')
     &      'hybgen, 10(-):',
     &      k,0.0,temp(i,j,k,n),saln(i,j,k,n),
     &          th3d(i,j,k,n)+thbase,theta(i,j,k)+thbase
          call flush(lp)
        endif !debug
      endif !too light
c
c --- massless or near-massless (thickness < dpthin) layers
c
      do k=kp+1,kk
        if (k.le.nhybrd) then
c ---     fill thin and massless layers on sea floor with fluid from above
          th3d(i,j,k,n)=th3d(i,j,k-1,n)
          saln(i,j,k,n)=saln(i,j,k-1,n)
          temp(i,j,k,n)=temp(i,j,k-1,n)
        elseif (th3d(i,j,k,n).ne.theta(i,j,k)) then
          if (hybflg.ne.2) then
c ---       fill with saln from above
            th3d(i,j,k,n)=max(theta(i,j,k), th3d(i,j,k-1,n))
            saln(i,j,k,n)=saln(i,j,k-1,n)
            temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
            saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
          else
c ---       fill with temp from above
            th3d(i,j,k,n)=max(theta(i,j,k), th3d(i,j,k-1,n))
            temp(i,j,k,n)=temp(i,j,k-1,n)
            saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
          endif
        endif
        do ktr= 1,ntracr
          tracer(i,j,k,n,ktr)=tracer(i,j,k-1,n,ktr)
        enddo
        if (mxlmy) then
          q2 (i,j,k,n)=q2( i,j,k-1,n)
          q2l(i,j,k,n)=q2l(i,j,k-1,n)
        endif
      enddo !k
c
c --- store one-dimensional arrays of -temp-, -saln-, and -p- for the 'old'
c --- vertical grid before attempting to restore isopycnic conditions
      pres(1)=p(i,j,1)
      do k=1,kk
        if     (hybflg.eq.0) then  !T&S
          s1d(k,1) = temp(i,j,k,n)
          s1d(k,2) = saln(i,j,k,n)
        elseif (hybflg.eq.1) then  !th&S
          s1d(k,1) = th3d(i,j,k,n)
          s1d(k,2) = saln(i,j,k,n)
        elseif (hybflg.eq.2) then  !th&T
          s1d(k,1) = th3d(i,j,k,n)
          s1d(k,2) = temp(i,j,k,n)
        endif
        do ktr= 1,ntracr
          s1d(k,2+ktr) = tracer(i,j,k,n,ktr)
        enddo
        if (mxlmy) then
          s1d(k,ntracr+3) = q2( i,j,k,n)
          s1d(k,ntracr+4) = q2l(i,j,k,n)
        endif
        pres(k+1)=p(i,j,k+1)
        dprs(k)  =pres(k+1)-pres(k)
        dpi( k)  =max(dprs(k),dpthin)
c
        if     (isopcm) then
          if     (k.le.fixlay) then
            lcm(k) = .false.  !fixed layers are never PCM
          else
c ---       thin and isopycnal layers remapped with PCM.
            lcm(k) = k.gt.nhybrd
     &               .or. dprs(k).le.dpthin
     &               .or. abs(th3d(i,j,k,n)-theta(i,j,k)).lt.hybiso
          endif !k<=fixlay:else
        endif !isopcm
      enddo !k
c
c --- try to restore isopycnic conditions by moving layer interfaces
c --- qhrlx(k) are relaxation coefficients (inverse baroclinic time steps)
c
      if (fixlay.ge.1) then
c
c ---   maintain constant thickness, layer k = 1
        k=1
        p_hat=p(i,j,k)+dp0ij(k)
        p(i,j,k+1)=min(p_hat,p(i,j,k+2))
      endif
c
      do k=2,nhybrd
c
        if (i.eq.itest .and. j.eq.jtest) then
          write(cinfo,'(a9,i2.2,1x)') '  do 88 k=',k
 109      format (i9,2i5,a,a/(i9,8x,a,a,i3,f9.2,f8.2,f9.2,f8.2))
          write(lp,109) nstep,itest+i0,jtest+j0,
     .      cinfo,':    othkns  odpth    nthkns  ndpth',
     .      (nstep,cinfo,':',ka,
     .      (pres(ka+1)-
     .       pres(ka)   )*qonem,
     .       pres(ka+1)  *qonem,
     .      (p(itest,jtest,ka+1)-
     .       p(itest,jtest,ka)   )*qonem,
     .       p(itest,jtest,ka+1)  *qonem,ka=1,kk)
          call flush(lp)
        endif !debug
c
        if (k.le.fixlay) then
c
c ---     maintain constant thickness, k <= fixlay
          if     (k.lt.kk) then  !p.kk+1 not changed
            p(i,j,k+1)=min(dp0cum(k+1),p(i,j,kk+1))
            if     (k.eq.fixlay) then
c ---         enforce interface order (may not be necessary).
              do ka= k+2,kk
                if     (p(i,j,ka).ge.p(i,j,k+1)) then
                  exit  ! usually get here quickly
                else
                  p(i,j,ka) = p(i,j,k+1)
                endif
              enddo !ka
            endif !k.eq.fixlay
          endif !k.lt.kk
c
          if (i.eq.itest .and. j.eq.jtest) then
            write(lp,'(a,i3.2,f8.2)') 'hybgen, fixlay :',
     &                                k+1,p(i,j,k+1)*qonem
            call flush(lp)
          endif !debug
        else
c
c ---     do not maintain constant thickness, k > fixlay
c
          if     (th3d(i,j,k,n).gt.theta(i,j,k)+epsil .and.
     &            k.gt.fixlay+1) then 
c
c ---       water in layer k is too dense
c ---       try to dilute with water from layer k-1
c ---       do not move interface if k = fixlay + 1
c
            if (th3d(i,j,k-1,n).ge.theta(i,j,k-1) .or.
     &          p(i,j,k).le.dp0cum(k)+onem .or.
     &          p(i,j,k+1)-p(i,j,k).le.p(i,j,k)-p(i,j,k-1)) then
c
c ---         if layer k-1 is too light, thicken the thinner of the two,
c ---         i.e. skip this layer if it is thicker.
c
              if (i.eq.itest .and. j.eq.jtest) then
                write(lp,'(a,3x,i2.2,1pe13.5)')
     &                'hybgen, too dense:',k,th3d(i,j,k,n)-theta(i,j,k)
              call flush(lp)
              endif !debug
c 
              if     ((theta(i,j,k)-th3d(i,j,k-1,n)).le.epsil) then
c               layer k-1 too dense, take entire layer
                p_hat =p(i,j,k-1)+dp0ij(k-1)
                p_hat0=p(i,j,k-1)-qqmn*dp0ij(k-1)
              else
                q=(theta(i,j,k)-th3d(i,j,k,  n))/
     &            (theta(i,j,k)-th3d(i,j,k-1,n))         ! -1 <= q < 0
                p_hat0=p(i,j,k)+q*(p(i,j,k+1)-p(i,j,k))  ! <p(i,j,k)
c               maintain minimum thickess of layer k-1.
                p_hat =p(i,j,k-1)+cushn(p_hat0-p(i,j,k-1),dp0ij(k-1))
              end if
                  if (i.eq.itest .and. j.eq.jtest) then
                    write(lp,'(a,i3.2,2f8.2)') 'hybgen,  p_hat(a)  :',
     &                    k,p_hat*qonem,(p(i,j,k-1)+dp0ij(k-1))*qonem
                    call flush(lp)
                  endif !debug
              p_hat=min(p_hat,p(i,j,kk+1))
                  if (i.eq.itest .and. j.eq.jtest) then
                    write(lp,'(a,i3.2,2f8.2)') 'hybgen,  p_hat(b)  :',
     &                    k,p_hat*qonem,(p(i,j,k-1)+dp0ij(k-1))*qonem
                    call flush(lp)
                  endif !debug
c
c ---         if isopycnic conditions cannot be achieved because of a blocking
c ---         layer in the interior ocean, move interface k-1 (and k-2 if
c ---         necessary) upward
c
              if     (k.le.fixlay+2) then
c ---           do nothing.
              else if (p_hat.ge.p(i,j,k) .and.
     &                 p(i,j,k-1).gt.dp0cum(k-1)+tenm .and.
     &                (p(i,j,kk+1)-p(i,j,k-1).lt.thkbop .or.
     &                 p(i,j,k-1) -p(i,j,k-2).gt.qqmx*dp0ij(k-2))) then ! k.gt.2
                p_hat2=p(i,j,k-2)+
     &                 cushn(p(i,j,k-1)-p_hat+p_hat0-p(i,j,k-2),
     &                       dp0ij(k-2))
                if (p_hat2.lt.p(i,j,k-1)-onemm) then
                  p(i,j,k-1)=(1.0-qhrlx(k-1))*p(i,j,k-1) +
     &                            qhrlx(k-1) *max(p_hat2,
     &                                    2.0*p(i,j,k-1)-p_hat)
                  if (i.eq.itest .and. j.eq.jtest) then
                    write(lp,'(a,i3.2,f8.2)') 'hybgen,  1blocking :',
     &                    k-1,p(i,j,k-1)*qonem
                    call flush(lp)
                  endif !debug
                  p_hat=p(i,j,k-1)+cushn(p_hat0-p(i,j,k-1),dp0ij(k-1))
                elseif (k.le.fixlay+3) then
c ---             do nothing.
                elseif (p(i,j,k-2).gt.dp0cum(k-2)+tenm .and.
     &                 (p(i,j,kk+1)-p(i,j,k-2).lt.thkbop .or.
     &                  p(i,j,k-2) -p(i,j,k-3).gt.qqmx*dp0ij(k-3))) then
                  p_hat3=p(i,j,k-3)+cushn(p(i,j,k-2)-p_hat+
     &                              p_hat0-p(i,j,k-3),
     &                              dp0ij(k-3))
                  if (p_hat3.lt.p(i,j,k-2)-onemm) then
                    p(i,j,k-2)=(1.0-qhrlx(k-2))*p(i,j,k-2) +
     &                              qhrlx(k-2)*max(p_hat3,
     &                                      2.0*p(i,j,k-2)-p(i,j,k-1))
                    if (i.eq.itest .and. j.eq.jtest) then
                      write(lp,'(a,i3.2,f8.2)') 'hybgen,  2blocking :',
     &                      k-2,p(i,j,k-2)*qonem
                      call flush(lp)
                    endif !debug
                    p_hat2=p(i,j,k-2)+cushn(p(i,j,k-1)-p_hat+
     &                                      p_hat0-p(i,j,k-2),
     &                                      dp0ij(k-2))
                    if (p_hat2.lt.p(i,j,k-1)-onemm) then
                      p(i,j,k-1)=(1.0-qhrlx(k-1))*p(i,j,k-1) +
     &                                qhrlx(k-1) *max(p_hat2,
     &                                        2.0*p(i,j,k-1)-p_hat)
                      if (i.eq.itest .and. j.eq.jtest) then
                        write(lp,'(a,i3.2,f8.2)')
     &                             'hybgen,  3blocking :',
     &                             k-1,p(i,j,k-1)*qonem
                        call flush(lp)
                      endif !debug
                      p_hat=p(i,j,k-1)+cushn(p_hat0-p(i,j,k-1),
     &                                       dp0ij(k-1))
                    endif !p_hat2
                  endif !p_hat3
                endif !p_hat2:blocking
              endif !blocking
                  if (i.eq.itest .and. j.eq.jtest) then
                    write(lp,'(a,i3.2,2f8.2)') 'hybgen,  p_hat(c)  :',
     &                    k,p_hat*qonem,(p(i,j,k-1)+dp0ij(k-1))*qonem
                    call flush(lp)
                  endif !debug
              p_hat=min(p_hat,p(i,j,kk+1))
                  if (i.eq.itest .and. j.eq.jtest) then
                    write(lp,'(a,i3.2,2f8.2)') 'hybgen,  p_hat(d)  :',
     &                    k,p_hat*qonem,(p(i,j,k-1)+dp0ij(k-1))*qonem
                    call flush(lp)
                  endif !debug
c
c ---         move upper interface up or down
              p(i,j,k)=min( (1.0-qhrlx(k))*p(i,j,k) +
     &                           qhrlx(k) *p_hat,
     &                      p(i,j,k+1) )
c
              if (i.eq.itest .and. j.eq.jtest) then
                write(lp,'(a,i3.2,f8.2)') 'hybgen, entrain(k) :',
     &                                    k,p(i,j,k)*qonem
                call flush(lp)
              endif !debug
c
            endif  !too-dense adjustment
c
          elseif (th3d(i,j,k,n).lt.theta(i,j,k)-epsil) then   ! layer too light
c
c ---       water in layer k is too light
c ---       try to dilute with water from layer k+1
c ---       do not entrain if layer k touches bottom
c
            if (p(i,j,k+1).lt.p(i,j,kk+1)) then  ! k<kk
              if (th3d(i,j,k+1,n).le.theta(i,j,k+1) .or.
     &            p(i,j,k+1).le.dp0cum(k+1)+onem    .or.
     &            p(i,j,k+1)-p(i,j,k).lt.p(i,j,k+2)-p(i,j,k+1)) then
c
c ---           if layer k+1 is too dense, thicken the thinner of the two,
c ---           i.e. skip this layer if it is not thinner than the other.
c
                if (i.eq.itest .and. j.eq.jtest) then
                  write(lp,'(a,3x,i2.2,1pe13.5)')
     &                 'hybgen, too light:',k,
     &                  theta(i,j,k)-th3d(i,j,k,n)
                  call flush(lp)
                endif !debug
c
                if     ((th3d(i,j,k+1,n)-theta(i,j,k)).le.epsil) then
c                 layer k-1 too light, take entire layer
                  p_hat=p(i,j,k+2)
                else
                  q=(th3d(i,j,k,  n)-theta(i,j,k))/
     &              (th3d(i,j,k+1,n)-theta(i,j,k))          !-1 <= q < 0
                  p_hat=p(i,j,k+1)+q*(p(i,j,k)-p(i,j,k+1))  !>p(i,j,k+1)
                endif
c
c ---           if layer k+1 does not touch the bottom then maintain minimum
c ---           thicknesses of layers k and k+1 as much as possible,
c ---           but permit layers to collapse to zero thickness at the bottom
c
                if     (p(i,j,k+2).lt.p(i,j,kk+1)) then
                  if     (p(i,j,kk+1)-p(i,j,k).gt.
     &                    dp0ij(k)+dp0ij(k+1)     ) then
                    p_hat=p(i,j,k+2)-cushn(p(i,j,k+2)-p_hat,dp0ij(k+1))
                  endif
                  p_hat=p(i,j,k)+max(p_hat-p(i,j,k),dp0ij(k))
                  p_hat=min(p_hat,
     &                      max(0.5*(p(i,j,k+1)+p(i,j,k+2)),
     &                               p(i,j,k+2)-dp0ij(k+1)))
                else
                  p_hat=min(p(i,j,k+2),p_hat)
                endif !p.k+2<p.kk+1
                if (p_hat.gt.p(i,j,k+1)) then
c ---             entrain layer k+1 water into layer k.
                  p(i,j,k+1)=(1.0-qhrlx(k+1))*p(i,j,k+1) +
     &                            qhrlx(k+1) *p_hat
                endif !entrain
c
                if (i.eq.itest .and. j.eq.jtest) then
                  write(lp,'(a,i3.2,f8.2)')
     &                 'hybgen, entrain(k+):',k,p(i,j,k+1)*qonem
                  call flush(lp)
                endif !debug
c
              endif !too-light adjustment
            endif !above bottom
          endif !too dense or too light
c
c ---     if layer above is still too thin, move interface down.
          p_hat0=min(p(i,j,k-1)+dp0ij(k-1),p(i,j,kk+1))
          if (p_hat0.gt.p(i,j,k)) then
            p_hat =(1.0-qhrlx(k-1))*p(i,j,k)+
     &                  qhrlx(k-1) *p_hat0
            p(i,j,k)=min(p_hat,p(i,j,k+1))
c
            if (i.eq.itest .and. j.eq.jtest) then
              write(lp,'(a,i3.2,f8.2)')
     &             'hybgen, min. thknss (k+):',k,p(i,j,k+1)*qonem
              call flush(lp)
            endif !debug
          endif
c
        endif !k.le.fixlay:else
c
      enddo !k  vertical coordinate relocation
c
c --- enforce layer thickness ratios
c
      if     (hybthk.gt.0.0) then
        p(i,j,2) = max( p(i,j,2),  p(i,j,1) )
        p(i,j,3) = max( p(i,j,3),  p(i,j,2) )
        do k= 2,kk-2
          p(i,j,k+2) = max( p(i,j,k+2),  p(i,j,k+1) )
          if     (th3d(i,j,k-1,n).lt.theta(i,j,k)-epsil .and.
     &            th3d(i,j,k+1,n).gt.theta(i,j,k)+epsil      ) then
            dptop = p(i,j,k)   - p(i,j,k-1)
            dpmid = p(i,j,k+1) - p(i,j,k)
            dpbot = p(i,j,k+2) - p(i,j,k+1)
            if     (p(i,j,k+2)+dpmid.lt.p(i,j,kk+1) .and.
     &              dpmid.lt.dp0k(k)                .and.
     &              dpmid.lt.hybthk*min(dptop,dpbot)     ) then
c
c ---         thicken layer k by taking fluid from above and below
c ---         in proportion so that the net is at k's target density
c
              q=(theta(i,j,k)    -th3d(i,j,k-1,n))/
     &          (th3d( i,j,k+1,n)-th3d(i,j,k-1,n))         ! 0<q<1
              dpinc = min(       dp0k(k)-dpmid  ,
     &                     (hybthk*dptop-dpmid)/
     &                     (1.0+(1.0-q)*hybthk) ,
     &                     (hybthk*dpbot-dpmid)/
     &                     (1.0+     q *hybthk) )
              dpinc = dpinc*min(qhrlx(k),qhrlx(k+1))
              p(i,j,k)   = p(i,j,k)   - (1.0-q)*dpinc
              p(i,j,k+1) = p(i,j,k+1) +      q *dpinc
c
              if (i.eq.itest .and. j.eq.jtest) then
                write(lp,'(a,i3.2,f10.3)')
     &               'hybgen, dptop:',k,dptop*qonem
                write(lp,'(a,i3.2,f10.3)')
     &               'hybgen, dpmid:',k,dpmid*qonem
                write(lp,'(a,i3.2,f10.3)')
     &               'hybgen, dpbot:',k,dpbot*qonem
                write(lp,'(a,i3.2,f10.3)')
     &               'hybgen, dpinc:',k,dpinc*qonem
                write(lp,'(a,i3.2,f10.3)')
     &               'hybgen,     q:',k,q
                write(lp,'(a,i3.2,f10.3)')
     &               'hybgen, qhrlx:',k,min(qhrlx(k),qhrlx(k+1))
                call flush(lp)
              endif !debug
            endif
          endif !stable density
        enddo !k
      endif !hybthk
c
c --- remap scalar field profiles from the 'old' vertical
c --- grid onto the 'new' vertical grid.
c
      if     (lconserve) then  !usually .false.
        do ktr=1,nums1d
          asum(ktr,1) = 0.d0
          asum(ktr,2) = 0.d0
          asum(ktr,3) = 0.d0
        enddo !ktr
      endif
c
      prsf(1) = p(i,j,1)
      do k=1,kk
        dp(i,j,k,n) = max( p(i,j,k+1)-prsf(k), 0.0 )  !enforce interface order
c ---   to avoid roundoff errors in -dpudpv- after restart, make sure p=p(dp)
        prsf(k+1)   = prsf(k) + dp(i,j,k,n)
        p(i,j,k+1)  = prsf(k+1)
      enddo
      if     (hybmap.eq.0) then !PCM
        call hybgen_pcm_remap(s1d,pres,dprs,
     &                        f1d,prsf,kk,kk,nums1d,dpthin)
      elseif (hybmap.eq.1) then !PLM (as in 2.1.08)
        call hybgen_plm_coefs(s1d,     dprs,lcm,c1d,
     &                                 kk,   nums1d,dpthin)
        call hybgen_plm_remap(s1d,pres,dprs,    c1d,
     &                        f1d,prsf,kk,kk,nums1d,dpthin)
      elseif (hybmap.eq.2) then !PPM
        call hybgen_ppm_coefs(s1d,     dpi, lcm,c1d,
     &                                 kk,   nums1d,dpthin)
        call hybgen_ppm_remap(s1d,pres,dprs,    c1d,
     &                        f1d,prsf,kk,kk,nums1d,dpthin)
      elseif (hybmap.eq.3) then !WENO-like
        call hybgen_weno_coefs(s1d,     dpi, lcm,c1d,
     &                                  kk,   nums1d,dpthin)
        call hybgen_weno_remap(s1d,pres,dprs,    c1d,
     &                         f1d,prsf,kk,kk,nums1d,dpthin)
      endif
      do k=1,kk
        if     (hybflg.eq.0) then  !T&S
          temp(i,j,k,n) = f1d(k,1)
          saln(i,j,k,n) = f1d(k,2)
          th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
        elseif (hybflg.eq.1) then  !th&S
          th3d(i,j,k,n) = f1d(k,1)
          saln(i,j,k,n) = f1d(k,2)
          temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,
     &                         saln(i,j,k,n))
*         saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,
*    &                         temp(i,j,k,n))
        elseif (hybflg.eq.2) then  !th&T
          th3d(i,j,k,n) = f1d(k,1)
          temp(i,j,k,n) = f1d(k,2)
          saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,
     &                         temp(i,j,k,n))
        endif
        do ktr= 1,ntracr
          tracer(i,j,k,n,ktr) = f1d(k,2+ktr)
        enddo
        if (mxlmy) then
          q2( i,j,k,n) = f1d(k,ntracr+3)
          q2l(i,j,k,n) = f1d(k,ntracr+4)
        endif
c
        if     (lconserve) then  !usually .false.
          zthk = dp(i,j,k,n)
          do ktr= 1,nums1d
            asum(ktr,1) = asum(ktr,1) + s1d(k,ktr)*dprs(k)
            asum(ktr,2) = asum(ktr,2) + f1d(k,ktr)*zthk
          enddo
        endif !lconserve
c
      enddo !k
c
cdiag if (i.eq.itest .and. j.eq.jtest) then
*       write (lp,'(i9,3a/(i9,3f23.17))')
*    &  nstep,
*    &  '                   dens',
*    &  '                  thkns',
*    &  '                   dpth',
*    &  (k,tthe(k,1),    dprs(k)*qonem,    pres(k+1)*qonem,
*    &   k,th3d(i,j,k,n),dp(i,j,k,n)*qonem,p(i,j,k+1)*qonem,
*    &  k=1,kk)
cdiag   write (lp,'(i9,3a/(i9,3f23.17))')
cdiag&  nstep,
cdiag&  '               tracer.1',
cdiag&  '                  thkns',
cdiag&  '                   dpth',
cdiag&  (k,ttrc(      k,1,1),  dprs(k)*qonem,   pres(k+1)*qonem,
cdiag&   k,tracer(i,j,k,n,1),dp(i,j,k,n)*qonem,p(i,j,k+1)*qonem,
cdiag&  k=1,kk)
cdiag   call flush(lp)
cdiag endif !debug
c
      if     (lconserve) then  !usually .false.
c
c ---   enforce water column conservation
c
        do ktr=1,nums1d
          q = asum(ktr,1)-asum(ktr,2)
          if     (q.eq.0.0) then
            offset(ktr) = 0.0
          elseif (abs(asum(ktr,2)).lt.2.0*abs(q)) then
            offset(ktr) = sign(zp5,q*asum(ktr,2))  !        -0.5 or  +0.5
          else
            offset(ktr) =          q/asum(ktr,2)   !between -0.5 and +0.5
          endif
        enddo !ktr
        do k=1,kk
          if     (hybflg.eq.0) then  !T&S
            temp(i,j,k,n)=temp(i,j,k,n)*(1.0+offset(1))
            saln(i,j,k,n)=saln(i,j,k,n)*(1.0+offset(2))
            th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
          elseif (hybflg.eq.1) then  !th&S
            th3d(i,j,k,n)=th3d(i,j,k,n)*(1.0+offset(1))
            saln(i,j,k,n)=saln(i,j,k,n)*(1.0+offset(2))
            temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,
     &                           saln(i,j,k,n))
          elseif (hybflg.eq.2) then  !th&T
            th3d(i,j,k,n)=th3d(i,j,k,n)*(1.0+offset(1))
            temp(i,j,k,n)=temp(i,j,k,n)*(1.0+offset(2))
            saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,
     &                           temp(i,j,k,n))
          endif
          do ktr= 1,ntracr
            tracer(i,j,k,n,ktr)=tracer(i,j,k,n,ktr)*(1.0+offset(ktr+2))
          enddo !ktr
          if (mxlmy) then
            q2( i,j,k,n)=q2( i,j,k,n)*(1.0+offset(ntracr+3))
            q2l(i,j,k,n)=q2l(i,j,k,n)*(1.0+offset(ntracr+4))
          endif
c
          if     (.false.) then !debugging
            zthk = dp(i,j,k,n)
            if     (hybflg.eq.0) then  !T&S
              asum(1,3) = asum(1,3) + temp(i,j,k,n)*zthk
              asum(2,3) = asum(2,3) + saln(i,j,k,n)*zthk
            elseif (hybflg.eq.1) then  !th&S
              asum(1,3) = asum(1,3) + th3d(i,j,k,n)*zthk
              asum(2,3) = asum(2,3) + saln(i,j,k,n)*zthk
            elseif (hybflg.eq.2) then  !th&T
              asum(1,3) = asum(1,3) + th3d(i,j,k,n)*zthk
              asum(2,3) = asum(2,3) + temp(i,j,k,n)*zthk
            endif
            do ktr= 1,ntracr
              asum(ktr+2,3) = asum(ktr+2,3) + tracer(i,j,k,n,ktr)*zthk
            enddo !ktr
            if (mxlmy) then
              asum(ntracr+3,3) = asum(ntracr+3,3) +  q2( i,j,k,n)*zthk
              asum(ntracr+4,3) = asum(ntracr+4,3) +  q2l(i,j,k,n)*zthk
            endif
          endif !debuging
        enddo !k
c
        if     (.false. .and. !debugging
     &          i.eq.itest .and. j.eq.jtest) then
          do ktr= 1,nums1d
            write(lp,'(a,1p4e16.8,i3)')
     &        'hybgen,sum:',
     &        asum(ktr,1)/p(i,j,kk+1),
     &        asum(ktr,2)/p(i,j,kk+1),
     &        asum(ktr,3)/p(i,j,kk+1),
     &        offset(ktr),ktr
          enddo !ktr
        endif !debugging .and. i.eq.itest .and. j.eq.jtest
        if     (.false. .and. !debugging
     &          j.eq.jtest) then
          ktr=1
*         if     (abs(offset(ktr)).gt.1.e-08) then
          if     (abs(offset(ktr)).gt.1.e-12) then
            write(lp,'(a,1p4e16.8,i3)')
     &        'hybgen,sum:',
     &        asum(ktr,1)/p(i,j,kk+1),
     &        asum(ktr,2)/p(i,j,kk+1),
     &        asum(ktr,3)/p(i,j,kk+1),
     &        offset(ktr),i
          endif !large offset
        endif !debugging .and. j.eq.jtest
      endif !lconserve
c
cdiag if (i.eq.itest .and. j.eq.jtest) then
*       write (lp,'(i9,3a/(i9,3f23.17))')
*    &  nstep,
*    &  '                   dens',
*    &  '                  thkns',
*    &  '                   dpth',
*    &  (k,tthe(k,1),    dprs(k)*qonem,    pres(k+1)*qonem,
*    &   k,th3d(i,j,k,n),dp(i,j,k,n)*qonem,p(i,j,k+1)*qonem,
*    &  k=1,kk)
cdiag   write (lp,'(i9,3a/(i9,3f23.17))')
cdiag&  nstep,
cdiag&  '               tracer.1',
cdiag&  '                  thkns',
cdiag&  '                   dpth',
cdiag&  (k,ttrc(      k,1,1),  dprs(k)*qonem,   pres(k+1)*qonem,
cdiag&   k,tracer(i,j,k,n,1),dp(i,j,k,n)*qonem,p(i,j,k+1)*qonem,
cdiag&  k=1,kk)
cdiag   call flush(lp)
cdiag endif
c
 103  format (i9,2i5,a/(33x,i3,2f8.3,f8.3,f9.3,f9.2))
      if (i.eq.itest .and. j.eq.jtest) then
       if     (hybflg.eq.0) then  !T&S
        write (lp,103) nstep,itest+i0,jtest+j0,
     &  '    hybgen, do 22:  temp    saln    dens     thkns     dpth',
     &  (k,s1d(k,1),s1d(k,2),0.0,
     &   (pres(k+1)-pres(k))*qonem,pres(k+1)*qonem,
     &   k,temp(i,j,k,n),saln(i,j,k,n),
     &   th3d(i,j,k,n)+thbase,dp(i,j,k,n)*qonem,
     &   p(i,j,k+1)*qonem,
     &  k=1,kk)
       elseif (hybflg.eq.1) then  !th&S
        write (lp,103) nstep,itest+i0,jtest+j0,
     &  '    hybgen, do 22:  temp    saln    dens     thkns     dpth',
     &  (k,0.0,s1d(k,2),s1d(k,1)+thbase,
     &   (pres(k+1)-pres(k))*qonem,pres(k+1)*qonem,
     &   k,temp(i,j,k,n),saln(i,j,k,n),
     &   th3d(i,j,k,n)+thbase,dp(i,j,k,n)*qonem,
     &   p(i,j,k+1)*qonem,
     &  k=1,kk)
       elseif (hybflg.eq.2) then  !th&T
        write (lp,103) nstep,itest+i0,jtest+j0,
     &  '    hybgen, do 22:  temp    saln    dens     thkns     dpth',
     &  (k,s1d(k,2),0.0,s1d(k,1)+thbase,
     &   (pres(k+1)-pres(k))*qonem,pres(k+1)*qonem,
     &   k,temp(i,j,k,n),saln(i,j,k,n),
     &   th3d(i,j,k,n)+thbase,dp(i,j,k,n)*qonem,
     &   p(i,j,k+1)*qonem,
     &  k=1,kk)
       endif
       call flush(lp)
      endif !debug
c
*     enddo !i
*     enddo !l
c
      return
      end subroutine hybgen

      subroutine hybgen_pcm_remap(si,pi,dpi,
     &                            so,po,ki,ko,ks,thin)
      implicit none
c
      integer ki,ko,ks
      real    si(ki,ks),pi(ki+1),dpi(ki),
     &        so(ko,ks),po(ko+1),thin
c
c-----------------------------------------------------------------------
c  1) remap from one set of vertical cells to another.
c     method: piecewise constant across each input cell
c             the output is the average of the interpolation
c             profile across each output cell.
c
c  2) input arguments:
c       si    - initial scalar fields in pi-layer space
c       pi    - initial layer interface depths (non-negative)
c                  pi(   1) is the surface
c                  pi(ki+1) is the bathymetry
c                  pi(k+1) >= pi(k)
c       dpi   - initial layer thicknesses (dpi(k)=pi(k+1)-pi(k))
c       ki    - number of  input layers
c       ko    - number of output layers
c       ks    - number of fields
c       po    - target interface depths (non-negative)
c                  po(   1) is the surface
c                  po(ko+1) is the bathymetry (== pi(ki+1))
c                  po(k+1) >= po(k)
c       thin  - layer thickness (>0) that can be ignored
c
c  3) output arguments:
c       so    - scalar fields in po-layer space
c
c  4) Tim Campbell, Mississippi State University, October 2002.
c     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
c-----------------------------------------------------------------------
c
      integer i,k,l,lb,lt
      real    dpb,dpt,xb,xt,zb,zt,zx,o
      real*8  sz
c
      zx=pi(ki+1) !maximum depth
      zb=max(po(1),pi(1))
      lb=1
      do while (pi(lb+1).lt.zb .and. lb.lt.ki)
        lb=lb+1
      enddo
      do k= 1,ko  !output layers
        zt = zb
        zb = min(po(k+1),zx)
*       write(lp,*) 'k,zt,zb = ',k,zt,zb
        lt=lb !top will always correspond to bottom of previous
              !find input layer containing bottom output interface
        do while (pi(lb+1).lt.zb .and. lb.lt.ki)
          lb=lb+1
        enddo
        if     (zb-zt.le.thin .or. zt.ge.zx) then
          if     (k.ne.1) then
c
c ---       thin or bottomed layer, values taken from layer above
c
            do i= 1,ks
              so(k,i) = so(k-1,i)
            enddo !i
          else !thin surface layer
            do i= 1,ks
              so(k,i) = si(k,i)
            enddo !i
          endif
        else
c
c         form layer averages.
c         use PPM-like logic (may not have minimum operation count)
c
*         if     (pi(lb).gt.zt) then
*           write(lp,*) 'bad lb = ',lb
*           stop
*         endif
          if     (lt.ne.lb) then  !multiple layers
            xt=(zt-pi(lt))/max(dpi(lt),thin)
            xb=(zb-pi(lb))/max(dpi(lb),thin)
            dpt=pi(lt+1)-zt
            dpb=zb-pi(lb)
            do i= 1,ks
              o  = si((lt+lb)/2,i)  !offset to reduce round-off
              sz = dpt*(si(lt,i)-o)
              do l=lt+1,lb-1
                sz = sz+dpi(l)*(si(l,i)-o)
              enddo !l
              sz = sz+dpb*(si(lb,i)-o)
              so(k,i) = o + sz/(zb-zt)  !zb-zt>=thin
            enddo !i
          else  !single layer
            do i= 1,ks
              so(k,i) = si(lt,i)
            enddo !i
          endif
        endif !thin:std layer
      enddo !k
      return
      end subroutine hybgen_pcm_remap

      subroutine hybgen_plm_coefs(si,dpi,lc,ci,kk,ks,thin)
      implicit none
c
      integer kk,ks
      logical lc(kk)
      real    si(kk,ks),dpi(kk),ci(kk,ks),thin
c
c-----------------------------------------------------------------------
c  1) coefficents for remaping from one set of vertical cells to another.
c     method: piecewise linear across each input cell with
c             monotonized central-difference limiter.
c             the output is the average of the interpolation
c             profile across each output cell.
c
c  2) input arguments:
c       si    - initial scalar fields in pi-layer space
c       dpi   - initial layer thicknesses (dpi(k)=pi(k+1)-pi(k))
c       lc    - use PCM for selected layers
c       kk    - number of layers
c       ks    - number of fields
c       thin  - layer thickness (>0) that can be ignored
c
c  3) output arguments:
c       ci    - coefficents (slopes) for hybgen_plm_remap
c                profile(y)=si+ci*(y-1),  0<=y<=1
c
c  4) Tim Campbell, Mississippi State University, October 2002.
c     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
c-----------------------------------------------------------------------
c
      integer k,i
      real    qcen,zbot,zcen,ztop
c
      do i= 1,ks
        ci(1, i) = 0.0
        ci(kk,i) = 0.0
      enddo !i
      do k= 2,kk-1
        if     (lc(k) .or. dpi(k).le.thin) then  !use PCM
          do i= 1,ks
            ci(k,i) = 0.0
          enddo !i
        else
c ---     use qcen in place of 0.5 to allow for non-uniform grid
          qcen = dpi(k)/(dpi(k)+0.5*(dpi(k-1)+dpi(k+1)))  !dpi(k)>thin
          do i= 1,ks
c ---       PLM (non-zero slope, but no new extrema)
c ---       layer value is si-0.5*ci at top    interface,
c ---                  and si+0.5*ci at bottom interface.
c
c ---       monotonized central-difference limiter (van Leer, 1977,
c ---       JCP 23 pp 276-299).  For a discussion of PLM limiters, see
c ---       Finite Volume Methods for Hyperbolic Problems by R.J. Leveque.
            ztop = 2.0*(si(k,  i)-si(k-1,i))
            zbot = 2.0*(si(k+1,i)-si(k,  i))
            zcen =qcen*(si(k+1,i)-si(k-1,i))
            if     (ztop*zbot.gt.0.0) then !ztop,zbot are the same sign
              ci(k,i)=sign(min(abs(zcen),abs(zbot),abs(ztop)),zbot)
            else
              ci(k,i)=0.0  !local extrema, so no slope
            endif
          enddo !i
        endif  !PCM:PLM
      enddo !k
      return
      end subroutine hybgen_plm_coefs

      subroutine hybgen_plm_remap(si,pi,dpi,ci,
     &                            so,po,ki,ko,ks,thin)
      implicit none
c
      integer ki,ko,ks
      real    si(ki,ks),pi(ki+1),dpi(ki),ci(ki,ks),
     &        so(ko,ks),po(ko+1),thin
c
c-----------------------------------------------------------------------
c  1) remap from one set of vertical cells to another.
c     method: piecewise linear across each input cell
c             the output is the average of the interpolation
c             profile across each output cell.
c
c  2) input arguments:
c       si    - initial scalar fields in pi-layer space
c       pi    - initial layer interface depths (non-negative)
c                  pi(   1) is the surface
c                  pi(ki+1) is the bathymetry
c                  pi(k+1) >= pi(k)
c       dpi   - initial layer thicknesses (dpi(k)=pi(k+1)-pi(k))
c       ci    - coefficents (slopes) from hybgen_plm_coefs
c                profile(y)=si+ci*(y-1),  0<=y<=1
c       ki    - number of  input layers
c       ko    - number of output layers
c       ks    - number of fields
c       po    - target interface depths (non-negative)
c                  po(   1) is the surface
c                  po(ko+1) is the bathymetry (== pi(ki+1))
c                  po(k+1) >= po(k)
c       thin  - layer thickness (>0) that can be ignored
c
c  3) output arguments:
c       so    - scalar fields in po-layer space
c
c  4) Tim Campbell, Mississippi State University, October 2002.
c     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
c-----------------------------------------------------------------------
c
      integer i,k,l,lb,lt
      real    c0,qb0,qb1,qt0,qt1,xb,xt,zb,zt,zx,o
      real*8  sz
c
      zx=pi(ki+1) !maximum depth
      zb=max(po(1),pi(1))
      lb=1
      do while (pi(lb+1).lt.zb .and. lb.lt.ki)
        lb=lb+1
      enddo
      do k= 1,ko  !output layers
        zt = zb
        zb = min(po(k+1),zx)
*       write(lp,*) 'k,zt,zb = ',k,zt,zb
        lt=lb !top will always correspond to bottom of previous
              !find input layer containing bottom output interface
        do while (pi(lb+1).lt.zb .and. lb.lt.ki)
          lb=lb+1
        enddo
        if     (zb-zt.le.thin .or. zt.ge.zx) then
          if     (k.ne.1) then
c
c ---       thin or bottomed layer, values taken from layer above
c
            do i= 1,ks
              so(k,i) = so(k-1,i)
            enddo !i
          else !thin surface layer
            do i= 1,ks
              so(k,i) = si(k,i)
            enddo !i
          endif
        else
c
c         form layer averages.
c         use PPM-like logic (may not have minimum operation count)
c
*         if     (pi(lb).gt.zt) then
*           write(lp,*) 'bad lb = ',lb
*           stop
*         endif
          xt=(zt-pi(lt))/max(dpi(lt),thin)
          xb=(zb-pi(lb))/max(dpi(lb),thin)
          if     (lt.ne.lb) then  !multiple layers
            qt0 = (1.0-xt)
            qt1 = (1.0-xt**2)*0.5
            qb0 =      xb
            qb1 =      xb**2 *0.5
            do i= 1,ks
              o  = si((lt+lb)/2,i)  !offset to reduce round-off
              c0 = si(lt,i) - o - 0.5*ci(lt,i)
              sz=  dpi(lt)*(c0*qt0 + ci(lt,i)*qt1)
              do l=lt+1,lb-1
                sz = sz+dpi(l)*(si(l,i) - o)
              enddo !l
              c0 = si(lb,i) - o - 0.5*ci(lb,i)
              sz = sz+dpi(lb)*(c0*qb0 + ci(lb,i)*qb1)
              so(k,i) = o + sz/(zb-zt)  !zb-zt>=thin
            enddo !i
          else  !single layer
            qt1 = (xb**2-xt**2 - (xb-xt))*0.5
            do i= 1,ks
              sz = dpi(lt)*(ci(lt,i)*qt1)
              so(k,i) = si(lt,i) + sz/(zb-zt)  !zb-zt>=thin
            enddo !i
          endif
        endif !thin:std layer
      enddo !k
      return
      end subroutine hybgen_plm_remap

      subroutine hybgen_ppm_coefs(s,dp,lc,ci,kk,ks,thin)
      implicit none
c
      integer kk,ks
      logical lc(kk)
      real    s(kk,ks),dp(kk),ci(kk,ks,3),thin
c
c-----------------------------------------------------------------------
c  1) coefficents for remaping from one set of vertical cells to another.
c     method: monotonic piecewise parabolic across each input cell
c             the output is the average of the interpolation
c             profile across each output cell.
c
c     Colella, P. & P.R. Woodward, 1984, J. Comp. Phys., 54, 174-201.
c
c  2) input arguments:
c       s     - initial scalar fields in pi-layer space
c       dp    - initial layer thicknesses (>=thin)
c       lc    - use PCM for selected layers
c       kk    - number of layers
c       ks    - number of fields
c       thin  - layer thickness (>0) that can be ignored
c
c  3) output arguments:
c       ci    - coefficents for hybgen_ppm_remap
c                profile(y)=ci.1+ci.2*y+ci.3*y^2, 0<=y<=1
c
c  4) Tim Campbell, Mississippi State University, October 2002.
c     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
c-----------------------------------------------------------------------
c
      integer j,i
      real    da,a6,slj,scj,srj
      real    as(kk),al(kk),ar(kk)
      real     dpjp(kk), dp2jp(kk), dpj2p(kk),
     &        qdpjp(kk),qdp2jp(kk),qdpj2p(kk),dpq3(kk),qdp4(kk)
c
      !compute grid metrics
      do j=1,kk-1
         dpjp( j) = dp(j)   + dp(j+1)
         dp2jp(j) = dp(j)   + dpjp(j)
         dpj2p(j) = dpjp(j) + dp(j+1)
        qdpjp( j) = 1.0/dpjp( j)
        qdp2jp(j) = 1.0/dp2jp(j)
        qdpj2p(j) = 1.0/dpj2p(j)
      enddo !j
         dpq3(2) = dp(2)/(dp(1)+dpjp(2))
      do j=3,kk-1
         dpq3(j) = dp(j)/(dp(j-1)+dpjp(j)) !dp(j)/      (dp(j-1)+dp(j)+dp(j+1))
         qdp4(j) = 1.0/(dpjp(j-2)+dpjp(j)) !1.0/(dp(j-2)+dp(j-1)+dp(j)+dp(j+1))
      enddo !j
c
      do i= 1,ks
        !Compute average slopes: Colella, Eq. (1.8)
        as(1)=0.
        do j=2,kk-1
          if     (lc(j) .or. dp(j).le.thin) then  !use PCM
            as(j) = 0.0
          else
            slj=s(j,  i)-s(j-1,i)
            srj=s(j+1,i)-s(j,  i)
            if (slj*srj.gt.0.) then
              scj=dpq3(j)*( dp2jp(j-1)*srj*qdpjp(j)
     &                     +dpj2p(j)  *slj*qdpjp(j-1) )
              as(j)=sign(min(abs(2.0*slj),abs(scj),abs(2.0*srj)),scj)
            else
              as(j)=0.
            endif
          endif  !PCM:PPM
        enddo !j
        as(kk)=0.
        !Compute "first guess" edge values: Colella, Eq. (1.6)
        al(1)=s(1,i)  !1st layer PCM
        ar(1)=s(1,i)  !1st layer PCM
        al(2)=s(1,i)  !1st layer PCM
        do j=3,kk-1
          al(j)=s(j-1,i)+dp(j-1)*(s(j,i)-s(j-1,i))*qdpjp(j-1)
     &         +qdp4(j)*(
     &            2.*dp(j)*dp(j-1)*qdpjp(j-1)*(s(j,i)-s(j-1,i))*
     &            ( dpjp(j-2)*qdp2jp(j-1)
     &             -dpjp(j)  *qdpj2p(j-1) )
     &            -dp(j-1)*as(j)  *dpjp(j-2)*qdp2jp(j-1)
     &            +dp(j)  *as(j-1)*dpjp(j)  *qdpj2p(j-1)
     &              )
          ar(j-1)=al(j)
        enddo !j
        ar(kk-1)=s(kk,i)  !last layer PCM
        al(kk)  =s(kk,i)  !last layer PCM
        ar(kk)  =s(kk,i)  !last layer PCM
        !Impose monotonicity: Colella, Eq. (1.10)
        do j=2,kk-1
          if     (lc(j) .or. dp(j).le.thin) then  !use PCM
            al(j)=s(j,i)
            ar(j)=s(j,i)
          elseif ((s(j+1,i)-s(j,i))*(s(j,i)-s(j-1,i)).le.0.) then !local extremum
            al(j)=s(j,i)
            ar(j)=s(j,i)
          else
            da=ar(j)-al(j)
            a6=6.0*s(j,i)-3.0*(al(j)+ar(j))
            if     (da*a6 .gt.  da*da) then !peak in right half of zone
              al(j)=3.0*s(j,i)-2.0*ar(j)
            elseif (da*a6 .lt. -da*da) then !peak in left half of zone
              ar(j)=3.0*s(j,i)-2.0*al(j)
            endif
          endif
        enddo !j
        !Set coefficients
        do j=1,kk
          if     (al(j).ne.ar(j)) then
            ci(j,i,1)=al(j)
            ci(j,i,2)=ar(j)-al(j)
            ci(j,i,3)=6.0*s(j,i)-3.0*(al(j)+ar(j))
          else !PCM
            ci(j,i,1)=al(j)
            ci(j,i,2)=0.0
            ci(j,i,3)=0.0
          endif
        enddo !j
      enddo !i
      return
      end subroutine hybgen_ppm_coefs

      subroutine hybgen_ppm_remap(si,pi,dpi,ci,
     &                            so,po,ki,ko,ks,thin)
      implicit none
c
      integer ki,ko,ks
      real    si(ki,ks),pi(ki+1),dpi(ki),ci(ki,ks,3),
     &        so(ko,ks),po(ko+1),thin
c
c-----------------------------------------------------------------------
c  1) remap from one set of vertical cells to another.
c     method: monotonic piecewise parabolic across each input cell
c             the output is the average of the interpolation
c             profile across each output cell.
c     Colella, P. & P.R. Woodward, 1984, J. Comp. Phys., 54, 174-201.
c
c  2) input arguments:
c       si    - initial scalar fields in pi-layer space
c       pi    - initial layer interface depths (non-negative)
c                  pi(   1) is the surface
c                  pi(ki+1) is the bathymetry
c                  pi(k+1) >= pi(k)
c       dpi   - initial layer thicknesses (dpi(k)=pi(k+1)-pi(k))
c       ci    - coefficents from hybgen_ppm_coefs
c                profile(y)=ci.1+ci.2*y+ci.3*y^2, 0<=y<=1
c       ki    - number of  input layers
c       ko    - number of output layers
c       ks    - number of fields
c       po    - target interface depths (non-negative)
c                  po(   1) is the surface
c                  po(ko+1) is the bathymetry (== pi(ki+1))
c                  po(k+1) >= po(k)
c       thin  - layer thickness (>0) that can be ignored
c
c  3) output arguments:
c       so    - scalar fields in po-layer space
c
c  4) Tim Campbell, Mississippi State University, October 2002.
c     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
c-----------------------------------------------------------------------
c
      integer i,k,l,lb,lt
      real    qb0,qb1,qb2,qt0,qt1,qt2,xb,xt,zb,zt,zx,o
      real*8  sz
c
      zx=pi(ki+1) !maximum depth
      zb=max(po(1),pi(1))
      lb=1
      do while (pi(lb+1).lt.zb .and. lb.lt.ki)
        lb=lb+1
      enddo
      do k= 1,ko  !output layers
        zt = zb
        zb = min(po(k+1),zx)
*       write(lp,*) 'k,zt,zb = ',k,zt,zb
        lt=lb !top will always correspond to bottom of previous
              !find input layer containing bottom output interface
        do while (pi(lb+1).lt.zb .and. lb.lt.ki)
          lb=lb+1
        enddo
        if     (zb-zt.le.thin .or. zt.ge.zx) then
          if     (k.ne.1) then
c
c ---       thin or bottomed layer, values taken from layer above
c
            do i= 1,ks
              so(k,i) = so(k-1,i)
            enddo !i
          else !thin surface layer
            do i= 1,ks
              so(k,i) = si(k,i)
            enddo !i
          endif
        else
c
c         form layer averages.
c
*         if     (pi(lb).gt.zt) then
*           write(lp,*) 'bad lb = ',lb
*           stop
*         endif
          xt=(zt-pi(lt))/max(dpi(lt),thin)
          xb=(zb-pi(lb))/max(dpi(lb),thin)
          if     (lt.ne.lb) then  !multiple layers
            qt0 = (1.0-xt)
            qt1 = (1.-xt**2)*0.5
            qt2 = (1.-xt**3)/3.0
            qb0 =     xb
            qb1 =     xb**2 *0.5
            qb2 =     xb**3 /3.0
            do i= 1,ks
              o  = si((lt+lb)/2,i)  !offset to reduce round-off
              sz = dpi(lt)*(  (ci(lt,i,1)-o)*qt0
     &                       +(ci(lt,i,2)+
     &                         ci(lt,i,3) ) *qt1
     &                        -ci(lt,i,3)   *qt2 )
              do l=lt+1,lb-1
                sz = sz+dpi(l)*(si(l,i)-o)
              enddo !l
              sz = sz+dpi(lb)*(  (ci(lb,i,1)-o)*qb0
     &                         +(ci(lb,i,2)+
     &                           ci(lb,i,3) )  *qb1
     &                          -ci(lb,i,3)    *qb2 )
              so(k,i) = o + sz/(zb-zt)  !zb-zt>=thin
            enddo !i
          else  !single layer
            qt0 = (xb-xt)
            qt1 = (xb**2-xt**2)*0.5
            qt2 = (xb**3-xt**3)/3.0
            do i= 1,ks
              sz = dpi(lt)*(  (ci(lt,i,1)-o)*qt0
     &                      +(ci(lt,i,2)+
     &                        ci(lt,i,3) )  *qt1
     &                       -ci(lt,i,3)    *qt2 )
              so(k,i) = o + sz/(zb-zt)  !zb-zt>=thin
            enddo !i
          endif
        endif !thin:std layer
      enddo !k
      return
      end subroutine hybgen_ppm_remap

      subroutine hybgen_weno_coefs(s,dp,lc,ci,kk,ks,thin)
      implicit none
c
      integer kk,ks
      logical lc(kk)
      real    s(kk,ks),dp(kk),ci(kk,ks,2),thin
c
c-----------------------------------------------------------------------
c  1) coefficents for remaping from one set of vertical cells to another.
c     method: monotonic WENO-like alternative to PPM across each input cell
c             a second order polynomial approximation of the profiles
c             using a WENO reconciliation of the slopes to compute the 
c             interfacial values 
c
c     REFERENCE?
c
c  2) input arguments:
c       s     - initial scalar fields in pi-layer space
c       dp    - initial layer thicknesses (>=thin)
c       lc    - use PCM for selected layers
c       kk    - number of layers
c       ks    - number of fields
c       thin  - layer thickness (>0) that can be ignored
c
c  3) output arguments:
c       ci    - coefficents for hybgen_weno_remap
c                ci.1 is value at interface above
c                ci.2 is value at interface below
c                profile(y)=ci.1+ci.2*y+ci.3*y^2, 0<=y<=1
c
c  4) Laurent Debreu, Grenoble.
c     Alan J. Wallcraft,  Naval Research Laboratory,  July 2008.
c-----------------------------------------------------------------------
c
      real, parameter :: dsmll=1.0e-8
c
      integer j,i
      real    q,q01,q02,q001,q002
      real    qdpjm(kk),qdpjmjp(kk),dpjm2jp(kk)
      real    zw(kk+1,3)

      !compute grid metrics
      do j=2,kk-1
        qdpjm(  j) = 1.0/(dp(j-1) +     dp(j))
        qdpjmjp(j) = 1.0/(dp(j-1) +     dp(j) + dp(j+1))
        dpjm2jp(j) =      dp(j-1) + 2.0*dp(j) + dp(j+1)
      enddo !j
      j=kk
        qdpjm(  j) = 1.0/(dp(j-1) +     dp(j))
c
      do i= 1,ks
        do j=2,kk
          zw(j,3) = qdpjm(j)*(s(j,i)-s(j-1,i))
        enddo !j
          j = 1  !PCM first layer
            ci(j,i,1) = s(j,i)
            ci(j,i,2) = s(j,i)
            zw(j,  1) = 0.0
            zw(j,  2) = 0.0
        do j=2,kk-1
          if     (lc(j) .or. dp(j).le.thin) then  !use PCM
            ci(j,i,1) = s(j,i)
            ci(j,i,2) = s(j,i)
            zw(j,  1) = 0.0
            zw(j,  2) = 0.0
          else
            q001 = dp(j)*zw(j+1,3)
            q002 = dp(j)*zw(j,  3)
            if (q001*q002 < 0.0) then
              q001 = 0.0
              q002 = 0.0
            endif
            q01 = dpjm2jp(j)*zw(j+1,3)
            q02 = dpjm2jp(j)*zw(j,  3)
            if     (abs(q001) > abs(q02)) then
              q001 = q02
            endif
            if     (abs(q002) > abs(q01)) then
              q002 = q01
            endif
            q    = (q001-q002)*qdpjmjp(j)
            q001 = q001-q*dp(j+1)
            q002 = q002+q*dp(j-1)

            ci(j,i,2) = s(j,i)+q001
            ci(j,i,1) = s(j,i)-q002
            zw(  j,1) = (2.0*q001-q002)**2
            zw(  j,2) = (2.0*q002-q001)**2
          endif  !PCM:WEND
        enddo !j
          j = kk  !PCM last layer
            ci(j,i,1) = s(j,i)
            ci(j,i,2) = s(j,i)
            zw(j,  1) = 0.0
            zw(j,  2) = 0.0

        do j=2,kk
          q002 = max(zw(j-1,2),dsmll)
          q001 = max(zw(j,  1),dsmll)
          zw(j,3) = (q001*ci(j-1,i,2)+q002*ci(j,i,1))/(q001+q002)
        enddo !j
          zw(   1,3) = 2.0*s( 1,i)-zw( 2,3)  !not used?
          zw(kk+1,3) = 2.0*s(kk,i)-zw(kk,3)  !not used?

        do j=2,kk-1
          if     (.not.(lc(j) .or. dp(j).le.thin)) then  !don't use PCM
            q01  = zw(j+1,3)-s(j,i)
            q02  = s(j,i)-zw(j,3)
            q001 = 2.0*q01
            q002 = 2.0*q02
            if     (q01*q02 < 0.0) then
              q01 = 0.0
              q02 = 0.0
            elseif (abs(q01) > abs(q002)) then
              q01 = q002
            elseif (abs(q02) > abs(q001)) then
              q02 = q001
            endif
            ci(j,i,1) = s(j,i)-q02
            ci(j,i,2) = s(j,i)+q01
          endif  !PCM:WEND
        enddo !j
      enddo !i
      return
      end subroutine hybgen_weno_coefs

      subroutine hybgen_weno_remap(si,pi,dpi,ci,
     &                             so,po,ki,ko,ks,thin)
      implicit none
c
      integer ki,ko,ks
      real    si(ki,ks),pi(ki+1),dpi(ki),ci(ki,ks,2),
     &        so(ko,ks),po(ko+1),thin
c
c-----------------------------------------------------------------------
c  1) remap from one set of vertical cells to another.
c     method: monotonic WENO-like alternative to PPM across each input cell
c             a second order polynomial approximation of the profiles
c             using a WENO reconciliation of the slopes to compute the 
c             interfacial values 
c             the output is the average of the interpolation
c             profile across each output cell.
c
c     REFERENCE?
c
c  2) input arguments:
c       si    - initial scalar fields in pi-layer space
c       pi    - initial layer interface depths (non-negative)
c                  pi(   1) is the surface
c                  pi(ki+1) is the bathymetry
c                  pi(k+1) >= pi(k)
c       dpi   - initial layer thicknesses (dpi(k)=pi(k+1)-pi(k))
c       ci    - coefficents from hybgen_weno_coefs
c                ci.1 is value at interface above
c                ci.2 is value at interface below
c                profile(y)=ci.1+ci.2*y+ci.3*y^2, 0<=y<=1
c       ki    - number of  input layers
c       ko    - number of output layers
c       ks    - number of fields
c       po    - target interface depths (non-negative)
c                  po(   1) is the surface
c                  po(ko+1) is the bathymetry (== pi(ki+1))
c                  po(k+1) >= po(k)
c       thin  - layer thickness (>0) that can be ignored
c
c  3) output arguments:
c       so    - scalar fields in po-layer space
c
c  4) Laurent Debreu, Grenoble.
c     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
c-----------------------------------------------------------------------
c
      integer i,k,l,lb,lt
      real    dpb,dpt,qb0,qb1,qb2,qt0,qt1,qt2,xb,xt,zb,zt,zx,o
      real*8  sz
c
      zx=pi(ki+1) !maximum depth
      zb=max(po(1),pi(1))
      lb=1
      do while (pi(lb+1).lt.zb .and. lb.lt.ki)
        lb=lb+1
      enddo
      do k= 1,ko  !output layers
        zt = zb
        zb = min(po(k+1),zx)
*       write(lp,*) 'k,zt,zb = ',k,zt,zb
        lt=lb !top will always correspond to bottom of previous
              !find input layer containing bottom output interface
        do while (pi(lb+1).lt.zb .and. lb.lt.ki)
          lb=lb+1
        enddo
        if     (zb-zt.le.thin .or. zt.ge.zx) then
          if     (k.ne.1) then
c
c ---       thin or bottomed layer, values taken from layer above
c
            do i= 1,ks
              so(k,i) = so(k-1,i)
            enddo !i
          else !thin surface layer
            do i= 1,ks
              so(k,i) = si(k,i)
            enddo !i
          endif
        else
c
c         form layer averages.
c
*         if     (pi(lb).gt.zt) then
*           write(lp,*) 'bad lb = ',lb
*           stop
*         endif
          xt=(zt-pi(lt))/max(dpi(lt),thin)
          xb=(zb-pi(lb))/max(dpi(lb),thin)
          if     (lt.ne.lb) then  !multiple layers
            dpt = pi(lt+1)-zt
            dpb = zb-pi(lb)
            qt1 = xt*(xt-1.0)
            qt2 = qt1+xt
            qt0 = 1.0-qt1-qt2
            qb1 = (xb-1.0)**2
            qb2 = qb1-1.0+xb
            qb0 = 1.0-qb1-qb2
            do i= 1,ks
              o = si((lt+lb)/2,i)  !offset to reduce round-off
              sz = dpt*(qt0*(si(lt,i)  -o) +
     &                  qt1*(ci(lt,i,1)-o) +
     &                  qt2*(ci(lt,i,2)-o)  )
              do l=lt+1,lb-1
                sz = sz+dpi(l)*(si(l,i) - o)
              enddo !l
              sz  = sz + dpb*(qb0*(si(lb,i)  -o) +
     &                        qb1*(ci(lb,i,1)-o) +
     &                        qb2*(ci(lb,i,2)-o)  )
              so(k,i) = o + sz/(zb-zt)  !zb-zt>=thin
            enddo !i
          else !single layer
            qt1 = xb**2 + xt**2 + xb*xt + 1.0 - 2.0*(xb+xt)
            qt2 = qt1 - 1.0 + (xb+xt)
            qt0 = 1.0 - qt1 - qt2
            do i= 1,ks
              sz=qt0*(si(lt,i)  -o) +
     &           qt1*(ci(lt,i,1)-o) +
     &           qt2*(ci(lt,i,2)-o) 
              so(k,i) = o + sz
            enddo !i
          endif !layers
        endif !thin:std layer
      enddo !k
      return
      end subroutine hybgen_weno_remap

c
c
c> Revision history:
c>
c> Feb. 2000 -- total rewrite to convert to 'newzp' approach
c> Jul. 2000 -- added hybgenj for OpenMP parallelization
c> Oct. 2000 -- added hybgenbj to simplify OpenMP logic
c> Nov. 2000 -- fill massless layers on sea floor with salinity from above
c> Nov. 2000 -- unmixing of deepest inflated layer uses th&T&S from above
c> Nov. 2000 -- ignored isopycnic variance is now 0.002
c> Nov. 2000 -- iterate to correct for cabbeling
c> Nov. 2000 -- allow for "blocking" interior layers
c> Nov. 2000 -- hybflg selects conserved fields (any two of T/S/th)
c> Nov. 2002 -- replace PCM remapping with PLM when non-isopycnal
c> Apr. 2003 -- added dp00i for thinner minimum layers away from the surface
c> Dec. 2003 -- fixed tracer bug when deepest inflated layer is too light
c> Dec. 2003 -- improved water column conservation
c> Dec. 2003 -- compile time option for explicit water column conservation
c> Dec. 2003 -- ignored isopycnic variance is now 0.0001
c> Jan. 2004 -- shifted qqmn,qqmx range now used in cushion function
c> Mar. 2004 -- minimum thickness no longer enforced in near-bottom layers
c> Mar. 2004 -- ignored isopycnic variance is now epsil (i.e. very small)
c> Mar. 2004 -- relaxation to isopycnic layers controled via hybrlx
c> Mar. 2004 -- relaxation removes the need to correct for cabbeling
c> Mar. 2004 -- modified unmixing selection criteria
c> Mar. 2004 -- added isotop (topiso) for isopycnal layer minimum depths
c> Jun. 2005 -- hybrlx (qhybrlx) now input via blkdat.input
c> Jan. 2007 -- hybrlx now only active below "fixed coordinate" surface layers
c> Aug. 2007 -- removed mxlkta logic
c> Sep. 2007 -- added hybmap and hybiso for PCM,PLM,PPM remaper selection
c> Jan. 2008 -- updated logic for two layers (one too dense, other too light)
c> Jul. 2008 -- Added WENO, bugfix to PPM for lcm==.true.
c> Aug. 2008 -- Use WENO-like (vs PPM) for most velocity remapping
c> Aug. 2008 -- Switch more thin near-isopycnal layers to PCM remapping
c> Aug. 2010 -- Added hybthk

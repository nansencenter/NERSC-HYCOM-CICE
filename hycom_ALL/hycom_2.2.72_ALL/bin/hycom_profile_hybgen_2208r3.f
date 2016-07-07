      PROGRAM HYCOM_PROFILE_HYBGEN
      IMPLICIT NONE
C
C  hycom_profile_hybgen - Usage: hycom_profile_hybgen archv.txt blkdat.input archg.txt
C
C                 apply the HYCOM hybrid grid generator to a HYCOM profile.
C                 note that the velocities are unchanged on exit.
C                 based on hybgen version 2.2.08
C
C   archv.txt    is assumed to be an HYCOM archive text profile file
C   blkdat.input is assumed to be an HYCOM parameter file (version 2.1.34)
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
      REAL          RIN(99),PIN(99+1)
      INTEGER       IOS,K,KDM
      LOGICAL       LDEBUG
C
      INTEGER       KDMBLK,NHYBRD,NSIGMA,HYBFLG
      REAL          HYBRLX,QHYBRLX
      REAL          DP00,DP00X,DP00F,DS00,DS00X,DS00F,DP00I,THETA(99)
      REAL          DP0K(99),DP0KP(99),DS0K(99),DSSK(99),DEPTHS,ISOTOP
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
      call blkini(kdmblk,'kdm   ')
      call blkini(nhybrd,'nhybrd')
      call blkini(nsigma,'nsigma')
      call blkinr(dp00,  'dp00  ','(a6," =",f10.4," m")')
      call blkinr(dp00x, 'dp00x ','(a6," =",f10.4," m")')
      call blkinr(dp00f, 'dp00f ','(a6," =",f10.4," ")') 
      call blkinr(ds00,  'ds00  ','(a6," =",f10.4," m")')
      call blkinr(ds00x, 'ds00x ','(a6," =",f10.4," m")')
      call blkinr(ds00f, 'ds00f ','(a6," =",f10.4," ")') 
      call blkinr(dp00i, 'dp00i ','(a6," =",f10.4," m")')
      call blkinr(isotop,'isotop','(a6," =",f10.4," m")')
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
c
      DO K= 1,18
        READ(99,*)
      ENDDO
c
c --- 'hybrlx' = HYBGEN: inverse relaxation coefficient (time steps)
c                 (1.0 for no relaxation)
c --- 'hybflg' = hybrid generator  flag (0=T&S, 1=th&S, 2=th&T)
c
      call blkinr(hybrlx,'hybrlx','(a6," =",f10.4," time steps")')
      call blkini(hybflg,'hybflg')
c
      qhybrlx = 1.0/hybrlx
C
      CLOSE(99)
C
C     INITIALIZE FIXED GRID.
C
      CALL GEOPAR(DP00,DP00X,DP00F,DS00,DS00X,DS00F,DP00I,
     +            NHYBRD,NSIGMA,KDMBLK,
     +            DP0K,DP0KP,DS0K,DSSK)
C
C     COPY PROFILE HEADER TO OUTPUT.
C
      DO K= 1,4
        READ( 11,'(a)') CLINE
        WRITE(21,'(a)') CLINE(1:LEN_TRIM(CLINE))
      ENDDO
      READ( 11,'(a)') CLINE
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
        PIN(K) = P(K)
      ENDDO
      PIN(KDM+1) = P(KDM+1)
C
      DEPTHS = P(KDM+1)*ONEM
      ISOTOP = ISOTOP*ONEM
C
      CALL HYBGEN(T,S,R,DP,THETA,KDM,
     +            NHYBRD,HYBFLG, QHYBRLX,
     +            DP00I,DP0K,DP0KP,DS0K,DSSK,DEPTHS,ISOTOP, LDEBUG)
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
C     OUTPUT THE TRANSFORMED PROFILE.
C
      WRITE(21,'(a,a)')
     &  '#   k',
     &  '    utot    vtot    temp    saln    dens    thkns      dpth'
      DO K= 1,KDM
        THK = P(K+1) - P(K)
        WRITE(21,'(i4,1x,2f8.2,3f8.3,f9.3,f10.3)')
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
     &                  dp0k,dp0kp,ds0k,dssk)
      implicit none
c
      integer nhybrd,nsigma,kdm
      real    dp00,dp00x,dp00f,ds00,ds00x,ds00f,dp00i
      real    dp0k(kdm),dp0kp(kdm),ds0k(kdm),dssk(kdm)
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
c --- logorithmic k-dependence of dp0 (deep z's)
      dp00 =onem*dp00
      dp00x=onem*dp00x
      dp00i=onem*dp00i
      if     (nhybrd.eq.0) then
*       dp0k(1)=thkmin*onem
        dp0k(1)=  20.0*onem
      else
        dp0k(1)=dp00
      endif
      dp0kp(1)=dp0k(1)+onem
      dpm  = dp0k(1)*qonem
      dpms = dpm
      write(lp,*)
      write(lp,135) 1,dp0k(1)*qonem,dpm,dpms
 135  format('dp0k(',i2,') =',f7.2,' m',
     &          '    thkns =',f7.2,' m',
     &          '    depth =',f8.2,' m')
c
      dp0kf=1.0
      do k=2,kdm
        dp0kf=dp0kf*dp00f
        if     (k.le.nhybrd) then
          dp0k(k)=min(dp00*dp0kf,dp00x)
        else
          dp0k(k)=0.0
        endif
        dp0kp(k)=dp0k(k)+onem
        dpm  = dp0k(k)*qonem
        dpms = dpms + dpm
        write(lp,135) k,dp0k(k)*qonem,dpm,dpms
        if     (mnproc.eq.-99) then  ! bugfix that prevents optimization
          write(6,*) 'geopar: dp0kf  = ',dp0kf,    mnproc
          write(6,*) 'geopar: dp0k   = ',dp0k(k),k,mnproc
        endif
      enddo
c
c --- logorithmic k-dependence of ds0 (shallow z-s)
      ds00 =onem*ds00
      ds00x=onem*ds00x
      if     (nhybrd.eq.0) then
*       ds0k(1)=thkmin*onem
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
 130  format('ds0k(',i2,') =',f7.2,' m',
     &          '    thkns =',f7.2,' m',
     &          '    depth =',f8.2,' m')
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
     &                  nhybrd,hybflg, qhybrlx,
     &                  dp00i,dp0k,dp0kp,ds0k,dssk,depths,topiso,
     &                  ldebug)
      implicit none
c
      logical ldebug
      integer kdm, nhybrd,hybflg
      real     temp(1,1,kdm,1),
     &         saln(1,1,kdm,1),
     &         th3d(1,1,kdm,1),
     &           dp(1,1,kdm,1),
     &        theta(1,1,kdm)
      real    qhybrlx
      real    dp00i,dp0k(kdm),dp0kp(kdm),ds0k(kdm),dssk(kdm),
     &        depths(1,1),topiso(1,1)
c
      integer     mnproc,lp
      common/xxx/ mnproc,lp
c
      integer klist(1,1)
      real    p(1,1,kdm+1),dpmixl(1,1,1),q2(1,1,1,1),q2l(1,1,1,1)
      real    tracer(1,1,1,1,1),trcflg(1)
      logical mxlkta,thermo
      integer kk,n, nstep, i0,j0,itest,jtest
      real    epsil,onem,tencm,onemm,qonem,thbase
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
      logical, parameter :: lpcm=.false.      !PCM, instead of PLM, remapping
      logical, parameter :: lconserve=.false. !explicitly conserve each column
c
      double precision tsum,ssum,thsum,trsum(mxtrcr),psum,
     &                 q2sum,q2lsum,rpsum
      double precision asum(  mxtrcr+5,3)
      real             offset(mxtrcr+5)
      real ttem(kdm+1,2),tsal(kdm+1,2),tthe(kdm+1,2),
     &     tq2( kdm+1,2),tq2l(kdm+1,2),ttrc(kdm+1,2,mxtrcr),
     &     dprs(kdm+1),pres(kdm+2),
     &     qhrlx( kdm+1), !relaxation coefficient, from qhybrlx
     &     dp0ij( kdm),   !minimum layer thickness
     &     dp0cum(kdm+1)  !minimum interface depth
      real pwidth,p_hat,p_hat0,p_hat2,p_hat3,hybrlx,
     &     delt,deltm,dels,delsm,thnew,q,qtr,qts,thkbop,
     &     qbot,qcen,qtop,zbot,zbox,zcen,ztop,zthk,dpthin
      integer i,j,k,k1,ka,kbot,kbox,ktop,ksubl,kp,kq,ktr,l,fixlay
ccc   real colint,colins,colout,colous
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
      tencm = onem/10.0
      onemm = onem/1000.0
c
      dpthin = 0.001*onemm
*     thkbop = thkbot*onem
      thkbop =   10.0*onem
      hybrlx = 1.0/qhybrlx
c
******do 1 l=1,isp(j)
c
******do 2 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
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
      enddo
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
c
      if (itest.gt.0 .and. jtest.gt.0) then
        write (lp,'(a/(i6,1x,2f8.3,2f9.3,f9.3))')
     .  'hybgen:   thkns  minthk     dpth  mindpth   hybrlx',
     .  (k,dp(itest,jtest,k,n)*qonem,   dp0ij(k)*qonem,
     .      p(itest,jtest,k+1)*qonem,dp0cum(k+1)*qonem,
     .      1.0/qhrlx(k+1),
     .   k=1,kk)
      endif
c
c --- mark the layer containing the mixed layer base by -ksubl- if the
c --- kraus-turner mixed layer model is selected. otherwise set ksubl to -1.
      if     (mxlkta) then
        ksubl=-1
        do k=1,kk
          if (p(i,j,k  ).lt.dpmixl(i,j,n) .and. 
     &        p(i,j,k+1).ge.dpmixl(i,j,n)+onemm) then
            ksubl=k
          endif
        enddo
        klist(i,j)=ksubl
      else
        ksubl=-1
      endif
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
      k=kp
c
      if     (lunmix .and. !usually .true.
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
              if (i.eq.itest .and. j.eq.jtest) then
                write(lp,'(a,i4,i3,5e12.3)')
     &            'hybgen, 10(+):',
     &            k,ktr,p_hat,p(i,j,k),p(i,j,k-1),
     &            qtr,tracer(i,j,k-1,n,ktr)
                call flush(lp)
              endif !debug
            endif !trcflg
          enddo !ktr
        endif !tracers
        if (mxlmy) then
          qtr=p_hat/(p(i,j,k)-p(i,j,k-1))  !ok because <1.0 and >0.0
          q2( i,j,k-1,n)=q2( i,j,k-1,n)+
     &                     qtr*(q2( i,j,k,n)-q2( i,j,k-1,n))
          q2l(i,j,k-1,n)=q2l(i,j,k-1,n)+
     &                     qtr*(q2l(i,j,k,n)-q2l(i,j,k-1,n))
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
c --- store one-dimensional arrays of -temp-, -saln-, and -p- for the 'old'
c --- vertical grid before attempting to restore isopycnic conditions
      pres(1)=p(i,j,1)
      kp=0
      do k=1,kk
        k1=k+kp
        if (k.ne.ksubl) then
          tthe(k1,1)=th3d(i,j,k,n)
          ttem(k1,1)=temp(i,j,k,n)
          tsal(k1,1)=saln(i,j,k,n)
          do ktr= 1,ntracr
            ttrc(k1,1,ktr)=tracer(i,j,k,n,ktr)
          enddo
          if (mxlmy) then
            tq2( k1,1)=q2( i,j,k,n)
            tq2l(k1,1)=q2l(i,j,k,n)
          endif
          pres(k1+1)=p(i,j,k+1)
          dprs(k1)  =pres(k1+1)-pres(k1)
        else				!  k = ksubl
c ---     expand layer into two sublayers, above and below mixed layer base
          tthe(k,1)=tthe(k-1,1)
          ttem(k,1)=ttem(k-1,1)
          tsal(k,1)=tsal(k-1,1)
          pres(k+1)=dpmixl(i,j,n)
          dprs(k)  =pres(k+1)-pres(k)
          q=(p(i,j,k)-dpmixl(i,j,n))/(p(i,j,k+1)-dpmixl(i,j,n))
          tthe(k+1,1)=th3d(i,j,k,n)+q*(tthe(k,1)-th3d(i,j,k,n))
          ttem(k+1,1)=temp(i,j,k,n)+q*(ttem(k,1)-temp(i,j,k,n))
          tsal(k+1,1)=saln(i,j,k,n)+q*(tsal(k,1)-saln(i,j,k,n))
          pres(k+2)=p(i,j,k+1)
          dprs(k+1)=pres(k+2)-pres(k+1)
          do ktr= 1,ntracr
            ttrc(k  ,1,ktr)=ttrc(k-1,1,ktr)
            ttrc(k+1,1,ktr)=tracer(i,j,k,n,ktr)+
     &                      q*(ttrc(k,1,ktr)-tracer(i,j,k,n,ktr))
          enddo
          kp=1
        endif !k.ne.ksubl:else
      enddo !k
c
      if     (lpcm) then  !usually .false.
c
c ---   PCM (zero slope, recovers original hybgen behaviour).
c
        do k= 1,kk+kp
          tthe(k,2)=0.0
          ttem(k,2)=0.0
          tsal(k,2)=0.0
          do ktr= 1,ntracr
            ttrc(k,2,ktr)=0.0
          enddo
          if (mxlmy) then
            tq2( k,2)=0.0
            tq2l(k,2)=0.0
          endif
        enddo
      else
c
c ---   PLM (non-zero slope, but no new extrema)
c ---   layer value is (:,1)-0.5*(:,2) at top    interface,
c ---              and (:,1)+0.5*(:,2) at bottom interface.
c ---   still use PCM for isopycnal layers, because we don't
c ---   want detrainment to change the density of these layers.
c
c ---   monotonized central-difference limiter (van Leer, 1977,
c ---   JCP 23 pp 276-299).  For a discussion of PLM limiters, see
c ---   Finite Volume Methods for Hyperbolic Problems by R.J. Leveque.
c
        do k= 1,kk+kp
          if     (k.eq.1 .or. k.eq.kk+kp .or.
     &            dprs(k).le.dpthin      .or.
     &            dprs(k).gt.dp0kp(k)        ) then
c ---       top, bottom, thin and isopycnal layers have zero slope.
            tthe(k,2)=0.0
            ttem(k,2)=0.0
            tsal(k,2)=0.0
            do ktr= 1,ntracr
              ttrc(k,2,ktr)=0.0
            enddo
            if (mxlmy) then
              tq2( k,2)=0.0
              tq2l(k,2)=0.0
            endif
          else
c ---       interior non-isopycnal layer
c ---       use qcen in place of 0.5 to allow for non-uniform grid
            qcen = dprs(k)/(dprs(k)+0.5*(dprs(k-1)+dprs(k+1)))
c
            ztop = 2.0*(tthe(k,  1)-tthe(k-1,1))
            zbot = 2.0*(tthe(k+1,1)-tthe(k,  1))
            zcen =qcen*(tthe(k+1,1)-tthe(k-1,1))
            if     (ztop*zbot.gt.0.0) then
              tthe(k,2)=sign(min(abs(zcen),abs(zbot),abs(ztop)),zbot)
            else
              tthe(k,2)=0.0  !local extrema, so no slope
            endif
c
            ztop = 2.0*(ttem(k,  1)-ttem(k-1,1))
            zbot = 2.0*(ttem(k+1,1)-ttem(k,  1))
            zcen =qcen*(ttem(k+1,1)-ttem(k-1,1))
            if     (ztop*zbot.gt.0.0) then
              ttem(k,2)=sign(min(abs(zcen),abs(zbot),abs(ztop)),zbot)
            else
              ttem(k,2)=0.0  !local extrema, so no slope
            endif
c
            ztop = 2.0*(tsal(k,  1)-tsal(k-1,1))
            zbot = 2.0*(tsal(k+1,1)-tsal(k,  1))
            zcen =qcen*(tsal(k+1,1)-tsal(k-1,1))
            if     (ztop*zbot.gt.0.0) then
              tsal(k,2)=sign(min(abs(zcen),abs(zbot),abs(ztop)),zbot)
            else
              tsal(k,2)=0.0  !local extrema, so no slope
            endif
c
            do ktr= 1,ntracr
              ztop = 2.0*(ttrc(k,  1,ktr)-ttrc(k-1,1,ktr))
              zbot = 2.0*(ttrc(k+1,1,ktr)-ttrc(k,  1,ktr))
              zcen =qcen*(ttrc(k+1,1,ktr)-ttrc(k-1,1,ktr))
              if     (ztop*zbot.gt.0.0) then
                ttrc(k,2,ktr)=sign(min(abs(zcen),abs(zbot),abs(ztop)),
     &                             zbot)
              else
                ttrc(k,2,ktr)=0.0  !local extrema, so no slope
              endif
            enddo
c
            if (mxlmy) then
              ztop = 2.0*(tq2( k,  1)-tq2( k-1,1))
              zbot = 2.0*(tq2( k+1,1)-tq2( k,  1))
              zcen =qcen*(tq2( k+1,1)-tq2( k-1,1))
              if     (ztop*zbot.gt.0.0) then
                tq2( k,2)=sign(min(abs(zcen),abs(zbot),abs(ztop)),zbot)
              else
                tq2( k,2)=0.0  !local extrema, so no slope
              endif
              ztop = 2.0*(tq2l(k,  1)-tq2l(k-1,1))
              zbot = 2.0*(tq2l(k+1,1)-tq2l(k,  1))
              zcen =qcen*(tq2l(k+1,1)-tq2l(k-1,1))
              if     (ztop*zbot.gt.0.0) then
                tq2l(k,2)=sign(min(abs(zcen),abs(zbot),abs(ztop)),zbot)
              else
                tq2l(k,2)=0.0  !local extrema, so no slope
              endif
            endif !mxlmy
          endif !top/bottom/thin/isopycnal or PLM layer
        enddo !k
      endif !PCM:PLM
c
ccc   colint=temp(i,j,1,n)*(p(i,j,2)-p(i,j,1))
ccc   colins=saln(i,j,1,n)*(p(i,j,2)-p(i,j,1))
c
c --- try to restore isopycnic conditions by moving layer interfaces
c
      do 88 k=2,nhybrd
            if (i.eq.itest .and. j.eq.jtest) then
              write(cinfo,'(a9,i2.2,1x)') '  do 88 k=',k
 109          format (i9,2i5,a,a/(i9,8x,a,a,i3,f9.2,f8.2,f9.2,f8.2))
              write(lp,109) nstep,itest+i0,jtest+j0,
     .         cinfo,':    othkns  odpth    nthkns  ndpth',
     .        (nstep,cinfo,':',k1,
     .         (pres(k1+1)-
     .          pres(k1)   )*qonem,
     .          pres(k1+1)  *qonem,
     .         (p(itest,jtest,k1+1)-
     .          p(itest,jtest,k1)   )*qonem,
     .          p(itest,jtest,k1+1)  *qonem,k1=1,kk)
              call flush(lp)
            endif !debug
c
ccc   colint=colint+temp(i,j,k,n)*(p(i,j,k+1)-p(i,j,k))
ccc   colins=colins+saln(i,j,k,n)*(p(i,j,k+1)-p(i,j,k))
c
c --- maintain constant thickness in layers 1 to fixlay
      if (k.le.fixlay+1) then
        p_hat=-999.0*dp0k(k-1)
        go to 9
      end if
c
c --- is density noticeably different from isopycnic reference value?
      if (abs(th3d(i,j,k,n)-theta(i,j,k)).lt.epsil) go to 8
c
      if (th3d(i,j,k,n).le.theta(i,j,k)) go to 7	!  layer too light
c
c --- water in layer k is too dense. try to dilute with water from layer k-1
c
c --- if layer k-1 is too light, thicken the thinner of the two,
c --- i.e. skip this layer if it is thicker.
      if (th3d(i,j,k-1,n).lt.theta(i,j,k-1) .and.
     &    p(i,j,k+1)-p(i,j,k).gt.p(i,j,k)-p(i,j,k-1)) go to 8
c
      if (i.eq.itest .and. j.eq.jtest) then
        write(lp,'(a,f8.4)') 'hybgen, too dense:',
     &                       th3d(i,j,k,n)-theta(i,j,k)
        call flush(lp)
      endif !debug
c
      if     ((theta(i,j,k)-th3d(i,j,k-1,n)).le.epsil) then
        p_hat=min(p(i,j,k-1),
     &            p(i,j,k)-999.0*dp0k(k-1))  ! take entire layer k-1
      else
        q=(theta(i,j,k)-th3d(i,j,k,n))/(theta(i,j,k)-th3d(i,j,k-1,n))  ! -ve
        p_hat=p(i,j,k)+q*(p(i,j,k+1)-p(i,j,k))
      end if
c
c --- maintain minimum layer thickess of layer k-1.
c
 9    continue
      p_hat=min(p_hat,2.0*p(i,j,max(k-1,3))) !layer no thicker than layers above
      p_hat0=p_hat
      p_hat=max(p(i,j,k-1)+dp0ij(k-1),
     &      min(p(i,j,k+1)-dp0ij(k),
     &          p(i,j,k-1)+cushn(p_hat0-p(i,j,k-1),dp0ij(k-1)) ))
      if (i.eq.itest .and. j.eq.jtest) then
        write(lp,'(a,3f8.2)') 'hybgen, 9: ',
     &   (p_hat0-p(i,j,k-1))*qonem,
     &   cushn(p_hat0-p(i,j,k-1),dp0ij(k-1))*qonem,p_hat*qonem
        call flush(lp)
      endif !debug
c
c --- if isopycnic conditions cannot be achieved because of a blocking
c --- layer in the interior ocean, move interface k-1 (and k-2 if
c --- necessary) upward
c
      if     (k.le.2) then
c ---   do nothing.
      else if (p_hat.ge.p(i,j,k) .and. 
     &         p(i,j,k-1).gt.dp0cum(k-1)+onem .and.
     &         (p(i,j,kk+1)-p(i,j,k-1).lt.thkbop .or.
     &          p(i,j,k-1) -p(i,j,k-2).gt.qqmx*dp0ij(k-2))) then  ! k.gt.2
        p_hat2=p(i,j,k-2)+cushn(p(i,j,k-1)-p_hat+p_hat0-p(i,j,k-2),
     &                          dp0ij(k-2))
        if (p_hat2.lt.p(i,j,k-1)-onemm) then
          p(i,j,k-1)=(1.0-qhrlx(k-1))*p(i,j,k-1) + 
     &                    qhrlx(k-1) *max(p_hat2,2.0*p(i,j,k-1)-p_hat)
          if (i.eq.itest .and. j.eq.jtest) then
            write(lp,'(a,i3.2,f8.2)') 'hybgen,  1blocking :',
     .                              k-1,p(i,j,k-1)*qonem
            call flush(lp)
          endif !debug
          p_hat=p(i,j,k-1)+cushn(p_hat0-p(i,j,k-1),dp0ij(k-1))
        else if (k.le.3) then
c ---     do nothing.
        else if (p(i,j,k-2).gt.dp0cum(k-2)+onem .and.
     &           (p(i,j,kk+1)-p(i,j,k-2).lt.thkbop .or.
     &            p(i,j,k-2) -p(i,j,k-3).gt.qqmx*dp0ij(k-3))) then  ! k.gt.3
          p_hat3=p(i,j,k-3)+cushn(p(i,j,k-2)-p_hat+p_hat0-p(i,j,k-3),
     &                            dp0ij(k-3))
          if (p_hat3.lt.p(i,j,k-2)-onemm) then
            p(i,j,k-2)=(1.0-qhrlx(k-2))*p(i,j,k-2) +
     &                 qhrlx(k-2)*max(p_hat3,2.0*p(i,j,k-2)-p(i,j,k-1))
            if (i.eq.itest .and. j.eq.jtest) then
              write(lp,'(a,i3.2,f8.2)') 'hybgen,  2blocking :',
     .                                k-2,p(i,j,k-2)*qonem
              call flush(lp)
            endif !debug
            p_hat2=p(i,j,k-2)+cushn(p(i,j,k-1)-p_hat+p_hat0-p(i,j,k-2),
     &                              dp0ij(k-2))
            if (p_hat2.lt.p(i,j,k-1)-onemm) then
              p(i,j,k-1)=(1.0-qhrlx(k-1))*p(i,j,k-1) +
     &                    qhrlx(k-1) *max(p_hat2,2.0*p(i,j,k-1)-p_hat)
              if (i.eq.itest .and. j.eq.jtest) then
                write(lp,'(a,i3.2,f8.2)') 'hybgen,  3blocking :',
     .                                  k-1,p(i,j,k-1)*qonem
                call flush(lp)
              endif !debug
              p_hat=p(i,j,k-1)+cushn(p_hat0-p(i,j,k-1),dp0ij(k-1))
            end if
          end if
        end if
      end if
c
      if (p_hat.le.p(i,j,k)) then
c
c --- upper intfc moves up. entrain layer k-1 water into layer k
c
        q = p_hat - p(i,j,k)
        p(i,j,k)=max(pres(k-1),min(pres(k+1),p(i,j,k)+qhrlx(k)*q))
c ---   qhrlx is relaxation coefficient (inverse baroclinic time steps)
        if (i.eq.itest .and. j.eq.jtest) then
          write(lp,'(a,i3.2,2f8.2)') 'hybgen, entrain(k) :',
     .                               k,p_hat*qonem,p(i,j,k)*qonem
          call flush(lp)
        endif !debug
c
        if     (k.ge.3) then
          do ka= k,kk
c ---       layer ka no thicker than sum of layers above
            if     (p(i,j,ka+1).le.2.0*p(i,j,ka)) then
              exit
            endif
            if (i.eq.itest .and. j.eq.jtest) then
              write(lp,'(a,i3.2,2f8.2)') 'hybgen, nothick(ka) :',
     .                                   ka,p(i,j,ka+1)*qonem,
     .                                   2.0*p(i,j,ka)*qonem
              call flush(lp)
            endif !debug
            q = 2.0*p(i,j,ka) - p(i,j,ka+1)
            p(i,j,ka+1)=p(i,j,ka+1)+qhrlx(ka+1)*q
          enddo !ka
        endif !k>=3
c
      else				!  p_hat > p(i,j,k)
c
c --- move upper interface down and entrain layer k water into layer k-1
c --- if maintenance of minimum thickness forces the interface below the
c --- bottom, then all interfaces below k must reside at the bottom
c
        p_hat0 = p_hat
        q      = p_hat - p(i,j,k)
        p_hat  = max(pres(k-1),min(pres(k+1),p(i,j,k)+qhrlx(k)*q))
c ---   qhrlx is relaxation coefficient (inverse baroclinic time steps)
        if (p_hat.gt.p(i,j,kk+1)) then
          do ka=k,kk
            p(i,j,k)=p(i,j,kk+1)
          enddo
        else
          p_hat=min(p_hat,p(i,j,k+1))
          p(i,j,k)=p_hat
        endif
        if (i.eq.itest .and. j.eq.jtest) then
          write(lp,'(a,i3.2,2f8.2)') 'hybgen, entrain(k-):',
     .                            k,p_hat0*qonem,p(i,j,k)*qonem
          call flush(lp)
        endif !debug
c
c --- do we need to inflate layer k by lowering  l o w e r  interface?
c
        if (k.lt.kk) then
          if (p(i,j,k+2).lt.p(i,j,kk+1)) then
            p_hat=p(i,j,k+1)
            go to 6
          endif
        endif
      endif
      go to 8
c
c --- water in layer k is too light. try to dilute with water from layer k+1
c
 7    continue
c
c --- is layer k touching the sea floor?
      if (p(i,j,k+1).eq.p(i,j,kk+1)) go to 8
c
c --- are we below any KT mixed layer?
      if (mxlkta .and. k.le.klist(i,j)) go to 8
c
c --- if layer k+1 is too dense, thicken the thinner of the two,
c --- i.e. skip this layer if it is thicker.
      if (th3d(i,j,k+1,n).gt.theta(i,j,k+1) .and.
     &    p(i,j,k+1).le.2.0*p(i,j,max(k,3)) .and.
     &    p(i,j,k+1)-p(i,j,k).gt.p(i,j,k+2)-p(i,j,k+1)) go to 8
c
      if (i.eq.itest .and. j.eq.jtest) then
        write(lp,'(a,f8.4)') 'hybgen, too light:',
     &                       theta(i,j,k)-th3d(i,j,k,n)
        call flush(lp)
      endif !debug
c
      if     ((th3d(i,j,k+1,n)-theta(i,j,k)).le.epsil) then
        p_hat=max(p(i,j,k+2),
     &            p(i,j,k+1)+999.0*dp0k(k))  ! take entire layer k+1
      else
        q=(th3d(i,j,k,n)-theta(i,j,k))/(th3d(i,j,k+1,n)-theta(i,j,k))  ! -ve
        p_hat=p(i,j,k+1)+q*(p(i,j,k)-p(i,j,k+1))
      endif
c
c --- if layer k+1 does not touch the bottom and does not contain the KT
c --- mixed layer base, then maintain minimum thicknesses of layers
c --- k and k+1 to the greatest extent possible, but permit layers to
c --- collapse to zero thickness at the bottom
 6    continue
      if (mxlkta .and. k.eq.klist(i,j)) go to 8
      if     (p(i,j,k+2).lt.p(i,j,kk+1)) then
        if     (p(i,j,kk+1)-p(i,j,k).gt.dp0ij(k)+dp0ij(k+1)) then
          p_hat=p(i,j,k+2)-cushn(p(i,j,k+2)-p_hat,dp0ij(k+1))
        endif
        p_hat=p(i,j,k)  +  max(p_hat    -p(i,j,k),dp0ij(k))
        p_hat=min(p_hat,
     &            2.0*p(i,j,max(k,3)), !layer no thicker than layers above
     &            max(0.5*(p(i,j,k+1)+p(i,j,k+2)),
     &                     p(i,j,k+2)-dp0ij(k+1)))
      else
        p_hat=min(p(i,j,k+2),p_hat)
      endif !p.k+2<p.kk+1
c
      if (p_hat.gt.p(i,j,k+1)) then
c
c ---   entrain layer k+1 water into layer k. 
        p(i,j,k+1)=(1.0-qhrlx(k+1))*p(i,j,k+1) + qhrlx(k+1)*p_hat
      if (i.eq.itest .and. j.eq.jtest) then
        write(lp,'(a,i3.2,f8.2)') 'hybgen, entrain(k+):',
     .                            k+1,p(i,j,k+1)*qonem
        call flush(lp)
      endif !debug
      elseif (p_hat.eq.2.0*p(i,j,max(k,3))) then
c
c ---   entrain layer k water into layer k+1. 
c ---   note that p_hat > p(i,j,k)
        p(i,j,k+1)=(1.0-qhrlx(k+1))*p(i,j,k+1) + qhrlx(k+1)*p_hat
      if (i.eq.itest .and. j.eq.jtest) then
        write(lp,'(a,i3.2,f8.2)') 'hybgen, nothick(k+):',
     .                            k+1,p(i,j,k+1)*qonem
        call flush(lp)
      endif !debug
      endif !entrain:no thicker
 8    continue
c
c --- if layer above is too thin, move interface down.
c --- dp0ij(k-1) is "fixed coordinate" layer thickness.
      p_hat0=p(i,j,k-1)+dp0ij(k-1)
      if (p_hat0.gt.p(i,j,k)   .and.
     &    p_hat0.lt.p(i,j,kk+1)     ) then
        p(i,j,k)=p_hat0
      endif !p_hat0
c
c --- enforce interface order (is this necessary?), usually inexpensive.
      do ka= k+1,kk
        if     (p(i,j,ka).ge.p(i,j,k)) then
          exit  ! usually get here quickly
        else
          p(i,j,ka) = p(i,j,k)
        endif
      enddo !ka
 88   continue !k
c
ccc   write (lp,'(a,1p,2e16.8,e9.1)') 'temp. column integral:',
ccc  .  colint,colout,(colout-colint)/colint
ccc   write (lp,'(a,1p,2e16.8,e9.1)') 'saln. column integral:',
ccc  .  colins,colous,(colous-colins)/colins
c
c --- remap scalar field profiles from the 'old' vertical
c --- grid onto the 'new' vertical grid, using PLM
c
      if     (lconserve) then  !usually .false.
        do ktr=1,ntracr+5
          asum(ktr,1) = 0.d0
          asum(ktr,2) = 0.d0
          asum(ktr,3) = 0.d0
        enddo !ktr
      endif
c
      zbot=0.0
      kbot=1
      do k=1,kk
        ztop=zbot  !top is bottom of previous layer
        ktop=kbot
        if     (ztop.ge.pres(ktop+1)) then
          ktop=ktop+1
        endif
            if (ktop.eq.kbot .and. i.eq.itest .and. j.eq.jtest) then
              write(lp,'(a,2i3,2f23.16)')
     &          'k,ktop,dp =',k,ktop,
     &          ztop*qonem,(pres(ktop+1)-ztop)*qonem
            endif !debug
        zbot=p(i,j,k+1)
        zthk=zbot-ztop
        dp(i,j,k,n)=zthk
            if (i.eq.itest .and. j.eq.jtest) then
              write(lp,'(a,2i3,3f23.16)')
     &          'k,z =',k,k,ztop*qonem,zbot*qonem,zthk*qonem
            endif !debug
        if     (zthk.gt.dpthin .and. ztop.lt.p(i,j,kk+1)) then
c ---     normal layer
          kbot=ktop
          do while (pres(kbot+1).lt.zbot.and.kbot.lt.kk+kp)
            kbot=kbot+1
          enddo
              if (i.eq.itest .and. j.eq.jtest) then
                write(lp,'(a,2i3,2f23.16)')
     &            'k,p =',ktop,kbot,
     &            pres(ktop)*qonem,pres(kbot+1)*qonem
              endif !debug
c
c ---     include thin adjacent layers in sum
          zbox=zbot
          do k1= k+1,kk
            if     (p(i,j,k1+1)-p(i,j,k1).gt.dpthin) then
              exit !thick layer
            else
              zbox=p(i,j,k1+1)  !include thin adjacent layers
              if     (zbox.eq.p(i,j,kk+1)) then
                exit !at bottom
              endif
            endif
          enddo
          zthk=zbox-ztop
c
          kbox=ktop
          do while (pres(kbox+1).lt.zbox.and.kbox.lt.kk+kp)
            kbox=kbox+1
          enddo
              if (i.eq.itest .and. j.eq.jtest) then
                write(lp,'(a,2i3,2f23.16)')
     &            'k,px=',ktop,kbox,
     &            pres(ktop)*qonem,pres(kbox+1)*qonem
              endif !debug
          if     (ktop.eq.kbox) then
c ---       single layer
            if     (p(i,j,k)  .ne.pres(kbox)   .or.
     &              p(i,j,k+1).ne.pres(kbox+1)     ) then
c ---         part of a single layer
              zcen = 0.5*(ztop+zbox)
              if     (dprs(kbox).gt.dpthin) then
                q = 0.5 - (pres(kbox+1)-zcen)/dprs(kbox)
              else
                q = 0.0
              endif
              if     (hybflg.eq.0) then  !T&S
                temp(i,j,k,n)=ttem(kbox,1)+q*ttem(kbox,2)
                saln(i,j,k,n)=tsal(kbox,1)+q*tsal(kbox,2)
                th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
              elseif (hybflg.eq.1) then  !th&S
                th3d(i,j,k,n)=tthe(kbox,1)+q*tthe(kbox,2)
                saln(i,j,k,n)=tsal(kbox,1)+q*tsal(kbox,2)
                temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,
     &                               saln(i,j,k,n))
              elseif (hybflg.eq.2) then  !th&T
                th3d(i,j,k,n)=tthe(kbox,1)+q*tthe(kbox,2)
                temp(i,j,k,n)=ttem(kbox,1)+q*ttem(kbox,2)
                saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,
     &                               temp(i,j,k,n))
              endif
              do ktr= 1,ntracr
                tracer(i,j,k,n,ktr)=ttrc(kbox,1,ktr)+q*ttrc(kbox,2,ktr)
              enddo
              if (mxlmy) then
                q2( i,j,k,n)=tq2( kbox,1)+q*tq2( kbox,2)
                q2l(i,j,k,n)=tq2l(kbox,1)+q*tq2l(kbox,2)
              endif
            else
c ---         all of a single layer
              temp(i,j,k,n)=ttem(kbox,1)
              saln(i,j,k,n)=tsal(kbox,1)
              th3d(i,j,k,n)=tthe(kbox,1)
              do ktr= 1,ntracr
                tracer(i,j,k,n,ktr)=ttrc(kbox,1,ktr)
              enddo
              if (mxlmy) then
                q2( i,j,k,n)=tq2( kbox,1)
                q2l(i,j,k,n)=tq2l(kbox,1)
              endif
              q = 0.0 !for debugging only
            endif !part:all of single layer
              if (i.eq.itest .and. j.eq.jtest) then
                write(lp,'(a,i3,3f23.16)')
     &            'k,q =',kbox,q,ttrc(kbox,1,1),tracer(i,j,k,n,1)
*    &            'k,q =',kbox,q,tthe(kbox,1),th3d(i,j,k,n)
              endif !debug
          else
c ---       multiple layers.
            if     (ktop.le.k .and. kbox.ge.k) then
              ka = k
            elseif (kbox-ktop.ge.3) then
              ka = (kbox+ktop)/2
            elseif (dprs(ktop).ge.dprs(kbox)) then
              ka = ktop
            else
              ka = kbox
            endif !choose ka
c ---       calculate as perturbation from layer ka (reduces roundoff)
            offset(1)=ttem(ka,1)
            offset(2)=tsal(ka,1)
            offset(3)=tthe(ka,1)
            if (mxlmy) then
              offset(4)=tq2( ka,1)
              offset(5)=tq2l(ka,1)
            endif
            do ktr= 1,ntracr
              offset(ktr+5)=ttrc(ka,1,ktr)
            enddo
c
            qtop = pres(ktop+1)-ztop !partial layer thickness
            zcen = 0.5*(ztop+pres(ktop+1))
            if     (dprs(ktop).gt.dpthin) then
              q = 0.5 - (pres(ktop+1)-zcen)/dprs(ktop)
            else
              q = 0.0
            endif
            tsum =((ttem(ktop,1)+q*ttem(ktop,2))-offset(1))*qtop
            ssum =((tsal(ktop,1)+q*tsal(ktop,2))-offset(2))*qtop
            thsum=((tthe(ktop,1)+q*tthe(ktop,2))-offset(3))*qtop
            do ktr= 1,ntracr
              trsum(ktr)=((  ttrc(ktop,1,ktr)+
     &                     q*ttrc(ktop,2,ktr) )-offset(ktr+5))*qtop
            enddo
            if (mxlmy) then
              q2sum =((tq2( ktop,1)+q*tq2( ktop,2))-offset(4))*qtop
              q2lsum=((tq2l(ktop,1)+q*tq2l(ktop,2))-offset(5))*qtop
            endif
              if (i.eq.itest .and. j.eq.jtest) then
                write(lp,'(a,i3,3f23.16)')
     &            'k,f =',ktop,qtop/max(dprs(ktop),dpthin),
     &                         qtop/zthk,trsum(1)/zthk
*    &            'k,f =',ktop,qtop/zthk,ttrc(ktop,1,1),trsum(1)/zthk
*    &            'k,f =',ktop,qtop/zthk,tthe(ktop,1),thsum/zthk
              endif !debug
c
            do k1= ktop+1,kbox-1
              tsum =tsum +(ttem(k1,1)-offset(1))*dprs(k1)
              ssum =ssum +(tsal(k1,1)-offset(2))*dprs(k1)
              thsum=thsum+(tthe(k1,1)-offset(3))*dprs(k1)
              do ktr= 1,ntracr
                trsum(ktr)=trsum(ktr)+
     &                       (ttrc(k1,1,ktr)-offset(ktr+5))*dprs(k1)
              enddo
              if (mxlmy) then
                q2sum =q2sum +(tq2( k1,1)-offset(4))*dprs(k1)
                q2lsum=q2lsum+(tq2l(k1,1)-offset(5))*dprs(k1)
              endif
              if (i.eq.itest .and. j.eq.jtest) then
                write(lp,'(a,i3,3f23.16)')
     &            'k,f =',k1,1.0,dprs(k1)/zthk,trsum(1)/zthk
*    &            'k,f =',k1,dprs(k1)/zthk,ttrc(k1,1,1),trsum(1)/zthk
*    &            'k,f =',k1,dprs(k1)/zthk,tthe(k1,1),thsum/zthk
              endif !debug
            enddo !k1
c
            qbot = zbox-pres(kbox) !partial layer thickness
            zcen = 0.5*(pres(kbox)+zbox)
            if     (dprs(kbox).gt.dpthin) then
              q = 0.5 - (pres(kbox+1)-zcen)/dprs(kbox)
            else
              q = 0.0
            endif
            tsum =tsum +((ttem(kbox,1)+q*ttem(kbox,2))-offset(1))*qbot
            ssum =ssum +((tsal(kbox,1)+q*tsal(kbox,2))-offset(2))*qbot
            thsum=thsum+((tthe(kbox,1)+q*tthe(kbox,2))-offset(3))*qbot
            do ktr= 1,ntracr
              trsum(ktr)=trsum(ktr)+
     &                     ((  ttrc(kbox,1,ktr)+
     &                       q*ttrc(kbox,2,ktr) )-offset(ktr+5))*qbot
            enddo
            if (mxlmy) then
              q2sum =q2sum +
     &                 ((tq2( kbox,1)+q*tq2( kbox,2))-offset(4))*qbot
              q2lsum=q2lsum+
     &                 ((tq2l(kbox,1)+q*tq2l(kbox,2))-offset(5))*qbot
            endif
              if (i.eq.itest .and. j.eq.jtest) then
                write(lp,'(a,i3,3f23.16)')
     &            'k,f =',kbox,qbot/max(dprs(kbox),dpthin),
     &                         qbot/zthk,trsum(1)/zthk
*    &            'k,f =',kbox,qbot/zthk,ttrc(kbox,1,1),trsum(1)/zthk
*    &            'k,f =',kbox,qbot/zthk,tthe(kbox,1),thsum/zthk
              endif !debug
c
            rpsum=1.0d0/zthk
            if     (hybflg.eq.0) then  !T&S
              temp(i,j,k,n)=offset(1)+tsum*rpsum
              saln(i,j,k,n)=offset(2)+ssum*rpsum
              th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
            elseif (hybflg.eq.1) then  !th&S
              th3d(i,j,k,n)=offset(3)+thsum*rpsum
              saln(i,j,k,n)=offset(2)+ ssum*rpsum
              temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
            elseif (hybflg.eq.2) then  !th&T
              th3d(i,j,k,n)=offset(3)+thsum*rpsum
              temp(i,j,k,n)=offset(1)+ tsum*rpsum
              saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
            endif
            do ktr= 1,ntracr
              tracer(i,j,k,n,ktr)=offset(ktr+5)+trsum(ktr)*rpsum
              if     (trcflg(ktr).ne.2) then !not a temperature tracer
c ---           correct for round-off below zero
                tracer(i,j,k,n,ktr)=max(tracer(i,j,k,n,ktr),0.0)
              endif
            enddo
            if (mxlmy) then
              q2( i,j,k,n)=max(dsmll,offset(4)+q2sum *rpsum)
              q2l(i,j,k,n)=max(dsmll,offset(5)+q2lsum*rpsum)
            endif
          endif !single or multiple layers
        else
c ---     thin or bottomed layer
            if (i.eq.itest .and. j.eq.jtest) then
              write(lp,'(a,i3)')
     &          'thin k =',k
            endif !debug
          if (k.le.nhybrd) then
c ---       fill with fluid from above
            th3d(i,j,k,n)=th3d(i,j,k-1,n)
            saln(i,j,k,n)=saln(i,j,k-1,n)
            temp(i,j,k,n)=temp(i,j,k-1,n)
          elseif (hybflg.ne.2) then
c ---       fill with saln from above
            th3d(i,j,k,n)=theta(i,j,k)
            saln(i,j,k,n)=saln(i,j,k-1,n)
            temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
            saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
          else
c ---       fill with temp from above
            th3d(i,j,k,n)=theta(i,j,k)
            temp(i,j,k,n)=temp(i,j,k-1,n)
            saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
          endif
          do ktr= 1,ntracr
            tracer(i,j,k,n,ktr)=tracer(i,j,k-1,n,ktr)
          enddo
          if (mxlmy) then
            q2( i,j,k,n)=q2( i,j,k-1,n)
            q2l(i,j,k,n)=q2l(i,j,k-1,n)
          endif
        endif !normal:thin layer
c
        if     (lconserve) then  !usually .false.
          asum(1,1) = asum(1,1) + ttem(    k,1)*dprs(k)
          asum(1,2) = asum(1,2) + temp(i,j,k,n)*zthk
          asum(2,1) = asum(2,1) + tsal(    k,1)*dprs(k)
          asum(2,2) = asum(2,2) + saln(i,j,k,n)*zthk
          asum(3,1) = asum(3,1) + tthe(    k,1)*dprs(k)
          asum(3,2) = asum(3,2) + th3d(i,j,k,n)*zthk
          if (mxlmy) then
            asum(4,1) = asum(4,1) + tq2(     k,1)*dprs(k)
            asum(4,2) = asum(4,2) +  q2( i,j,k,n)*zthk
            asum(5,1) = asum(5,1) + tq2l(    k,1)*dprs(k)
            asum(5,2) = asum(5,2) +  q2l(i,j,k,n)*zthk
          endif
          do ktr= 1,ntracr
            asum(ktr+5,1) = asum(ktr+5,1) + ttrc(      k,1,ktr)*dprs(k)
            asum(ktr+5,2) = asum(ktr+5,2) + tracer(i,j,k,n,ktr)*zthk
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
        do ktr=1,ntracr+5
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
            saln(i,j,k,n)=saln(i,j,k,n)*(1.0+offset(2))
            th3d(i,j,k,n)=th3d(i,j,k,n)*(1.0+offset(3))
            temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,
     &                           saln(i,j,k,n))
          elseif (hybflg.eq.2) then  !th&T
            temp(i,j,k,n)=temp(i,j,k,n)*(1.0+offset(1))
            th3d(i,j,k,n)=th3d(i,j,k,n)*(1.0+offset(3))
            saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,
     &                           temp(i,j,k,n))
          endif
          if (mxlmy) then
            q2( i,j,k,n)=q2( i,j,k,n)*(1.0+offset(4))
            q2l(i,j,k,n)=q2l(i,j,k,n)*(1.0+offset(5))
          endif
          do ktr= 1,ntracr
            tracer(i,j,k,n,ktr)=tracer(i,j,k,n,ktr)*(1.0+offset(ktr+5))
          enddo !ktr
c
          if     (.false.) then !debugging
            zthk = dp(i,j,k,n)
            asum(1,3) = asum(1,3) + temp(i,j,k,n)*zthk
            asum(2,3) = asum(2,3) + saln(i,j,k,n)*zthk
            asum(3,3) = asum(3,3) + th3d(i,j,k,n)*zthk
            if (mxlmy) then
              asum(4,3) = asum(4,3) +  q2( i,j,k,n)*zthk
              asum(5,3) = asum(5,3) +  q2l(i,j,k,n)*zthk
            endif
            do ktr= 1,ntracr
              asum(ktr+5,3) = asum(ktr+5,3) + tracer(i,j,k,n,ktr)*zthk
            enddo !ktr
          endif !debuging
        enddo !k
c
        if     (.false. .and. !debugging
     &          i.eq.itest .and. j.eq.jtest) then
          do ktr= 1,ntracr+5
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
          ktr=6
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
cdiag if (i.eq.itest .and. j.eq.jtest) then
cdiag   write (lp,103) nstep,itest+i0,jtest+j0,
cdiag&  '    hybgen, do 22:  temp    saln    dens     thkns     dpth',
cdiag&  (k,ttem(k,1),tsal(k,1),tthe(k,1)+thbase,
cdiag&   (pres(k+1)-pres(k))*qonem,pres(k+1)*qonem,
cdiag&   k,temp(i,j,k,n),saln(i,j,k,n),
cdiag&   th3d(i,j,k,n)+thbase,dp(i,j,k,n)*qonem,
cdiag&   p(i,j,k+1)*qonem,
cdiag&  k=1,kk),
cdiag&  (k,ttem(k,1),tsal(k,1),tthe(k,1)+thbase,
cdiag&   (pres(k+1)-pres(k))*qonem,pres(k+1)*qonem,
cdiag&   k=kk+1,kp)
cdiag   call flush(lp)
cdiag endif !debug
c
 2    continue  !i
c
c --- to avoid roundoff errors in -dpudpv- after restart, make sure p=p(dp)
c
      do 1 k=1,kk
*     do 1 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
 1    continue !i;k;l
c
*     if(mxlkta) then
*       do 71 l=1,isp(j)
*       do 71 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
c
*       do 71 k=1,kk
*       if(dpmixl(i,j,n).gt.p(i,j,k  ) .and.
*    &     dpmixl(i,j,n).le.p(i,j,k+1)      ) then
*         t1sav(i,j,n)=tmix(i,j)
*         s1sav(i,j,n)=smix(i,j)
*         tmlb(i,j,n)=temp(i,j,k,n)
*         smlb(i,j,n)=saln(i,j,k,n)
*         nmlb(i,j,n)=k
*       end if
c
*71     continue !k;i;l
*     end if
c
      return
      end
c
c
c> hybgen Revision history:
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

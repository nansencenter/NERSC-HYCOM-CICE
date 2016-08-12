      subroutine hybgen(temp,saln,th3d,dp,theta,kdm,
     &                  nhybrd,hybflg,
     &                  dp00i,dp0k,dp0kp,ds0k,dssk,depths, ldebug)
      implicit none
c
      logical ldebug
      integer kdm, nhybrd,hybflg
      real     temp(1,1,kdm,1),
     &         saln(1,1,kdm,1),
     &         th3d(1,1,kdm,1),
     &           dp(1,1,kdm,1),
     &        theta(1,1,kdm)
      real    dp00i,dp0k(kdm),dp0kp(kdm),ds0k(kdm),dssk(kdm),depths(1,1)
c
      integer     mnproc,lp
      common/xxx/ mnproc,lp
c
c --- ---------------------
c --- hybrid grid generator
c --- ---------------------
c
      integer i,j,kk,n, lp,nstep, i0,j0,itest,jtest
      real    epsil,onem,onemm,qonem,thbase
      real    p(1,1,kdm+1)
c
      double precision tsum,ssum,thsum,
     &                 q2sum,q2lsum,rpsum
      real ttem(kdm+1,2),tsal(kdm+1,2),tthe(kdm+1,2),
     &     dprs(kdm+1),pres(kdm+2),dp0ij(kdm),dp0cum(kdm+1)
      real pwidth,p_hat,p_hat0,p_hat2,p_hat3,
     &     delt,deltm,dels,delsm,thnew,q,qtr,qts,thkbop,
     &     qbot,qcen,qtop,zbot,zcen,ztop,zthk,dpthin
      integer k,k1,ka,kbot,ktop,ksubl,kp,ktr,l,iter
      character*12 cinfo
c
      double precision dsmll
      parameter       (dsmll=1.d-8)
c
c --- c u s h i o n   function (from Bleck & Benjamin, 1992):
c --- if delp >= qqmx*dp0 >>  dp0, -cushn- returns -delp-
c --- if delp <= qqmn*dp0 << -dp0, -cushn- returns  -dp0-
c
      real       qqmn,qqmx,cusha,cushb
      parameter (qqmn=-2.0, qqmx=4.0)  ! traditional range
*     parameter (qqmn=-4.0, qqmx=6.0)  ! somewhat wider range
      parameter (cusha=qqmn**2*(qqmx-1.0)/(qqmx-qqmn)**2)
      parameter (cushb=1.0/qqmn)
c
      real qq,cushn,delp,dp0
*     include 'stmt_fns.h'
c-----------------------------------------------------------------------------
      real sig,dsigdt,dsigds,tofsig,sofsig
c
      real    a0,a1,a2,cubr,cubq,cuban,cubrl,cubim
      real    r,s,t,prs,ylat
c
      real       ahalf,athird,afourth
      parameter (ahalf  =1./2.)
      parameter (athird =1./3.)
      parameter (afourth=1./4.)
c
c --- coefficients for sigma-0 (based on Brydon & Sun fit)
      real       c1,c2,c3,c4,c5,c6,c7
      parameter (c1=-1.36471E-01, c2= 4.68181E-02, c3= 8.07004E-01,
     &           c4=-7.45353E-03, c5=-2.94418E-03,
     &           c6= 3.43570E-05, c7= 3.48658E-05)
c
c --- coefficients for sigma-2 (based on Brydon & Sun fit)
csig2 real       c1,c2,c3,c4,c5,c6,c7
csig2 parameter (c1= 9.77093E+00, c2=-2.26493E-02, c3= 7.89879E-01,
csig2&           c4=-6.43205E-03, c5=-2.62983E-03,
csig2&           c6= 2.75835E-05, c7= 3.15235E-05)
c
c --- auxiliary statements for finding root of 3rd degree polynomial
      a0(s)=(c1+c3*s)/c6
      a1(s)=(c2+c5*s)/c6
      a2(s)=(c4+c7*s)/c6
      cubq(s)=athird*a1(s)-(athird*a2(s))**2
      cubr(r,s)=athird*(0.5*a1(s)*a2(s)-1.5*(a0(s)-r/c6))
     &           -(athird*a2(s))**3
c --- if q**3+r**2>0, water is too dense to yield real root at given
c --- salinitiy. setting q**3+r**2=0 in that case is equivalent to
c --- lowering sigma until a double real root is obtained.
      cuban(r,s)=athird*atan2(sqrt(max(0.,-(cubq(s)**3+cubr(r,s)**2))),
     &                        cubr(r,s))
      cubrl(r,s)=sqrt(-cubq(s))*cos(cuban(r,s))
      cubim(r,s)=sqrt(-cubq(s))*sin(cuban(r,s))
c
c --- -----------------
c --- equation of state
c --- -----------------
c
c --- sigma-theta as a function of temp (deg c) and salinity (mil)
c --- (friedrich-levitus 3rd degree polynomial fit)
c
      sig(t,s)=(c1+c3*s+t*(c2+c5*s+t*(c4+c7*s+c6*t)))
c
c --- d(sig)/dt
      dsigdt(t,s)=(c2+c5*s+2.*t*(c4+c7*s+1.5*c6*t))
c
c --- d(sig)/ds
      dsigds(t,s)=(c3+t*(c5+t*c7))
c
c --- temp (deg c) as a function of sigma and salinity (mil)
      tofsig(r,s)=-cubrl(r,s)+sqrt(3.)*cubim(r,s)-athird*a2(s)
c
c --- salinity (mil) as a function of sigma and temperature (deg c)
      sofsig(r,t)=(r-c1-t*(c2+t*(c4+c6*t)))/(c3+t*(c5+c7*t))
c-----------------------------------------------------------------------------
      qq(   delp,dp0)=max(qqmn, min(qqmx, delp/dp0))
      cushn(delp,dp0)=dp0*
     &                (1.0+cusha*(1.0-cushb*qq(delp,dp0))**2)*
     &                max(1.0, delp/(dp0*qqmx))
c
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
      onemm = onem/1000.0
c
      dpthin = 0.1*onemm
*     thkbop = thkbot*onem
      thkbop =   10.0*onem
c
*     do 1 l=1,isp(j)
c
*     do 2 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
c
      dp0cum(1)=0.0
      dp0ij( 1)=min( dp0k(1), max( ds0k(1), dssk(1)*depths(i,j) ) )
      dp0cum(2)=dp0cum(1)+dp0ij(1)
      p(i,j, 2)=p(i,j,1)+dp(i,j,1,n)
      do k=2,kk
        if     (dp0k(k).le.dp00i) then
          q =      dp0k(k)
        else
          q = max( dp00i,
     &             dp0k(k) * dp0k(k)/
     &                       max( dp0k(k),
     &                            p(i,j,k)-dp0cum(k) ) )
        endif
        dp0ij( k)  =min( q,max( ds0k(k), dssk(k)*depths(i,j) ) )
        dp0cum(k+1)=dp0cum(k)+dp0ij(k)
        p(i,j, k+1)=p(i,j,k)+dp(i,j,k,n)
      enddo
c
c --- mark the layer containing the mixed layer base by -ksubl- if the
c --- kraus-turner mixed layer model is selected. otherwise set ksubl to -1.
*     if     (mxlkta .and. thermo) then
*       ksubl=-1
*       do k=1,kk
*         if (p(i,j,k  ).lt.dpmixl(i,j,n) .and. 
*    &        p(i,j,k+1).ge.dpmixl(i,j,n)+onemm) then
*           ksubl=k
*         endif
*       enddo
*       klist(i,j)=ksubl
*     else
        ksubl=-1
*     endif
c
c --- does layer touch sea floor?
      do 10 k=3,kk
      if (p(i,j,k+1).gt.p(i,j,kk+1)-dpthin) then
        if (dp(i,j,k,n).le.dpthin) then
          if (k.le.nhybrd) then
c ---       fill massless layers on sea floor with fluid from above
            th3d(i,j,k,n)=th3d(i,j,k-1,n)
            saln(i,j,k,n)=saln(i,j,k-1,n)
            temp(i,j,k,n)=temp(i,j,k-1,n)
          elseif (hybflg.ne.2) then
c ---       fill massless layers on sea floor with saln from above
            th3d(i,j,k,n)=max(theta(i,j,k), th3d(i,j,k-1,n))
            saln(i,j,k,n)=saln(i,j,k-1,n)
            temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
            saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
          else
c ---       fill massless layers on sea floor with temp from above
            th3d(i,j,k,n)=max(theta(i,j,k), th3d(i,j,k-1,n))
            temp(i,j,k,n)=temp(i,j,k-1,n)
            saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
          endif
*       else if (k.eq.-99) then  ! always .false. at run time
        else if (k.le.nhybrd                              .and.
     &           (theta(i,j,k) -th3d(i,j,k,  n)).gt.0.002 .and.
     &           (th3d(i,j,k,n)-th3d(i,j,k-1,n)).gt.0.002      ) then
c
c ---     water in deepest inflated layer is too light.
c ---     split layer into 2 sublayers, one near the desired density
c ---     and one exactly matching the T&S properties of layer k-1.
c ---     To prevent "runaway" T or S, the result satisfies either
c ---       abs(T.k - T.k-1) <= abs(T.k-2 - T.k-1) or
c ---       abs(S.k - S.k-1) <= abs(S.k-2 - S.k-1) 
c ---     It is also limited to a 50% change in layer thickness.
c
          delsm=abs(saln(i,j,k-2,n)-saln(i,j,k-1,n))
          dels =abs(saln(i,j,k-1,n)-saln(i,j,k,  n))
          deltm=abs(temp(i,j,k-2,n)-temp(i,j,k-1,n))
          delt =abs(temp(i,j,k-1,n)-temp(i,j,k,  n))
c ---     sanity check on deltm and delsm
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
     &      (theta(i,j,k)-th3d(i,j,k-1,n))
          q=min(q,qts/(1.0+qts))  ! upper sublayer <= 50% of total
          p_hat=q*(p(i,j,k+1)-p(i,j,k))
          p(i,j,k)=p(i,j,k)+p_hat
c
          if     (hybflg.eq.0) then  !T&S
            temp(i,j,k,n)=temp(i,j,k,n)+(q/(1.0-q))*(temp(i,j,k,n)  -
     &                                               temp(i,j,k-1,n) )
            saln(i,j,k,n)=saln(i,j,k,n)+(q/(1.0-q))*(saln(i,j,k,n)  -
     &                                               saln(i,j,k-1,n) )
            th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
          elseif (hybflg.eq.1) then  !th&S
            th3d(i,j,k,n)=th3d(i,j,k,n)+(q/(1.0-q))*(th3d(i,j,k,n)  -
     &                                               th3d(i,j,k-1,n) )
            saln(i,j,k,n)=saln(i,j,k,n)+(q/(1.0-q))*(saln(i,j,k,n)  -
     &                                               saln(i,j,k-1,n) )
            temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
          elseif (hybflg.eq.2) then  !th&T
            th3d(i,j,k,n)=th3d(i,j,k,n)+(q/(1.0-q))*(th3d(i,j,k,n)  -
     &                                               th3d(i,j,k-1,n) )
            temp(i,j,k,n)=temp(i,j,k,n)+(q/(1.0-q))*(temp(i,j,k,n)  -
     &                                               temp(i,j,k-1,n) )
            saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
          endif
          if (i.eq.itest .and. j.eq.jtest) then
            write(lp,'(a,i3,f6.3,5f8.3)')
     &        'hybgen, 10(+):',
     &        k,q,temp(i,j,k,n),saln(i,j,k,n),
     &            th3d(i,j,k,n)+thbase,theta(i,j,k)+thbase
            call flush(lp)
          endif
        else
          if (i.eq.itest .and. j.eq.jtest) then
            write(lp,'(a,i3,f6.3,5f8.3)')
     &        'hybgen, 10(-):',
     &        k,0.0,temp(i,j,k,n),saln(i,j,k,n),
     &            th3d(i,j,k,n)+thbase,theta(i,j,k)+thbase
            call flush(lp)
          endif
        endif
      endif
 10   continue
c
c --- store one-dimensional arrays of -temp-, -saln-, and -p- for the 'old'
c --- vertical grid before attempting to restore isopycnic conditions
      pres(1)=p(i,j,1)
      kp=0
      do k=1,kk
        k1=k+kp
*       if (k.ne.ksubl) then
          tthe(k1,1)=th3d(i,j,k,n)
          ttem(k1,1)=temp(i,j,k,n)
          tsal(k1,1)=saln(i,j,k,n)
          pres(k1+1)=p(i,j,k+1)
          dprs(k1)  =pres(k1+1)-pres(k1)
*       else				!  k = ksubl
c ---     expand layer into two sublayers, above and below mixed layer base
*         tthe(k,1)=tthe(k-1,1)
*         ttem(k,1)=ttem(k-1,1)
*         tsal(k,1)=tsal(k-1,1)
*         pres(k+1)=dpmixl(i,j,n)
*         dprs(k)  =pres(k+1)-pres(k)
*         q=(p(i,j,k)-dpmixl(i,j,n))/(p(i,j,k+1)-dpmixl(i,j,n))
*         tthe(k+1,1)=th3d(i,j,k,n)+q*(tthe(k,1)-th3d(i,j,k,n))
*         ttem(k+1,1)=temp(i,j,k,n)+q*(ttem(k,1)-temp(i,j,k,n))
*         tsal(k+1,1)=saln(i,j,k,n)+q*(tsal(k,1)-saln(i,j,k,n))
*         pres(k+2)=p(i,j,k+1)
*         dprs(k+1)=pres(k+2)-pres(k+1)
*         kp=1
*       endif !k.ne.ksubl:else
      enddo !k
c
      if     (.false.) then
c
c ---   PCM (zero slope, recovers original hybgen behaviour).
c
        do k= 1,kk+kp
          tthe(k,2)=0.0
          ttem(k,2)=0.0
          tsal(k,2)=0.0
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
 109          format (19x,a,a/(i9,8x,a,a,i3,f9.2,f8.2,f9.2,f8.2))
              write(lp,109) 
     .         cinfo,':    othkns  odpth    nthkns  ndpth',
     .        (nstep,cinfo,':',k1,
     .         (pres(k1+1)-
     .          pres(k1)   )*qonem,
     .          pres(k1+1)  *qonem,
     .         (p(itest,jtest,k1+1)-
     .          p(itest,jtest,k1)   )*qonem,
     .          p(itest,jtest,k1+1)  *qonem,k1=1,kk)
              call flush(lp)
            endif
c
ccc   colint=colint+temp(i,j,k,n)*(p(i,j,k+1)-p(i,j,k))
ccc   colins=colins+saln(i,j,k,n)*(p(i,j,k+1)-p(i,j,k))
c
c --- maintain constant thickness in layer 1
      if (k.eq.2) then
        p_hat=-999.0*dp0k(1)
        go to 9
      end if
c
c --- are we dealing with a massless layer on the sea floor?
      if (p(i,j,k).gt.p(i,j,kk+1)-dpthin) then
        if (k.le.nhybrd) then
c ---     fill massless layers on sea floor with fluid from above
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
      endif
c
c --- is density noticeably different from isopycnic reference value?
      if (abs(th3d(i,j,k,n)-theta(i,j,k)).lt.0.002) go to 8
c
      if (th3d(i,j,k,n).le.theta(i,j,k)) go to 7	!  layer too light
c
c --- water in layer k is too dense. try to dilute with water from layer k-1
c
c --- if layer k-1 is too light, perform this adjustment every other pair
c --- of time steps beneath the z-coordinate domain. otherwise, layers
c --- k-1 and k may not be able to reach their target densities
      if (mod(nstep-1,4).le.1 .and.
     &    p(i,j,k).gt.dp0cum(k)+onem .and.
     &    th3d(i,j,k-1,n).lt.theta(i,j,k-1)) go to 8
c
      if (i.eq.itest .and. j.eq.jtest) then
        write(lp,'(a,f8.4)') 'hybgen, too dense:',
     &                       th3d(i,j,k,n)-theta(i,j,k)
        call flush(lp)
      endif
c
      if     ((theta(i,j,k)-th3d(i,j,k-1,n)).le.epsil) then
        p_hat=min(p(i,j,k-1),
     &            p(i,j,k)-999.0*dp0k(k-1))  ! take entire layer k-1
      else
        q=(theta(i,j,k)-th3d(i,j,k,n))/(theta(i,j,k)-th3d(i,j,k-1,n))  ! -ve
        p_hat=p(i,j,k)+q*(p(i,j,k+1)-p(i,j,k))
c
c ---   correct for cabbeling by performing an iterative procedure to make
c ---   sure that the target density is achieved to a tolerance of 0.0001.
c ---   if the water becomes more dense, do not move the interface
c
        if     (hybflg.eq.0 .and.
     &          p_hat.gt.p(i,j,k-1) .and. 
     &          p(i,j,k).gt.dp0cum(k)+onem) then
          do 58 iter=1,5
            q=(p(i,j,k)-p_hat)/max(p(i,j,k+1)-p_hat,onemm)
            delt=q*(temp(i,j,k-1,n)-temp(i,j,k,n))
            dels=q*(saln(i,j,k-1,n)-saln(i,j,k,n))
            thnew=sig(temp(i,j,k,n)+delt,saln(i,j,k,n)+dels)-thbase
            if (thnew.gt.th3d(i,j,k,n)) then  ! more dense
              p_hat=p(i,j,k)
              go to 158
            endif
            if (abs(thnew-theta(i,j,k)).le.0.0001) then
              go to 158
            endif
            q=(theta(i,j,k)-thnew)/(theta(i,j,k)-th3d(i,j,k-1,n))
            p_hat=min(p(i,j,k),p_hat+q*(p(i,j,k+1)-p_hat))
  58      continue
 158      continue
        end if
      end if
c
c --- maintain minimum layer thickess of layer k-1.
c
 9    continue
      p_hat0=p_hat
      p_hat=max(p(i,j,k-1)+dp0ij(k-1),
     &      min(p(i,j,k+1)-dp0ij(k),
     &          p(i,j,k-1)+cushn(p_hat0-p(i,j,k-1),dp0ij(k-1)) ))
      if (i.eq.itest .and. j.eq.jtest) then
        write(lp,'(a,2f8.2)') 'hybgen, 9: ',
     &   (p_hat0-p(i,j,k-1))*qonem,
     &   cushn(p_hat0-p(i,j,k-1),dp0ij(k-1))*qonem
        call flush(lp)
      endif
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
          p(i,j,k-1)=max(p_hat2,2.0*p(i,j,k-1)-p_hat)
          if (i.eq.itest .and. j.eq.jtest) then
            write(lp,'(a,i3.2,f8.2)') 'hybgen,  1blocking :',
     .                              k-1,p(i,j,k-1)*qonem
            call flush(lp)
          endif
          p_hat=p(i,j,k-1)+cushn(p_hat0-p(i,j,k-1),dp0ij(k-1))
        else if (k.le.3) then
c ---     do nothing.
        else if (p(i,j,k-2).gt.dp0cum(k-2)+onem .and.
     &           (p(i,j,kk+1)-p(i,j,k-2).lt.thkbop .or.
     &            p(i,j,k-2) -p(i,j,k-3).gt.qqmx*dp0ij(k-3))) then  ! k.gt.3
          p_hat3=p(i,j,k-3)+cushn(p(i,j,k-2)-p_hat+p_hat0-p(i,j,k-3),
     &                            dp0ij(k-3))
          if (p_hat3.lt.p(i,j,k-2)-onemm) then
            p(i,j,k-2)=max(p_hat3,2.0*p(i,j,k-2)-p(i,j,k-1))
            if (i.eq.itest .and. j.eq.jtest) then
              write(lp,'(a,i3.2,f8.2)') 'hybgen,  2blocking :',
     .                                k-2,p(i,j,k-2)*qonem
              call flush(lp)
            endif
            p_hat2=p(i,j,k-2)+cushn(p(i,j,k-1)-p_hat+p_hat0-p(i,j,k-2),
     &                              dp0ij(k-2))
            if (p_hat2.lt.p(i,j,k-1)-onemm) then
              p(i,j,k-1)=max(p_hat2,2.0*p(i,j,k-1)-p_hat)
              if (i.eq.itest .and. j.eq.jtest) then
                write(lp,'(a,i3.2,f8.2)') 'hybgen,  3blocking :',
     .                                  k-1,p(i,j,k-1)*qonem
                call flush(lp)
              endif
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
        p(i,j,k)=p_hat
        if (i.eq.itest .and. j.eq.jtest) then
          write(lp,'(a,i3.2,f8.2)') 'hybgen, entrain(k) :',
     .                            k,p(i,j,k)*qonem
          call flush(lp)
        endif
c
      else				!  p_hat > p(i,j,k)
c
c --- move upper interface down and entrain layer k water into layer k-1
c
        p_hat=min(p_hat,p(i,j,k+1))
        p(i,j,k)=p_hat
        if (i.eq.itest .and. j.eq.jtest) then
          write(lp,'(a,i3.2,f8.2)') 'hybgen, entrain(k-):',
     .                            k,p(i,j,k)*qonem
          call flush(lp)
        endif
c
c --- do we need to inflate layer k by lowering  l o w e r  interface?
c
        if (k.lt.kk) then
          p_hat=p(i,j,k+1)
          go to 6
        end if
      end if
      go to 8
c
c --- water in layer k is too light. try to dilute with water from layer k+1
c
 7    continue
c
c --- is this layer touching the sea floor?
      if (p(i,j,k+1).gt.p(i,j,kk+1)-dpthin .and. k.gt.2) go to 8
c
c --- are we below any KT mixed layer?
*     if (mxlkta .and. k.le.klist(i,j)) go to 8
c
c --- if layer k+1 is too dense, perform this adjustment every other pair
c --- of time steps beneath the z-coordinate domain. otherwise, layers
c --- k and k+1 may not be able to reach their target densities
      if (mod(nstep-1,4).ge.2 .and.
     &    p(i,j,k+1).gt.dp0cum(k+1)+onem   .and.
     &    th3d(i,j,k+1,n).gt.theta(i,j,k+1)     ) go to 8
c
      if (i.eq.itest .and. j.eq.jtest) then
        write(lp,'(a,f8.4)') 'hybgen, too light:',
     &                       theta(i,j,k)-th3d(i,j,k,n)
        call flush(lp)
      endif
c
      if     ((th3d(i,j,k+1,n)-theta(i,j,k)).le.epsil) then
        p_hat=max(p(i,j,k+2),
     &            p(i,j,k+1)+999.0*dp0k(k))  ! take entire layer k+1
      else
        q=(th3d(i,j,k,n)-theta(i,j,k))/(th3d(i,j,k+1,n)-theta(i,j,k))  ! +ve
        p_hat=p(i,j,k+1)+q*(p(i,j,k)-p(i,j,k+1))
c
c ---   correct for cabbeling by performing an iterative procedure to make
c ---   sure that the target density is achieved to a tolerance of 0.0001.
c ---   if the water becomes less dense, do not move the interface
c
        if     (hybflg.eq.0 .and.
     &          p_hat.lt.p(i,j,k+2) .and. 
     &          p(i,j,k+1).gt.dp0cum(k+1)) then
          do 59 iter=1,5
            q=(p_hat-p(i,j,k+1))/max(p_hat-p(i,j,k),onemm)
            delt=q*(temp(i,j,k+1,n)-temp(i,j,k,n))
            dels=q*(saln(i,j,k+1,n)-saln(i,j,k,n))
            thnew=sig(temp(i,j,k,n)+delt,saln(i,j,k,n)+dels)-thbase
            if (thnew.lt.th3d(i,j,k,n)) then  ! less dense
              p_hat=p(i,j,k+1)
              go to 159
            endif
            if (abs(thnew-theta(i,j,k)).le.0.0001) then
              go to 159
            endif
            q=(thnew-theta(i,j,k))/(th3d(i,j,k+1,n)-theta(i,j,k))
            p_hat=max(p(i,j,k+1),p_hat+q*(p(i,j,k)-p_hat))
  59      continue
 159      continue
        endif
      endif
c
c --- maintain minimum thickess of layers k and k+1
 6    continue
*     if (mxlkta .and. k.eq.klist(i,j)) go to 8
      if     (k.lt.nhybrd) then
        p_hat=p(i,j,k+2)-cushn(p(i,j,k+2)-p_hat,dp0ij(k+1))
      endif
      p_hat=p(i,j,k)  +  max(p_hat    -p(i,j,k),dp0ij(k))
      p_hat=min(p_hat,
     &          max(0.5*(p(i,j,k+1)+p(i,j,k+2)),
     &              p(i,j,k+2)-dp0ij(k+1)))
c
      if (p_hat.gt.p(i,j,k+1)) then
c
c --- entrain layer k+1 water into layer k
c
        p(i,j,k+1)=p_hat
        if (i.eq.itest .and. j.eq.jtest) then
          write(lp,'(a,i3.2,f8.2)') 'hybgen, entrain(k+):',
     .                            k+1,p(i,j,k+1)*qonem
          call flush(lp)
        endif
      end if
 8    continue
c
c --- if layer above is too thin, move upper interface down
      p_hat0=p(i,j,k-1)+dp0ij(k-1)
      if (p_hat0.gt.p(i,j,k)   .and.
     &    p_hat0.lt.p(i,j,kk+1)     ) then
        p(i,j,k)=p_hat0
      endif
c
c --- inforce interface order (is this necessary?), usually inexpensive.
      do ka= k+1,kk
        if     (p(i,j,ka).ge.p(i,j,k)) then
          exit  ! usually get here quickly
        else
          p(i,j,ka) = p(i,j,k)
        endif
      enddo
 88   continue
c
ccc   write (lp,'(a,1p,2e16.8,e9.1)') 'temp. column integral:',
ccc  .  colint,colout,(colout-colint)/colint
ccc   write (lp,'(a,1p,2e16.8,e9.1)') 'saln. column integral:',
ccc  .  colins,colous,(colous-colins)/colins
c
c --- remap scalar field profiles from the 'old' vertical
c --- grid onto the 'new' vertical grid, using PLM
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
         write(lp,'(a,2i4,f9.3,1pe15.6)')
     &     'k,ktop,dp = ',k,ktop,
     &     ztop*qonem,(pres(ktop+1)-ztop)*qonem
       endif
        zbot=p(i,j,k+1)
        zthk=zbot-ztop
        dp(i,j,k,n)=zthk
       if (i.eq.itest .and. j.eq.jtest) then
         write(lp,'(a,2i4,3f9.3)')
     &     'k,z = ',k,k,ztop*qonem,zbot*qonem,zthk*qonem
       endif
        if     (zthk.gt.dpthin .and. ztop.lt.p(i,j,kk+1)) then
c ---     normal layer
          kbot=ktop
          do while (pres(kbot+1).lt.zbot.and.kbot.lt.kk+kp)
            kbot=kbot+1
          enddo
         if (i.eq.itest .and. j.eq.jtest) then
           write(lp,'(a,2i4,2f9.3)')
     &            'k,p = ',ktop,kbot,
     &            pres(ktop)*qonem,pres(kbot+1)*qonem
         endif
          if     (ktop.eq.kbot) then
c ---       single layer
            if     (p(i,j,k)  .ne.pres(kbot)   .or.
     &              p(i,j,k+1).ne.pres(kbot+1)     ) then
c ---         part of a single layer
              zcen = 0.5*(ztop+zbot)
              if     (dprs(kbot).gt.dpthin) then
                q = 0.5 - (pres(kbot+1)-zcen)/dprs(kbot)
              else
                q = 0.0
              endif
              if     (hybflg.eq.0) then  !T&S
                temp(i,j,k,n)=ttem(kbot,1)+q*ttem(kbot,2)
                saln(i,j,k,n)=tsal(kbot,1)+q*tsal(kbot,2)
                th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
              else
                th3d(i,j,k,n)=tthe(kbot,1)+q*tthe(kbot,2)
                if     (abs(th3d(i,j,k,n)-theta(i,j,k)).gt.0.001) then  !T&S
                  temp(i,j,k,n)=ttem(kbot,1)+q*ttem(kbot,2)
                  saln(i,j,k,n)=tsal(kbot,1)+q*tsal(kbot,2)
                  th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
                elseif (hybflg.eq.1) then  !th&S, layer is approx. isopycnal
                  saln(i,j,k,n)=tsal(kbot,1)+q*tsal(kbot,2)
                  temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,
     &                                 saln(i,j,k,n))
                elseif (hybflg.eq.2) then  !th&T, layer is approx. isopycnal
                  temp(i,j,k,n)=ttem(kbot,1)+q*ttem(kbot,2)
                  saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,
     &                                 temp(i,j,k,n))
                endif
              endif
            endif !part of single layer
          else
c ---       multiple layers.
            qtop = pres(ktop+1)-ztop !partial layer thickness
            zcen = 0.5*(ztop+pres(ktop+1))
            if     (dprs(ktop).gt.dpthin) then
              q = 0.5 - (pres(ktop+1)-zcen)/dprs(ktop)
            else
              q = 0.0
            endif
            thsum=(tthe(ktop,1)+q*tthe(ktop,2))*qtop
            tsum =(ttem(ktop,1)+q*ttem(ktop,2))*qtop
            ssum =(tsal(ktop,1)+q*tsal(ktop,2))*qtop
c
            do k1= ktop+1,kbot-1
              thsum=thsum+tthe(k1,1)*dprs(k1)
              tsum =tsum +ttem(k1,1)*dprs(k1)
              ssum =ssum +tsal(k1,1)*dprs(k1)
            enddo !k1
c
            qbot = zbot-pres(kbot) !partial layer thickness
            zcen = 0.5*(pres(kbot)+zbot)
            if     (dprs(kbot).gt.dpthin) then
              q = 0.5 - (pres(kbot+1)-zcen)/dprs(kbot)
            else
              q = 0.0
            endif
            thsum=thsum+(tthe(kbot,1)+q*tthe(kbot,2))*qbot
            tsum =tsum +(ttem(kbot,1)+q*ttem(kbot,2))*qbot
            ssum =ssum +(tsal(kbot,1)+q*tsal(kbot,2))*qbot
c
            rpsum=1.0d0/zthk
            if     (hybflg.eq.0) then  !T&S
              temp(i,j,k,n)=tsum*rpsum
              saln(i,j,k,n)=ssum*rpsum
              th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
            else
              th3d(i,j,k,n)=thsum*rpsum
*             if     (abs(th3d(i,j,k,n)-theta(i,j,k)).gt.0.001) then  !T&S
              if     (abs(th3d(i,j,k,n)-theta(i,j,k)).gt.9.999) then  !T&S
                temp(i,j,k,n)=tsum*rpsum
                saln(i,j,k,n)=ssum*rpsum
                th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
              elseif (hybflg.eq.1) then  !th&S, layer is approx. isopycnal
                saln(i,j,k,n)=ssum*rpsum
                temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
              elseif (hybflg.eq.2) then  !th&T, layer is approx. isopycnal
                temp(i,j,k,n)=tsum*rpsum
                saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
              endif
            endif
          endif !single or multiple layers
        else
c ---     thin or bottomed layer
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
        endif !normal:thin layer
      enddo !k
c
 103  format (i9,2i5,a/(33x,i3,2f8.3,f8.3,f9.3,f9.2))
      if (i.eq.itest .and. j.eq.jtest) then
        write (lp,103) nstep,itest+i0,jtest+j0,
     &  '    hybgen, do 22:  temp    saln    dens     thkns     dpth',
     &  (k,ttem(k,1),tsal(k,1),tthe(k,1)+thbase,
     &   (pres(k+1)-pres(k))*qonem,pres(k+1)*qonem,
     &   k,temp(i,j,k,n),saln(i,j,k,n),
     &   th3d(i,j,k,n)+thbase,dp(i,j,k,n)*qonem,
     &   p(i,j,k+1)*qonem,
     &  k=1,kk),
     &  (k,ttem(k,1),tsal(k,1),tthe(k,1)+thbase,
     &   (pres(k+1)-pres(k))*qonem,pres(k+1)*qonem,
     &   k=kk+1,kp)
        call flush(lp)
      endif
c
 2    continue  !i
c
c --- to avoid roundoff errors in -dpudpv- after restart, make sure p=p(dp)
c
      do 1 k=1,kk
*     do 1 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
 1    continue
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
*71     continue
*     end if
c
      return
      end

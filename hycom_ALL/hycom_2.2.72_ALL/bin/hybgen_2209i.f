      subroutine hybgen(m,n)
      use mod_xc    ! HYCOM communication interface
      use mod_pipe  ! HYCOM debugging interface
c
c --- hycom version 1.0
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n
c
c --- ---------------------
c --- hybrid grid generator
c --- ---------------------
c
      logical, parameter :: lpipe_hybgen=.false.  !for debugging
c
      integer   i,j,k,l
      character text*12
c
 103  format (i9,2i5,a/(33x,i3,2f8.3,f8.3,f9.3,f9.2))
cdiag if (itest.gt.0 .and. jtest.gt.0) then
cdiag   write (lp,103) nstep,itest+i0,jtest+j0,
cdiag.  '  entering hybgen:  temp    saln    dens     thkns     dpth',
cdiag.  (k,temp(itest,jtest,k,n),saln(itest,jtest,k,n),
cdiag.  th3d(itest,jtest,k,n)+thbase,dp(itest,jtest,k,n)*qonem,
cdiag.  p(itest,jtest,k+1)*qonem,k=1,kk)
cdiag endif
c
      call xctilr(dpmixl( 1-nbdy,1-nbdy,  n),1, 1, 1,1, halo_ps)
c
      margin = 0  ! no horizontal derivatives
c
!$OMP PARALLEL DO PRIVATE(j)
!$OMP&             SHARED(m,n)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        call hybgenaj(m,n, j)
      enddo
!$OMP END PARALLEL DO
c
c --- vertical momentum flux across moving interfaces (the s-dot term in the
c --- momentum equation) - required to locally conserve momentum when hybgen
c --- moves vertical coordinates first, store old interface pressures in
c --- -pu-, -pv-
c
!$OMP PARALLEL DO PRIVATE(j,l,i,k)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do l=1,isu(j)
          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            pu(i,j,1)=0.0
            pu(i,j,2)=dpu(i,j,1,n)
          enddo
          do k=2,kk
            do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
              pu(i,j,k+1)=pu(i,j,k)+dpu(i,j,k,n)
            enddo
          enddo
        enddo
        do l=1,isv(j)
          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            pv(i,j,1)=0.0
            pv(i,j,2)=dpv(i,j,1,n)
          enddo
          do k=2,kk
            do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
              pv(i,j,k+1)=pv(i,j,k)+dpv(i,j,k,n)
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
c
c --- update layer thickness at -u,v- points.
      call dpudpv(dpu(1-nbdy,1-nbdy,1,n),
     &            dpv(1-nbdy,1-nbdy,1,n),
     &            p,depthu,depthv, margin)  ! p's halo extended by dpudpv
c
      if     (lpipe .and. lpipe_hybgen) then
c ---   compare two model runs.
        do k= 1,kk+1
          write (text,'(a9,i3)') 'p      k=',k
          call pipe_compare(p( 1-nbdy,1-nbdy,k),ip,text)
          write (text,'(a9,i3)') 'pu     k=',k
          call pipe_compare(pu(1-nbdy,1-nbdy,k),iu,text)
          write (text,'(a9,i3)') 'pv     k=',k
          call pipe_compare(pv(1-nbdy,1-nbdy,k),iv,text)
        enddo
        do k= 1,kk
          write (text,'(a9,i3)') 'dp     k=',k
          call pipe_compare(dp( 1-nbdy,1-nbdy,k,n),ip,text)
          write (text,'(a9,i3)') 'dpu    k=',k
          call pipe_compare(dpu(1-nbdy,1-nbdy,k,n),iu,text)
          write (text,'(a9,i3)') 'dpv    k=',k
          call pipe_compare(dpv(1-nbdy,1-nbdy,k,n),iv,text)
        enddo
      endif
c
!$OMP PARALLEL DO PRIVATE(j)
!$OMP&             SHARED(m,n)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        call hybgenbj(m,n, j)
      enddo
!$OMP END PARALLEL DO
c
      return
      end subroutine hybgen

      subroutine hybgenaj(m,n,j )
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n, j
c
c --- --------------------------------------------
c --- hybrid grid generator, single j-row (part A).
c --- --------------------------------------------
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
     &        qhrlx( kdm),         !relaxation coefficient, from qhybrlx
     &        dp0ij( kdm),         !minimum layer thickness
     &        dp0cum(kdm+1)        !minimum interface depth
      real    p_hat,p_hat0,p_hat2,p_hat3,hybrlx,
     &        delt,deltm,dels,delsm,q,qtr,qts,thkbop,
     &        zthk,dpthin
      integer i,k,ka,kp,ktr,l,fixlay,nums1d
ccc   real colint,colins,colout,colous
cdiag      character*12 cinfo
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
      include 'stmt_fns.h'
      qq(   delp,dp0)=max(qqmn, min(qqmx, delp/dp0))
      cushn(delp,dp0)=dp0*
     &                (1.0+cusha*(1.0-cushb*qq(delp,dp0))**2)*
     &                max(1.0, delp/(dp0*qqmx))
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
        do k=1,kk
          lcm(k) = .false.  !use same remapper for all layers
        enddo !k
      endif
c
      do 1 l=1,isp(j)
c
      do 2 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
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
cdiag if (itest.gt.0 .and. jtest.gt.0) then
cdiag   write (lp,'(a/(i6,1x,2f8.3,2f9.3,f9.3))')
cdiag.  'hybgen:   thkns  minthk     dpth  mindpth   hybrlx',
cdiag.  (k,dp(itest,jtest,k,n)*qonem,   dp0ij(k)*qonem,
cdiag.      p(itest,jtest,k+1)*qonem,dp0cum(k+1)*qonem,
cdiag.      1.0/qhrlx(k+1),
cdiag.   k=1,kk)
cdiag endif
c
c --- identify the deepest layer kp with significant thickness (> dpthin)
c
      kp = 2  !minimum allowed value
      do k=kk,3,-1
        if (p(i,j,k+1)-p(i,j,k).gt.dpthin) then
          kp=k
          exit
        endif
      enddo
c
c --- massless or near-massless (thickness <= dpthin) layers
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
        if (mxlmy) then
          qtr=p_hat/(p(i,j,k)-p(i,j,k-1))  !ok because <1.0 and >0.0
          q2( i,j,k-1,n)=q2( i,j,k-1,n)+
     &                     qtr*(q2( i,j,k,n)-q2( i,j,k-1,n))
          q2l(i,j,k-1,n)=q2l(i,j,k-1,n)+
     &                     qtr*(q2l(i,j,k,n)-q2l(i,j,k-1,n))
        endif
cdiag        if (i.eq.itest .and. j.eq.jtest) then
cdiag          write(lp,'(a,i3,f6.3,5f8.3)')
cdiag     &      'hybgen, 10(+):',
cdiag     &      k,q,temp(i,j,k,n),saln(i,j,k,n),
cdiag     &          th3d(i,j,k,n)+thbase,theta(i,j,k)+thbase
cdiag          call flush(lp)
cdiag        endif !debug
cdiag        if (i.eq.itest .and. j.eq.jtest) then
cdiag          write(lp,'(a,i3,f6.3,5f8.3)')
cdiag     &      'hybgen, 10(-):',
cdiag     &      k,0.0,temp(i,j,k,n),saln(i,j,k,n),
cdiag     &          th3d(i,j,k,n)+thbase,theta(i,j,k)+thbase
cdiag          call flush(lp)
cdiag        endif !debug
      endif !too light
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
c ---     top, bottom, thin and isopycnal layers remapped with PCM.
          lcm(k) = k.eq.1 .or. k.eq.kk .or.
     &             dprs(k).le.dpthin   .or.
     &             dprs(k).gt.dp0kp(k)
        endif !isopcm
      enddo !k
c
ccc   colint=temp(i,j,1,n)*(p(i,j,2)-p(i,j,1))
ccc   colins=saln(i,j,1,n)*(p(i,j,2)-p(i,j,1))
c
c --- try to restore isopycnic conditions by moving layer interfaces
c
      do 88 k=2,nhybrd
cdiag            if (i.eq.itest .and. j.eq.jtest) then
cdiag              write(cinfo,'(a9,i2.2,1x)') '  do 88 k=',k
cdiag 109          format (i9,2i5,a,a/(i9,8x,a,a,i3,f9.2,f8.2,f9.2,f8.2))
cdiag              write(lp,109) nstep,itest+i0,jtest+j0,
cdiag     .         cinfo,':    othkns  odpth    nthkns  ndpth',
cdiag     .        (nstep,cinfo,':',ka,
cdiag     .         (pres(ka+1)-
cdiag     .          pres(ka)   )*qonem,
cdiag     .          pres(ka+1)  *qonem,
cdiag     .         (p(itest,jtest,ka+1)-
cdiag     .          p(itest,jtest,ka)   )*qonem,
cdiag     .          p(itest,jtest,ka+1)  *qonem,ka=1,kk)
cdiag              call flush(lp)
cdiag            endif !debug
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
c --- if layer k-1 is too light, perform this adjustment every other pair
c --- of time steps beneath the z-coordinate domain. otherwise, layers
c --- k-1 and k may not be able to reach their target densities
      if (mod(nstep-1,4).le.1 .and.
     &    p(i,j,k).gt.dp0cum(k)+onem .and.
     &    th3d(i,j,k-1,n).lt.theta(i,j,k-1)) go to 8
c
cdiag      if (i.eq.itest .and. j.eq.jtest) then
cdiag        write(lp,'(a,f8.4)') 'hybgen, too dense:',
cdiag     &                       th3d(i,j,k,n)-theta(i,j,k)
cdiag        call flush(lp)
cdiag      endif !debug
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
      p_hat0=p_hat
      p_hat=max(p(i,j,k-1)+dp0ij(k-1),
     &      min(p(i,j,k+1)-dp0ij(k),
     &          p(i,j,k-1)+cushn(p_hat0-p(i,j,k-1),dp0ij(k-1)) ))
cdiag      if (i.eq.itest .and. j.eq.jtest) then
cdiag        write(lp,'(a,3f8.2)') 'hybgen, 9: ',
cdiag     &   (p_hat0-p(i,j,k-1))*qonem,
cdiag     &   cushn(p_hat0-p(i,j,k-1),dp0ij(k-1))*qonem,p_hat*qonem
cdiag        call flush(lp)
cdiag      endif !debug
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
cdiag          if (i.eq.itest .and. j.eq.jtest) then
cdiag            write(lp,'(a,i3.2,f8.2)') 'hybgen,  1blocking :',
cdiag     .                              k-1,p(i,j,k-1)*qonem
cdiag            call flush(lp)
cdiag          endif !debug
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
cdiag            if (i.eq.itest .and. j.eq.jtest) then
cdiag              write(lp,'(a,i3.2,f8.2)') 'hybgen,  2blocking :',
cdiag     .                                k-2,p(i,j,k-2)*qonem
cdiag              call flush(lp)
cdiag            endif !debug
            p_hat2=p(i,j,k-2)+cushn(p(i,j,k-1)-p_hat+p_hat0-p(i,j,k-2),
     &                              dp0ij(k-2))
            if (p_hat2.lt.p(i,j,k-1)-onemm) then
              p(i,j,k-1)=(1.0-qhrlx(k-1))*p(i,j,k-1) +
     &                    qhrlx(k-1) *max(p_hat2,2.0*p(i,j,k-1)-p_hat)
cdiag              if (i.eq.itest .and. j.eq.jtest) then
cdiag                write(lp,'(a,i3.2,f8.2)') 'hybgen,  3blocking :',
cdiag     .                                  k-1,p(i,j,k-1)*qonem
cdiag                call flush(lp)
cdiag              endif !debug
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
cdiag        if (i.eq.itest .and. j.eq.jtest) then
cdiag          write(lp,'(a,i3.2,28.2)') 'hybgen, entrain(k) :',
cdiag     .                            k,p_hat*qonem,p(i,j,k)*qonem
cdiag          call flush(lp)
cdiag        endif !debug
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
cdiag        if (i.eq.itest .and. j.eq.jtest) then
cdiag          write(lp,'(a,i3.2,2f8.2)') 'hybgen, entrain(k-):',
cdiag     .                            k,p_hat0*qonem,p(i,j,k)*qonem
cdiag          call flush(lp)
cdiag        endif !debug
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
c --- if layer k+1 is too dense, perform this adjustment every other pair
c --- of time steps beneath the z-coordinate domain. otherwise, layers
c --- k and k+1 may not be able to reach their target densities
      if (mod(nstep-1,4).ge.2 .and.
     &    p(i,j,k+1).gt.dp0cum(k+1)+onem   .and.
     &    th3d(i,j,k+1,n).gt.theta(i,j,k+1)     ) go to 8
c
cdiag      if (i.eq.itest .and. j.eq.jtest) then
cdiag        write(lp,'(a,f8.4)') 'hybgen, too light:',
cdiag     &                       theta(i,j,k)-th3d(i,j,k,n)
cdiag        call flush(lp)
cdiag      endif !debug
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
      if     (p(i,j,k+2).lt.p(i,j,kk+1)) then
        if     (p(i,j,kk+1)-p(i,j,k).gt.dp0ij(k)+dp0ij(k+1)) then
          p_hat=p(i,j,k+2)-cushn(p(i,j,k+2)-p_hat,dp0ij(k+1))
        endif
        p_hat=p(i,j,k)  +  max(p_hat    -p(i,j,k),dp0ij(k))
        p_hat=min(p_hat,
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
cdiag      if (i.eq.itest .and. j.eq.jtest) then
cdiag        write(lp,'(a,i3.2,f8.2)') 'hybgen, entrain(k+):',
cdiag   .                            k+1,p(i,j,k+1)*qonem
cdiag        call flush(lp)
cdiag      endif !debug
      endif !entrain
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
c --- grid onto the 'new' vertical grid.
c
      if     (lconserve) then  !usually .false.
        do ktr=1,ntracr+4
          asum(ktr,1) = 0.d0
          asum(ktr,2) = 0.d0
          asum(ktr,3) = 0.d0
        enddo !ktr
      endif
c
      prsf(1) = p(i,j,1)
      do k=1,kk
        prsf(k+1)   = p(i,j,k+1)
        dp(i,j,k,n) = prsf(k+1)-prsf(k)
      enddo
      if     (hybmap.eq.0) then !PCM
        call hybgen_pcm_remap(s1d,pres,dpi,
     &                        f1d,prsf,kk,nums1d,dpthin)
      elseif (hybmap.eq.1) then !PLM
        call hybgen_plm_coefs(s1d,     dpi,lcm,c1d,
     &                                 kk,nums1d,dpthin)
        call hybgen_plm_remap(s1d,pres,dpi,    c1d,
     &                        f1d,prsf,kk,nums1d,dpthin)
      elseif (hybmap.eq.3) then !PLM (as in 2.1.08)
        call hybgen_plm_coefx(s1d,     dpi,lcm,c1d,
     &                                 kk,nums1d,dpthin)
        call hybgen_plm_remap(s1d,pres,dpi,    c1d,
     &                        f1d,prsf,kk,nums1d,dpthin)
      elseif (hybmap.eq.2) then !PPM
        call hybgen_ppm_coefs(s1d,     dpi,lcm,c1d,
     &                                 kk,nums1d,dpthin)
        call hybgen_ppm_remap(s1d,pres,dpi,    c1d,
     &                        f1d,prsf,kk,nums1d,dpthin)
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
cdiag if (i.eq.itest .and. j.eq.jtest) then
cdiag  if     (hybflg.eq.0) then  !T&S
cdiag   write (lp,103) nstep,itest+i0,jtest+j0,
cdiag&  '    hybgen, do 22:  temp    saln    dens     thkns     dpth',
cdiag&  (k,s1d(k,1),s1d(k,2),0.0,
cdiag&   (pres(k+1)-pres(k))*qonem,pres(k+1)*qonem,
cdiag&   k,temp(i,j,k,n),saln(i,j,k,n),
cdiag&   th3d(i,j,k,n)+thbase,dp(i,j,k,n)*qonem,
cdiag&   p(i,j,k+1)*qonem,
cdiag&  k=1,kk)
cdiag  elseif (hybflg.eq.1) then  !th&S
cdiag   write (lp,103) nstep,itest+i0,jtest+j0,
cdiag&  '    hybgen, do 22:  temp    saln    dens     thkns     dpth',
cdiag&  (k,0.0,s1d(k,2),s1d(k,1)+thbase,
cdiag&   (pres(k+1)-pres(k))*qonem,pres(k+1)*qonem,
cdiag&   k,temp(i,j,k,n),saln(i,j,k,n),
cdiag&   th3d(i,j,k,n)+thbase,dp(i,j,k,n)*qonem,
cdiag&   p(i,j,k+1)*qonem,
cdiag&  k=1,kk)
cdiag  elseif (hybflg.eq.2) then  !th&T
cdiag   write (lp,103) nstep,itest+i0,jtest+j0,
cdiag&  '    hybgen, do 22:  temp    saln    dens     thkns     dpth',
cdiag&  (k,s1d(k,2),0.0,s1d(k,1)+thbase,
cdiag&   (pres(k+1)-pres(k))*qonem,pres(k+1)*qonem,
cdiag&   k,temp(i,j,k,n),saln(i,j,k,n),
cdiag&   th3d(i,j,k,n)+thbase,dp(i,j,k,n)*qonem,
cdiag&   p(i,j,k+1)*qonem,
cdiag&  k=1,kk)
cdiag  endif
cdiag   call flush(lp)
cdiag endif !debug
c
 2    continue  !i
c
c --- to avoid roundoff errors in -dpudpv- after restart, make sure p=p(dp)
c
      do 1 k=1,kk
      do 1 i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
 1    continue !i;k;l
c
      return
      end subroutine hybgenaj

      subroutine hybgenbj(m,n, j)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'common_blocks.h'
c
      integer m,n, j
c
c --- --------------------------------------------
c --- hybrid grid generator, single j-row (part B).
c --- --------------------------------------------
c
      integer i,k,l
      logical lcm(kdm)      !use PCM for some layers?
      real    s1d(kdm),     !original scalar fields
     &        f1d(kdm),     !final    scalar fields
     &        c1d(kdm,3),   !interpolation coefficients
     &        dpi( kdm),    !original layer thicknesses, >= dpthin
     &        dprs(kdm),    !original layer thicknesses
     &        pres(kdm+1),  !original layer interfaces
     &        prsf(kdm+1)   !final    layer interfaces
      real    dpthin
c
c --- vertical momentum flux across moving interfaces (the s-dot term in the
c --- momentum equation) - required to locally conserve momentum when hybgen
c --- moves vertical coordinates.
c
      dpthin = 0.001*onemm
c
      do 412 l=1,isu(j)
      do 412 i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
c
c --- store one-dimensional arrays of -u- and -p- for the 'old' vertical grid
      pres(1)=pu(i,j,1)
      do k=1,kk
        s1d(k)   =u(i,j,k,n)
        pres(k+1)=pu(i,j,k+1)
        dprs(k)  =pres(k+1)-pres(k)
        dpi( k)  =max(dprs(k),dpthin)
c       
        if     (isopcm) then
c ---     top, bottom, thin and isopycnal layers remapped with PCM.
          lcm(k) = k.eq.1 .or. k.eq.kk .or.
     &             dprs(k).le.dpthin   .or.
     &             dprs(k).gt.dp0kp(k)
        endif !isopcm
      enddo !k
c
c --- remap -u- profiles from the 'old' vertical grid onto the
c --- 'new' vertical grid.
c
      prsf(1) = pu(i,j,1)
      do k=1,kk
        pu(i,j,k+1) = pu(i,j,k) + dpu(i,j,k,n)  !new vertical grid
        prsf(k+1)   = pu(i,j,k+1)
      enddo
      if     (hybmap.eq.0) then !PCM
        call hybgen_pcm_remap(s1d,pres,dpi,
     &                        f1d,prsf,kk,1,dpthin)
      elseif (hybmap.eq.1) then !PLM
        call hybgen_plm_coefs(s1d,     dpi,lcm,c1d,
     &                                 kk,1,dpthin)
        call hybgen_plm_remap(s1d,pres,dpi,    c1d,
     &                        f1d,prsf,kk,1,dpthin)
      elseif (hybmap.eq.2) then !PPM
        call hybgen_ppm_coefs(s1d,     dpi,lcm,c1d,
     &                                 kk,1,dpthin)
        call hybgen_ppm_remap(    pres,dpi,    c1d,
     &                        f1d,prsf,kk,1,dpthin)
      endif
      do k=1,kk
        u(i,j,k,n) = f1d(k)
      enddo !k
c
 104  format (i9,2i5,a/(33x,i3,f8.3,f9.3,f9.2))
cdiag if (i.eq.itest .and. j.eq.jtest) then
cdiag   write (lp,104) nstep,itest+i0,jtest+j0,
cdiag&  '   hybgen, do 412:  u       thkns     dpth',
cdiag&  (k,s1d(k,1),
cdiag&   (pres(k+1)-pres(k))*qonem,pres(k+1)*qonem,
cdiag&   k,u(i,j,k,n),
cdiag&   dpu(i,j,k,n)*qonem,pu(i,j,k+1)*qonem,
cdiag&   k=1,kk)
cdiag endif !debug
c
 412  continue !l;i
c
      do 512 l=1,isv(j)
      do 512 i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
c
c --- store one-dimensional arrays of -v- and -p- for the 'old' vertical grid
      pres(1)=pv(i,j,1)
      do k=1,kk
        s1d(k)   =v(i,j,k,n)
        pres(k+1)=pv(i,j,k+1)
        dprs(k)  =pres(k+1)-pres(k)
        dpi( k)  =max(dprs(k),dpthin)
c       
        if     (isopcm) then
c ---     top, bottom, thin and isopycnal layers remapped with PCM.
          lcm(k) = k.eq.1 .or. k.eq.kk .or.
     &             dprs(k).le.dpthin   .or.
     &             dprs(k).gt.dp0kp(k)
        endif !isopcm
      enddo !k
c
c --- remap -v- profiles from the 'old' vertical grid onto the
c --- 'new' vertical grid.
c
      prsf(1) = pv(i,j,1)
      do k=1,kk
        pv(i,j,k+1) = pv(i,j,k) + dpv(i,j,k,n)  !new vertical grid
        prsf(k+1)   = pv(i,j,k+1)
      enddo
      if     (hybmap.eq.0) then !PCM
        call hybgen_pcm_remap(s1d,pres,dpi,
     &                        f1d,prsf,kk,1,dpthin)
      elseif (hybmap.eq.1) then !PLM
        call hybgen_plm_coefs(s1d,     dpi,lcm,c1d,
     &                                 kk,1,dpthin)
        call hybgen_plm_remap(s1d,pres,dpi,    c1d,
     &                        f1d,prsf,kk,1,dpthin)
      elseif (hybmap.eq.2) then !PPM
        call hybgen_ppm_coefs(s1d,     dpi,lcm,c1d,
     &                                 kk,1,dpthin)
        call hybgen_ppm_remap(    pres,dpi,    c1d,
     &                        f1d,prsf,kk,1,dpthin)
      endif
      do k=1,kk
        v(i,j,k,n) = f1d(k)
      enddo !k
c
cdiag if (i.eq.itest .and. j.eq.jtest) then
cdiag   write (lp,104) nstep,itest+i0,jtest+j0,
cdiag&  '   hybgen, do 512:  v       thkns     dpth',
cdiag&  (k,s1d(k,1),
cdiag&   (pres(k+1)-pres(k))*qonem,pres(k+1)*qonem,
cdiag&   k,v(i,j,k,n),
cdiag&   dpv(i,j,k,n)*qonem,pv(i,j,k+1)*qonem,
cdiag&   k=1,kk)
cdiag endif !debug
c
 512  continue !l;i
c
      return
      end subroutine hybgenbj

      subroutine hybgen_pcm_remap(si,pi,dpi,
     &                            so,po,kk,ks,thin)
      implicit none
c
      integer kk,ks
      real    si(kk,ks),pi(kk+1),dpi(kk),
     &        so(kk,ks),po(kk+1),thin
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
c                  pi(kk+1) is the bathymetry
c                  pi(k+1) >= pi(k)
c       dpi   - initial layer thicknesses (>= thin)
c       kk    - number of layers
c       ks    - number of fields
c       po    - target interface depths (non-negative)
c                  po(   1) is the surface
c                  po(kk+1) is the bathymetry (== pi(kk+1))
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
      integer i,k,l,lf
      real    q,zb,zt,o(ks)
      real*8  sok(ks)
c
      lf=1
      zb=po(1)
      do k= 1,kk
        zt = zb
        zb = po(k+1)
*       write(lp,*) 'k,zt,zb = ',k,zt,zb
        if     (zb-zt.le.thin .or. zt.ge.pi(kk+1)) then
c
c ---     thin or bottomed layer, values taken from layer above
c
          do i= 1,ks
            so(k,i) = so(k-1,i)
          enddo !i
        else
c
c         form layer averages.
c
*         if     (pi(lf).gt.zt) then
*           write(lp,*) 'bad lf = ',lf
*           stop
*         endif
          do i= 1,ks
            sok(i) = 0.0
              o(i) = si(lf,i)  !offset to reduce round-off
          enddo !i
          do l= lf,kk
            if     (pi(l).gt.zb) then
*             write(lp,*) 'l,lf= ',l,lf,l-1
              lf = l-1
              exit
            elseif (pi(l).ge.zt .and. pi(l+1).le.zb) then
c
c             the input layer is completely inside the output layer
c
              q   = dpi(l)/(zb-zt)
              do i= 1,ks
                sok(i) = sok(i) + q*(si(l,i) - o(i))
              enddo !i
*             write(lp,*) 'L,q = ',l,q
            else
c
c             the input layer is partially inside the output layer
c
              q   = max(min(pi(l+1),zb)-max(pi(l),zt),thin)/(zb-zt)
              do i= 1,ks
                sok(i) = sok(i) + q*(si(l,i) - o(i))
              enddo !i
*             write(lp,*) 'l,q = ',l,q
            endif
          enddo !l
          do i= 1,ks
            so(k,i) = o(i) + sok(i)
          enddo !i
        endif  !thin:else
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
c     method: piecewise linear across each input cell
c             the output is the average of the interpolation
c             profile across each output cell.
c     minimod limiter, followed by a discontinuity minimization pass.
c
c  2) input arguments:
c       si    - initial scalar fields in pi-layer space
c       dpi   - initial layer thicknesses (>= thin)
c       lc    - use PCM for selected layers
c       kk    - number of layers
c       ks    - number of fields
c       thin  - layer thickness (>0) that can be ignored
c
c  3) output arguments:
c       ci    - coefficents (slopes) for hybgen_plm_remap
c
c  4) Tim Campbell, Mississippi State University, October 2002.
c     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
c-----------------------------------------------------------------------
c
      integer k,i
      real    ql,qr,sl,sr
c
      do i= 1,ks
        ci(1, i) = 0.0  !layer 1  is PCM
        ci(kk,i) = 0.0  !layer kk is PCM
      enddo !i
      do k= 2,kk-1
        if     (lc(k) .or. dpi(k).le.thin) then  !use PCM
          do i= 1,ks
            ci(k,i) = 0.0
          enddo !i
        else
          ql = 2.0*dpi(k)/(dpi(k-1)+dpi(k))
          qr = 2.0*dpi(k)/(dpi(k+1)+dpi(k))
          do i= 1,ks
c ---       PLM (non-zero slope, but no new extrema)
c ---       layer value is si-0.5*ci at top    interface,
c ---                  and si+0.5*ci at bottom interface.
c
c ---       nonuniform minimod limiter.  For a discussion of PLM limiters, see
c ---       Finite Volume Methods for Hyperbolic Problems by R.J. Leveque.
            sl=ql*(si(k,  i)-si(k-1,i))
            sr=qr*(si(k+1,i)-si(k,  i))
            if (sl*sr.gt.0.) then
              ci(k,i)=sign(min(abs(sl),abs(sr)),sl)
            else
              ci(k,i)=0.0  !local extrema, so no slope
            endif
          enddo !i
        endif  !PCM:PLM
      enddo !k
c
c --- Minimize discontinuities between zones
c
      do k= 2,kk-1
        do i= 1,ks
          if (ci(k,i).ne.0.0) then
            ql=(si(k,  i)-si(k-1,i)) - 0.5*(ci(k,  i)+ci(k-1,i))
            qr=(si(k+1,i)-si(k,  i)) - 0.5*(ci(k+1,i)+ci(k,  i))
            ci(k,i)=ci(k,i)+2.0*sign(min(abs(ql),abs(qr)),ql)
          endif
        enddo !i
      enddo !k
      return
      end subroutine hybgen_plm_coefs

      subroutine hybgen_plm_coefx(si,dpi,lc,ci,kk,ks,thin)
      implicit none
c
      integer kk,ks
      logical lc(kk)
      real    si(kk,ks),dpi(kk),ci(kk,ks),thin
c
c-----------------------------------------------------------------------
c  1) coefficents for remaping from one set of vertical cells to another.
c     method: piecewise linear across each input cell
c             the output is the average of the interpolation
c             profile across each output cell.
c     2.2.08 version: that uses monotonized central-difference limiter.
c
c  2) input arguments:
c       si    - initial scalar fields in pi-layer space
c       dpi   - initial layer thicknesses (>= thin)
c       lc    - use PCM for selected layers
c       kk    - number of layers
c       ks    - number of fields
c       thin  - layer thickness (>0) that can be ignored
c
c  3) output arguments:
c       ci    - coefficents (slopes) for hybgen_plm_remap
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
          qcen = dpi(k)/(dpi(k)+0.5*(dpi(k-1)+dpi(k+1)))
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
      end subroutine hybgen_plm_coefx

      subroutine hybgen_plm_remap(si,pi,dpi,ci,
     &                            so,po,kk,ks,thin)
      implicit none
c
      integer kk,ks
      real    si(kk,ks),pi(kk+1),dpi(kk),ci(kk,ks),
     &        so(kk,ks),po(kk+1),thin
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
c                  pi(kk+1) is the bathymetry
c                  pi(k+1) >= pi(k)
c       dpi   - initial layer thicknesses (>= thin)
c       ci    - coefficents (slopes) from hybgen_plm_coefs
c       kk    - number of layers
c       ks    - number of fields
c       po    - target interface depths (non-negative)
c                  po(   1) is the surface
c                  po(kk+1) is the bathymetry (== pi(kk+1))
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
      integer i,k,l,lf,lo
      real    q,qc,zb,zc,zt,o(ks)
      real*8  sok(ks)
c
      lf=1
      zb=po(1)
      do k= 1,kk
        zt = zb
        zb = po(k+1)
*       write(lp,*) 'k,zt,zb = ',k,zt,zb
        if     (zb-zt.le.thin .or. zt.ge.pi(kk+1)) then
c
c ---     thin or bottomed layer, values taken from layer above
c
          do i= 1,ks
            so(k,i) = so(k-1,i)
          enddo !i
        else
c
c         form layer averages.
c
*         if     (pi(lf).gt.zt) then
*           write(lp,*) 'bad lf = ',lf
*           stop
*         endif
          do i= 1,ks
            sok(i) = 0.0
              o(i) = si(lf,i)  !offset to reduce round-off
          enddo !i
          do l= lf,kk
            if     (pi(l).gt.zb) then
c
c             the input layer is below the output layer
c
*             write(lp,*) 'l,lf= ',l,lf,l-1
              lf = l-1
              exit
            elseif (pi(l).ge.zt .and. pi(l+1).le.zb) then
c
c             the input layer is completely inside the output layer
c
              q   = dpi(l)/(zb-zt)
              do i= 1,ks
                sok(i) = sok(i) + q*(si(l,i)-o(i))
              enddo !i
*             write(lp,*) 'L,q = ',l,q
            else
c
c             the input layer is partially inside the output layer
c             average of linear profile is its center value
c
              q   = max( min(pi(l+1),zb)-max(pi(l),zt), thin )/(zb-zt)
              zc  = 0.5*(min(pi(l+1),zb)+max(pi(l),zt))
              qc  = (zc-pi(l))/dpi(l) - 0.5
              do i= 1,ks
                sok(i) = sok(i) + q*(si(l,i) - o(i) + qc*ci(l,i))
              enddo !i
*             write(lp,*) 'l,q = ',l,q
            endif
          enddo !l
          do i= 1,ks
            so(k,i) = o(i) + sok(i)
          enddo !i
        endif  !thin:else
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
c       dp    - initial layer thicknesses (>= thin)
c       lc    - use PCM for selected layers
c       kk    - number of layers
c       ks    - number of fields
c       thin  - layer thickness (>0) that can be ignored
c
c  3) output arguments:
c       ci    - coefficents for hybgen_ppm_remap
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
      do j=1,kk
         dpq3( j) = dp(j)/(dp(j-1)+dp(j)+dp(j+1))
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
          if ((s(j+1,i)-s(j,i))*(s(j,i)-s(j-1,i)).le.0.) then !local extremum
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
          ci(j,i,1)=al(j)
          ci(j,i,2)=ar(j)-al(j)
          ci(j,i,3)=6.0*s(j,i)-3.0*(al(j)+ar(j))
        enddo !j
      enddo !i
      return
      end subroutine hybgen_ppm_coefs

      subroutine hybgen_ppm_remap(si,pi,dpi,ci,
     &                            so,po,kk,ks,thin)
      implicit none
c
      integer kk,ks
      real    si(kk,ks),pi(kk+1),dpi(kk),ci(kk,ks,3),
     &        so(kk,ks),po(kk+1),thin
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
c                  pi(kk+1) is the bathymetry
c                  pi(k+1) >= pi(k)
c       dpi   - initial layer thicknesses (>= thin)
c       ci    - coefficents from hybgen_ppm_coefs
c       kk    - number of layers
c       ks    - number of fields
c       po    - target interface depths (non-negative)
c                  po(   1) is the surface
c                  po(kk+1) is the bathymetry (== pi(kk+1))
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
      real    xb,xt,zb,zt,o
      real*8  sz
c
      lb=1
      zb=po(1)
      do k= 1,kk
        zt = zb
        zb = po(k+1)
*       write(lp,*) 'k,zt,zb = ',k,zt,zb
        if     (zb-zt.le.thin .or. zt.ge.pi(kk+1)) then
c
c ---     thin or bottomed layer, values taken from layer above
c
          do i= 1,ks
            so(k,i) = so(k-1,i)
          enddo !i
        else
c
c         form layer averages.
c
*         if     (pi(lb).gt.zt) then
*           write(lp,*) 'bad lb = ',lb
*           stop
*         endif
          lt=lb !top will always correspond to bottom of previous
          lb=lt !find input layer containing bottom output interface
          do while (pi(lb+1).lt.zb.and.lb.lt.kk)
            lb=lb+1
          enddo
          xt=(zt-pi(lt))/dpi(lt)
          xb=(zb-pi(lb))/dpi(lb)
          do i= 1,ks
            o = si((lt+lb)/2,i)  !offset to reduce round-off
            if     (lt.ne.lb) then
              sz=   dpi(lt)*(     (ci(lt,i,1)-o)*(1.-xt)
     &                       +0.5*(ci(lt,i,2)+
     &                             ci(lt,i,3) ) *(1.-xt**2)
     &                            -ci(lt,i,3)   *(1.-xt**3)/3.0 )
              do l=lt+1,lb-1
                sz=sz+dpi(l)*(si(l,i)-o)
              enddo !l
              sz=sz+dpi(lb)*(     (ci(lb,i,1)-o)*    xb
     &                       +0.5*(ci(lb,i,2)+
     &                             ci(lb,i,3) ) *    xb**2
     &                            -ci(lb,i,3)   *    xb**3 /3.0 )
            else
              sz=dpi(lt)*(     (ci(lt,i,1)-o)*(xb-xt)
     &                    +0.5*(ci(lt,i,2)+
     &                          ci(lt,i,3) ) *(xb**2-xt**2)
     &                         -ci(lt,i,3)   *(xb**3-xt**3)/3.0 )
            endif
            so(k,i) = o + sz/(zb-zt)
          enddo !i
        endif !thin:std layer
      enddo !k
      return
      end subroutine hybgen_ppm_remap

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
c> Aug. 2007 -- added hybmap and isopcm for PCM,PLM,PPM remaper selection

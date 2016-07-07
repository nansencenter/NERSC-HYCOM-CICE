c-----------------------------------------------------------------------------
      real sig,dsigdt,dsigds,tofsig,sofsig,kappaf,kappaf1,kappafs,
     &     sigloc,dsiglocdt,dsiglocds,tofsigloc
c
      real    a0,a1,a2,cubr,cubq,cuban,cubrl,cubim
      real    c1p,c2p,c3p,c4p,c5p,c6p,c7p
      real    r,s,t,prs,ylat
      integer kkf
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
      real       pref
      parameter (pref=0.0)
c
c --- coefficients for sigma-2 (based on Brydon & Sun fit)
csig2 real       c1,c2,c3,c4,c5,c6,c7
csig2 parameter (c1= 9.77093E+00, c2=-2.26493E-02, c3= 7.89879E-01,
csig2&           c4=-6.43205E-03, c5=-2.62983E-03,
csig2&           c6= 2.75835E-05, c7= 3.15235E-05)
csig2 real       pref
csig2 parameter (pref=2000.e4)
c
c --- coefficients for kappa^(theta)
c --- new values (w.r.t. t-toff,s-soff,prs) from Shan Sun, 2/8/01
c --- southern (1) and northern (2) hemispheres
      real, parameter, dimension(2) ::
     &  toff = (/  0.0,          3.0         /)
     & ,soff = (/ 34.0,         35.0         /)
     & ,qt   = (/ -2.89196E-01, -2.61829E-01 /)
     & ,qs   = (/ -1.08670E-01, -1.05131E-01 /)
     & ,qtt  = (/  4.56626E-03,  4.29277E-03 /)
     & ,qst  = (/  7.90504E-04,  7.71097E-04 /)
     & ,qttt = (/ -3.03869E-05, -3.03869E-05 /)
     & ,qpt  = (/  1.07106E-09,  1.00638E-09 /)
     & ,qpst = (/  1.41542E-11,  1.48599E-11 /)
     & ,qptt = (/ -1.31384E-11, -1.31384E-11 /)
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
c
c --- thermobaric compressibility coefficient (integral from prs to pref)
c --- Sun et.al. (1999) JPO 29 pp 2719-2729
c --- kappaf1 used internally to simplify offsetting T and S.
c --- always invoke via kappaf.
c --- t: potential temperature; s: psu; prs: pressure; ylat: latitude
      kappafs(ylat)=max(0.0,min(1.0,(ylat+30.0)/60.0))  !0 at 30S, 1 at 30N
      kappaf1(t,s,prs,kkf)=(1.e-11*qthref)*(prs-pref)*
     &  ( s*( qs(kkf)+t* qst(kkf) ) +
     &    t*( qt(kkf)+t*(qtt(kkf)+t*qttt(kkf))+
     &        0.5*(prs+pref)*(qpt(kkf)+s*qpst(kkf)+t*qptt(kkf)) ) )
      kappaf(t,s,prs,ylat)=
     &     (1.0-kappafs(ylat))*
     &     kappaf1(max(-2.0,min(32.0,t))-toff(1),
     &             max(30.0,min(38.0,s))-soff(1),
     &             prs,1) +
     &          kappafs(ylat) *
     &     kappaf1(max(-2.0,min(32.0,t))-toff(2),
     &             max(30.0,min(38.0,s))-soff(2),
     &             prs,2)
c
c --- locally referenced sigma, a fit towards Jackett & McDougall (1995)
c --- t: potential temperature; s: psu; prs: pressure
      c1p(prs)=alphap(1)+1.e-5*prs*(betap(1)+1.e-5*prs*gammap(1))
      c2p(prs)=alphap(2)+1.e-5*prs*(betap(2)+1.e-5*prs*gammap(2))
      c3p(prs)=alphap(3)+1.e-5*prs*(betap(3)+1.e-5*prs*gammap(3))
      c4p(prs)=alphap(4)+1.e-5*prs*(betap(4)+1.e-5*prs*gammap(4))
      c5p(prs)=alphap(5)+1.e-5*prs*(betap(5)+1.e-5*prs*gammap(5))
      c6p(prs)=alphap(6)+1.e-5*prs*(betap(6)+1.e-5*prs*gammap(6))
      c7p(prs)=alphap(7)+1.e-5*prs*(betap(7)+1.e-5*prs*gammap(7))
      sigloc(t,s,prs)=c1p(prs)+c3p(prs)*s+
     &       t*(c2p(prs)+c5p(prs)*s+t*(c4p(prs)+c7p(prs)*s+c6p(prs)*t))
      dsiglocdt(t,s,prs)=(c2p(prs)+c5p(prs)*s+
     &       2.*t*(c4p(prs)+c7p(prs)*s+1.5*c6p(prs)*t))
      dsiglocds(t,s,prs)=(c3p(prs)+t*(c5p(prs)+t*c7p(prs)))
c
c> Revision history
c>
c> May  2000 - conversion to SI units
c> Jul  2000 - removed rarely used functions, constants via parameter
c> Jan  2002 - removed geometery functions
c> Dec  2002 - new thermobaricity fit with toff=0.0,soff=34.0
c> Jun  2003 - removed sigma4
c> Jun  2003 - two thermobaricity reference states, based on latitude
c> Jun  2003 - added locally referenced sigma
c-----------------------------------------------------------------------------

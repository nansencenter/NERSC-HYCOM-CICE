c-----------------------------------------------------------------------------
      real sig,dsigdt,dsigds,tofsig,sofsig,kappaf,kappaf1
c
      real a0,a1,a2,cubr,cubq,cuban,cubrl,cubim
      real r,s,t,prs
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
c --- coefficients for sigma-4 (based on Brydon & Sun fit)
csig4 real       c1,c2,c3,c4,c5,c6,c7
csig4 parameter (c1= 1.92362E+01, c2=-8.82080E-02, c3= 7.73552E-01,
csig4&           c4=-5.46858E-03, c5=-2.31866E-03,
csig4&           c6= 2.11306E-05, c7= 2.82474E-05)
csig4 real       pref
csig4 parameter (pref=4000.e4)
c
c --- coefficients for kappa^(theta)
c --- new values (w.r.t. t-toff,s-soff,prs-poff) from Shan Sun, 2/8/01
      real       toff,soff,poff,prefo
      parameter (toff=3.0,soff=35.0,poff=2000.e4,prefo=pref-poff)
c
      real       qt,qs,qtt,qst,qttt,qpt,qpst,qptt
      parameter (qt  =-2.41701E-01,  qs  =-1.05131E-01,
     &           qtt = 4.03001E-03,  qst = 1.06829E-03,
     &           qttt=-3.03869E-05,  qpt = 1.00638E-09,
     &           qpst= 1.48599E-11,  qptt=-1.31384E-11)
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
c --- kappaf1 used internally to simplify offsetting T,S,prs.
c --- always invoke via kappaf.
      kappaf1(t,s,prs)=(1.e-11*qthref)*(prs-prefo)*
     &  ( s*( qs+t* qst ) +
     &    t*( qt+t*(qtt+t*qttt)+
     &        0.5*(prs+prefo)*(qpt+s*qpst+t*qptt) ) )
      kappaf(t,s,prs)=kappaf1(max(-2.0,min(32.0,t))-toff,
     &                        max(30.0,min(38.0,s))-soff,
     &                        prs-poff)
c
c> Revision history
c>
c> May  2000 - conversion to SI units
c> Jul  2000 - removed rarely used functions, constants via parameter
c> Jan  2002 - removed geometery functions
c> Aug  2002 - new thermobaricity fit with toff=3.0
c-----------------------------------------------------------------------------

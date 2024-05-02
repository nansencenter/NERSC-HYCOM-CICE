c-----------------------------------------------------------------------------
      integer, parameter ::
     &  sigver=4  !9-term sigma-2
csig0&  sigver=3  !9-term sigma-0
c
      real    sig,dsigdt,dsigds,tofsig,sofsig,kappaf,
     &        sigloc,dsiglocdt,dsiglocds,tofsigloc
c
      real    kappaf1
      real    sofsig_a,sofsig_b,sofsig_c
      real    a0,a1,a2,cubr,cubq,cuban,cubrl,cubim
      real    c1l,c2l,c3l,c4l,c5l,c6l,c7l
c
      real    r,s,t,prs
      integer kkf
c
      real, parameter ::
     &   ahalf=1.0/2.0,
     &   a3rd =1.0/3.0, athird =a3rd,
     &   a4th =1.0/4.0, afourth=a4th
c
c --- sigma-theta as a function of temp (deg c) and salinity (psu)
c --- (9-term polynomial fit to T:[-2:30],S:[18:38])
c
c --- coefficients for sigma-0.
csig0 real, parameter ::
csig0&   c1=-4.311829E-02,  !const. coefficent
csig0&   c2= 5.429948E-02,  !T      coefficent
csig0&   c3= 8.011774E-01,  !   S   coefficent
csig0&   c4=-7.641336E-03,  !T^2    coefficent
csig0&   c5=-3.258442E-03,  !T  S   coefficent
csig0&   c6= 3.757643E-05,  !T^3    coefficent
csig0&  rc6=1.0/c6,
csig0&   c7= 3.630361E-05,  !T^2S   coefficent
csig0&   c8= 8.675546E-05,  !   S^2 coefficent
csig0&   c9= 3.995086E-06,  !T  S^2 coefficent
csig0& pref= 0.0           ! reference pressure, Pascals
c --- coefficients for sigma-2.
      real, parameter ::
     &   c1= 9.903308E+00,  !const. coefficent
     &   c2=-1.618075E-02,  !T      coefficent
     &   c3= 7.819166E-01,  !   S   coefficent
     &   c4=-6.593939E-03,  !T^2    coefficent
     &   c5=-2.896464E-03,  !T  S   coefficent
     &   c6= 3.038697E-05,  !T^3    coefficent
     &  rc6= 1.0/c6,
     &   c7= 3.266933E-05,  !T^2S   coefficent
     &   c8= 1.180109E-04,  !   S^2 coefficent
     &   c9= 3.399511E-06,  !T  S^2 coefficent
     & pref= 2000.0e4       ! reference pressure, Pascals
c
c --- HYCOM pressure to bar, for locally referenced equations
      real, parameter :: prs2pb=1.e-5       !Pascals to bar
c
c --- coefficients for kappa^(theta)
c --- new values (w.r.t. t-toff,s-soff,prs) from Shan Sun, Sep.2004
c --- 1=Arctic/Antarctic; 2=Atlantic; 3=Mediterranean
      real, parameter ::
     &   rhoref=1.e3  !rhoref=qthref kg/m^3
      real, parameter ::
     &   sclkap=1.e-11
      real, parameter, dimension(3) ::
     &  toff = (/  0.0,             3.0,            13.0 /)
     & ,soff = (/ 34.5,            35.0,            38.5 /)
     & ,qttt = (/ -3.03869354E-05, -3.03869352E-05, -3.03869353E-05 /)
     & ,qtt  = (/  4.56625601E-03,  4.29277358E-03,  3.38116552E-03 /)
     & ,qt   = (/ -2.88801209E-01, -2.61828868E-01, -1.81335007E-01 /)
     & ,qs   = (/ -1.08670290E-01, -1.05131061E-01, -9.33336309E-02 /)
     & ,qst  = (/  7.90503772E-04,  7.71096940E-04,  1.07270585E-03 /)
     & ,qpt  = (/  1.07813750E-09,  1.00638435E-09,  7.57239852E-10 /)
     & ,qpst = (/  1.41541548E-11,  1.48598578E-11,  3.89226107E-12 /)
     & ,qptt = (/ -1.31383708E-11, -1.31383707E-11, -1.31383708E-11 /)
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
c --- auxiliary statements for finding root of cubic polynomial
      a0(s,r)=(c1+s*(c3+s*c8)-r)*rc6  !constant  coefficient
      a1(s)  =(c2+s*(c5+s*c9)  )*rc6  !linear    coefficient
      a2(s)  =(c4+s* c7        )*rc6  !quadratic coefficient
                                      !cubic     coefficient is c6*rc6=1.0
      cubq(s)  =a3rd*     a1(s)                   -(a3rd*a2(s))**2
      cubr(r,s)=a3rd*(0.5*a1(s)*a2(s)-1.5*a0(s,r))-(a3rd*a2(s))**3
c --- if q**3+r**2>0, water is too dense to yield real root at given
c --- salinitiy. setting q**3+r**2=0 in that case is equivalent to
c --- lowering sigma until a double real root is obtained.
      cuban(r,s)=a3rd*atan2(sqrt(max(0.0,-(cubq(s)**3+cubr(r,s)**2))),
     &                      cubr(r,s))
      cubrl(r,s)=sqrt(-cubq(s))*cos(cuban(r,s))
      cubim(r,s)=sqrt(-cubq(s))*sin(cuban(r,s))
c
c --- -----------------
c --- equation of state
c --- -----------------
c
c --- sigma-theta as a function of temp (deg c) and salinity (psu)
c --- (polynomial fit that is cubic in T and quadratic in S)
c
      sig(t,s)=(c1+s*(c3+s* c8)+
     &             t*(c2+s*(c5+s*c9)+t*(c4+s*c7+t*c6)))
c
c --- d(sig)/dt
      dsigdt(t,s)=(c2+s*(c5+s*c9)+2.0*t*(c4+s*c7+1.5*t*c6))
c
c --- d(sig)/ds
      dsigds(t,s)=(c3+t*(c5+t*c7)+2.0*s* c8)
c
c --- temp (deg c) as a function of sigma and salinity (psu)
c --- find a cubic polynominal root of t**3+a2*t**2+a1*t+a0=0
      tofsig(r,s)=-cubrl(r,s)+sqrt(3.0)*cubim(r,s)-a3rd*a2(s)
c
c --- salinity (psu) as a function of sigma and temperature (deg c)
c --- find a quadratic polynominal root of a*s**2+b*s+c=0
      sofsig_a(r,t)=(c8+t* c9)                !quadratic coefficient
      sofsig_b(r,t)=(c3+t*(c5+t* c7))         !linear    coefficient
      sofsig_c(r,t)=(c1+t*(c2+t*(c4+t*c6))-r) !constant  coefficient
      sofsig(r,t)=(2.0*sofsig_c(r,t))/
     &            (-sofsig_b(r,t)
     &             -sign(sqrt(sofsig_b(r,t)**2-
     &                        4.0*sofsig_a(r,t)*sofsig_c(r,t)),
     &                   sofsig_b(r,t)))
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
c --- thermobaric compressibility coefficient (integral from prs to pref)
c ---     Sun et.al. (1999) JPO 29 pp 2719-2729.
c --- kappaf1 used internally to simplify offsetting T and S,
c --- always invoke via kappaf.
c --- offset limits based on stability estimates from:
c ---     Hallberg (2005) Ocean Modelling 8 pp 279-300.
c --- t: potential temperature (degC); s: salinity (psu);
c --- r: potential density (sigma); prs: pressure; kkf: ref.state
c ---     example: kappaf(4.5,34.5,36.406,1.e7,1) = -0.12301201
c ---     example: kappaf(4.5,34.5,36.406,1.e7,2) = -0.03356404
c ---     example: kappaf(4.5,34.5,36.406,1.e7,3) =  0.05201003
      kappaf1(t,s,r,prs,kkf)=(r+rhoref)*
     &  (exp(sclkap*(prs-pref)*
     &        ( s*( qs(kkf)+t* qst(kkf) ) +
     &          t*( qt(kkf)+t*(qtt(kkf)+t*qttt(kkf))+
     &              0.5*(prs+pref)*
     &              (qpt(kkf)+s*qpst(kkf)+t*qptt(kkf)) ) ) )
     &   -1.0)
      kappaf(t,s,r,prs,kkf)=
     &     kappaf1(max(-1.2,         t-toff(kkf) ),  !Hallberg,T-only: -1.8,0.9
     &             max(-3.0,min(1.5, s-soff(kkf))),  !Hallberg,S-only: -4.2,2.1
     &             r,prs,kkf)
c
c> Revision history
c>
c> May  2000 - conversion to SI units
c> Jul  2000 - removed rarely used functions, constants via parameter
c> Jan  2002 - removed geometery functions
c> Dec  2002 - new thermobaricity fit with toff=0.0,soff=34.0
c> Jun  2003 - removed sigma4
c> Jun  2003 - added locally referenced sigma
c> Sep  2004 - added kkf to kappaf, select one of three reference states
c> Aug  2006 - more restrictive kappaf1 offset limits
c> Sep  2006 - 9-term polynominal fit to T:[-2:30],S:[18:38]
c> May  2007 - added sigver
c> Mar  2009 - modified limits in kappaf
c> Mar  2009 - more accurate kappaf, with potential density
c-----------------------------------------------------------------------------

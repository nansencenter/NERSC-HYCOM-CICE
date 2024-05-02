c-----------------------------------------------------------------------------
      integer, parameter ::
     &  sigver=7  !12-term sigma-0
csig2&  sigver=8  !12-term sigma-2
c
      real    sig,dsigdt,dsigds,tofsig,sofsig,kappaf,
     &        sigloc,dsiglocdt,dsiglocds
c
      real    sig_n,sig_d,sig_q, dsigdt_n,dsigdt_d, dsigds_n,dsigds_d
      real    tofsig_a,tofsig_b,tofsig_c
      real    sofsig_a,sofsig_b,sofsig_c
      real    kappaf1
      real    sigloc_n,sigloc_d,sigloc_q,
     &        dsiglocdt_n,dsiglocdt_d, dsiglocds_n,dsiglocds_d
c
      real    r,s,t,pdb,prs
      integer kkf
c
      real, parameter ::
     &   aone =1.0,
     &   ahalf=1.0/2.0,
     &   a3rd =1.0/3.0, athird =a3rd,
     &   a4th =1.0/4.0, afourth=a4th
c
      real, parameter ::
     &   c1= 1.0e-01,      !not used, required to compile mxkrtm.f
     &   c2= 1.0e-02,      !not used, required to compile mxkrtm.f
     &   c3= 1.0e-03,      !not used, required to compile mxkrtm.f
     &   c4= 1.0e-04,      !not used, required to compile mxkrtm.f
     &   c5= 1.0e-05,      !not used, required to compile mxkrtm.f
     &   c6= 1.0e-06,      !not used, required to compile mxkrtm.f
     &   c7= 1.0e-07       !not used, required to compile mxkrtm.f
c
c --- REFERENCE?
c
c --- coefficients for 18-term rational function sigloc().
      real, parameter ::
     & c001=-1.4627567840659594d-01,  !num. constant    coefficent
     & c002= 6.4247392832635697d-02,  !num.    T        coefficent
     & c003= 8.1213979591704621d-01,  !num.       S     coefficent
     & c004=-8.1321489441909698d-03,  !num.    T^2      coefficent
     & c005= 4.5199845091090296d-03,  !num.    T  S     coefficent
     & c006= 4.6347888132781394d-04,  !num.       S^2   coefficent
     & c007= 5.0879498675039621d-03,  !num. P           coefficent
     & c008= 1.6333913018305079d-05,  !num. P  T        coefficent
     & c009= 4.3899924880543972d-06   !num. P     S     coefficent
      real, parameter ::
     & c011= 1.0000000000000000d+00,  !den. constant    coefficent
     & c012= 1.0316374535350838d-02,  !den.    T        coefficent
     & c013= 8.9521792365142522d-04,  !den.       S     coefficent
     & c014=-2.8438341552142710d-05,  !den.    T^2      coefficent
     & c015=-1.1887778959461776d-05,  !den.    T  S     coefficent
     & c016=-4.0163964812921489d-06,  !den.       S^2   coefficent
     & c017= 1.1995545126831476d-05,  !den. P           coefficent
     & c018= 5.5234008384648383d-08,  !den. P  T        coefficent
     & c019= 8.4310335919950873d-09   !den. P     S     coefficent
      real, parameter ::
     & c004x2=c004*2.d0,              !for dsigdt and dsiglocdt
     & c014x2=c014*2.d0,              !for dsigdt and dsiglocdt
     & c006x2=c006*2.d0,              !for dsigds and dsiglocds
     & c016x2=c016*2.d0               !for dsigds and dsiglocds
      real, parameter ::
     & sqrmin=0.d0,                   !sqrt arg can't be negative
     & sofmin=0.d0                    !salinity can't be negative
c --- reference pressure.
      real, parameter :: prs2pdb=1.d-4     !Pascals to dbar
      real, parameter :: pref=   0.d0      !ref. pressure in Pascals, sigma0
csig2 real, parameter :: pref=2000.d4      !ref. pressure in Pascals, sigma2
      real, parameter :: rpdb=pref*prs2pdb !ref. pressure in dbar
c --- coefficients for 12-term rational function sig() at rpdb.
      real, parameter ::
     & c101=c001+rpdb*c007,           !num. constant    coefficent
     & c102=c002+rpdb*c008,           !num.    T        coefficent
     & c103=c003+rpdb*c009            !num.       S     coefficent
      real, parameter ::
     & c111=c011+rpdb*c017,           !num. constant    coefficent
     & c112=c012+rpdb*c018,           !num.    T        coefficent
     & c113=c013+rpdb*c019            !num.       S     coefficent
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
c --- -----------------
c --- equation of state
c --- -----------------
c
c --- sigma at rpdb (dbar) as a function of temp (deg c) and salinity (psu)
c
      sig_n(t,s) = c101+(c102+c004*t+c005*s)*t  +
     &                  (c103+       c006*s)*s 
      sig_d(t,s) = c111+(c112+c014*t+c015*s)*t  +
     &                  (c113       +c016*s)*s 
      sig_q(t,s) = aone/sig_d(t,s)
      sig(  t,s) = sig_n(t,s)*sig_q(t,s)
c
c --- d(sig)/dt
      dsigdt_n(t,s) = c102+c004x2*t+c005*s
      dsigdt_d(t,s) = c112+c014x2*t+c015*s
      dsigdt(  t,s) = (dsigdt_n(t,s)-
     &                 dsigdt_d(t,s)*sig_n(t,s)*sig_q(t,s))*sig_q(t,s)
c
c --- d(sig)/ds
      dsigds_n(t,s) = c103+c005*t+c006x2*s
      dsigds_d(t,s) = c113+c015*t+c016x2*s
      dsigds(  t,s) = (dsigds_n(t,s)-
     &                 dsigds_d(t,s)*sig_n(t,s)*sig_q(t,s))*sig_q(t,s)
c
c --- temp (deg c) as a function of sigma and salinity (psu)
c --- find a quadratic polynominal root of a*t**2+b*t+c=0
      tofsig_a(r,s)=(   c004 -
     &               r* c014   )                  !quadratic coefficient
      tofsig_b(r,s)=(  (c102+      c005*s) -
     &               r*(c112+      c015*s)  )     !linear    coefficient
      tofsig_c(r,s)=(  (c101+(c103+c006*s)*s) -
     &               r*(c111+(c113+c016*s)*s)  )  !constant  coefficient
      tofsig(r,s)=( -tofsig_b(r,s)
     &              -sqrt(max(sqrmin,
     &                            tofsig_b(r,s)**2 -
     &                        4.0*tofsig_a(r,s)*tofsig_c(r,s))) ) /
     &            (2.0*tofsig_a(r,s))
c
c --- salinity (psu) as a function of sigma and temperature (deg c)
c --- find a quadratic polynominal root of a*s**2+b*s+c=0
      sofsig_a(r,t)=(   c006 -
     &               r* c016   )                  !quadratic coefficient
      sofsig_b(r,t)=(  (c103+      c005*t) -
     &               r*(c113+      c015*t)  )     !linear    coefficient
      sofsig_c(r,t)=(  (c101+(c102+c004*t)*t) -
     &               r*(c111+(c112+c014*t)*t)  )  !constant  coefficient
      sofsig(r,s)=max(sofmin,
     &                ( -sofsig_b(r,s)
     &                  +sqrt(max(sqrmin,
     &                                sofsig_b(r,s)**2 -
     &                            4.0*sofsig_a(r,s)*sofsig_c(r,s))) ) /
     &                (2.0*sofsig_a(r,s)) )
c
c --- locally referenced sigma, using the 18-term equation of state.
c --- t: potential temperature; s: psu; prs: pressure
c
      sigloc_n(t,s,pdb) = c001+(c002+c004*t+c005*s)*t  +
     &                         (c003+       c006*s)*s  +
     &                         (c007+c008*t+c009*s)*pdb
      sigloc_d(t,s,pdb) = c011+(c012+c014*t+c015*s)*t  +
     &                         (c013       +c016*s)*s  +
     &                         (c017+c018*t+c019*s)*pdb
      sigloc_q(t,s,pdb) = aone/sigloc_d(t,s,pdb)
      sigloc(  t,s,prs)=sigloc_n(t,s,prs*prs2pdb)*
     &                  sigloc_q(t,s,prs*prs2pdb)
c
c --- d(sig)/dt
      dsiglocdt_n(t,s,pdb) = c002+c004x2*t+c005*s+c008*pdb
      dsiglocdt_d(t,s,pdb) = c012+c014x2*t+c015*s+c018*pdb
      dsiglocdt(  t,s,prs)=(dsiglocdt_n(t,s,prs*prs2pdb)-
     &                      dsiglocdt_d(t,s,prs*prs2pdb)*
     &                         sigloc_n(t,s,prs*prs2pdb)*
     &                         sigloc_q(t,s,prs*prs2pdb) ) *
     &                         sigloc_q(t,s,prs*prs2pdb)
c
c --- d(sig)/ds
      dsiglocds_n(t,s,pdb) = c003+c005*t+c006x2*s+c009*pdb
      dsiglocds_d(t,s,pdb) = c013+c015*t+c016x2*s+c019*pdb
      dsiglocds(  t,s,prs)=(dsiglocds_n(t,s,prs*prs2pdb)-
     &                      dsiglocds_d(t,s,prs*prs2pdb)*
     &                         sigloc_n(t,s,prs*prs2pdb)*
     &                         sigloc_q(t,s,prs*prs2pdb) ) *
     &                         sigloc_q(t,s,prs*prs2pdb)
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
c> Oct  2010 - 17-term rational function equation of state
c> Oct  2010 - 12-term rational function equation of state
c-----------------------------------------------------------------------------

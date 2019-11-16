import numpy

# Sigma class. NB: This is the standard 7-term approach
class Sigma(object) :

   DTHIRD=1.0/3.0


   def __init__(self,flag) :
# --- coefficients for sigma-0 (based on Brydon & Sun fit)
      self._sigma=flag
      if flag == 0 :
         self.C1=-1.36471E-01
         self.C2= 4.68181E-02
         self.C3= 8.07004E-01
         self.C4=-7.45353E-03
         self.C5=-2.94418E-03
         self.C6= 3.43570E-05
         self.C7= 3.48658E-05
         self.pref=0.       # reference pressure, Pascals
      elif flag == 2 :
         self.C1= 9.77093E+00#  !const. coefficent
         self.C2=-2.26493E-02#  !T      coefficent
         self.C3= 7.89879E-01#  !   S   coefficent
         self.C4=-6.43205E-03#  !T^2    coefficent
         self.C5=-2.62983E-03#  !T  S   coefficent
         self.C6= 2.75835E-05#  !T^3    coefficent
         self.C7= 3.15235E-05#  !T^2S   coefficent
         self.pref=2000.0e4 # reference pressure, Pascals
      else :
         raise ValueError,"flag<>0 not implemented"

   def A0(self,S) :
      return (self.C1+self.C3*S)/self.C6

   def A1(self,S) :
      return (self.C2+self.C5*S)/self.C6

   def A2(self,S) :
      return (self.C4+self.C7*S)/self.C6

   def CUBQ(self,S) :
      return self.DTHIRD*A1(S)-(self.DTHIRD*A2(S))**2

   def CUBR(self,R,S) :
      return self.DTHIRD*(0.50*A1(S)*A2(S)-1.50*(A0(S)-R/self.C6)) -(self.DTHIRD*A2(S))**3

   def CUBAN(self,R,S) :
      return self.DTHIRD*ATAN2(SQRT(MAX(DZERO,-(self.CUBQ(S)**3+CUBR(R,S)**2))),CUBR(R,S))

   def CUBRL(self,R,S) :
      return SQRT(-self.CUBQ(S))*COS(self.CUBAN(R,S))

   def CUBIM(self,R,S) :
      return SQRT(-self.CUBQ(S))*SIN(self.CUBAN(R,S))

# --- temp (deg c) as a function of sigma and salinity (mil)
   def TOFSIG(self,R,S) :
      return -self.CUBRL(R,S)+SQRT(3.)*self.CUBIM(R,S)-self.DTHIRD*self.A2(S)

# --- salinity (mil) as a function of sigma and temperature (deg c)
   def SOFSIG(self,R,T) :
      return (R-self.C1-T*(self.C2+T*(self.C4+self.C6*T)))/(self.C3+T*(self.C5+self.C7*T))

# --- sigma-theta as a function of temp (deg c) and salinity (mil)
# --- (friedrich-levitus 3rd degree polynomial fit)
   def sig(self,T,S) :
      return (self.C1+self.C3*S+T*(self.C2+self.C5*S+T*(self.C4+self.C7*S+self.C6*T)))


# --- auxiliary statements for finding root of 3rd degree polynomial
#     A0(S)=(C1+C3*S)/C6
#     A1(S)=(C2+C5*S)/C6
#     A2(S)=(C4+C7*S)/C6
#     CUBQ(S)=DTHIRD*A1(S)-(DTHIRD*A2(S))**2
#     CUBR(R,S)=DTHIRD*(0.5D0*A1(S)*A2(S)-1.5D0*(A0(S)-R/C6))
#    .   -(DTHIRD*A2(S))**3
# --- if q**3+r**2>0, water is too dense to yield real root at given
# --- salinitiy. setting q**3+r**2=0 in that case is equivalent to
# --- lowering sigma until a double real root is obtained.
#     CUBAN(R,S)=DTHIRD*ATAN2(SQRT(MAX(DZERO,
#    .   -(CUBQ(S)**3+CUBR(R,S)**2))),CUBR(R,S))
#     CUBRL(R,S)=SQRT(-CUBQ(S))*COS(CUBAN(R,S))
#     CUBIM(R,S)=SQRT(-CUBQ(S))*SIN(CUBAN(R,S))
#
# --- temp (deg c) as a function of sigma and salinity (mil)
#     TOFSIG(R,S)=-CUBRL(R,S)+SQRT(3.)*CUBIM(R,S)-DTHIRD*A2(S)
#
# --- salinity (mil) as a function of sigma and temperature (deg c)
#     SOFSIG(R,T)=(R-C1-T*(C2+T*(C4+C6*T)))/(C3+T*(C5+C7*T))
#
# --- sigma-theta as a function of temp (deg c) and salinity (mil)
# --- (friedrich-levitus 3rd degree polynomial fit)
#     SIG(T,S)=(C1+C3*S+T*(C2+C5*S+T*(C4+C7*S+C6*T)))

   @property
   def sigma(self) : return self._sigma




class Sigma12Term(object) :
     _c001=-1.4627567840659594e-01 #num. constant    coefficent
     _c002= 6.4247392832635697e-02 #num.    T        coefficent
     _c003= 8.1213979591704621e-01 #num.       S     coefficent
     _c004=-8.1321489441909698e-03 #num.    T^2      coefficent
     _c005= 4.5199845091090296e-03 #num.    T  S     coefficent
     _c006= 4.6347888132781394e-04 #num.       S^2   coefficent
     _c007= 5.0879498675039621e-03 #num. P           coefficent
     _c008= 1.6333913018305079e-05 #num. P  T        coefficent
     _c009= 4.3899924880543972e-06 #num. P     S     coefficent
     #
     _c011= 1.0000000000000000e+00 #den. constant    coefficent
     _c012= 1.0316374535350838e-02 #den.    T        coefficent
     _c013= 8.9521792365142522e-04 #den.       S     coefficent
     _c014=-2.8438341552142710e-05 #den.    T^2      coefficent
     _c015=-1.1887778959461776e-05 #den.    T  S     coefficent
     _c016=-4.0163964812921489e-06 #den.       S^2   coefficent
     _c017= 1.1995545126831476e-05 #den. P           coefficent
     _c018= 5.5234008384648383e-08 #den. P  T        coefficent
     _c019= 8.4310335919950873e-09 #den. P     S     coefficent
     #
     _c004x2=_c004*2.e0 #for dsigdt and dsiglocdt
     _c014x2=_c014*2.e0 #for dsigdt and dsiglocdt
     _c006x2=_c006*2.e0 #for dsigds and dsiglocds
     _c016x2=_c016*2.e0 #for dsigds and dsiglocds
     #
     _sqrmin=0.e0 #sqrt arg can't be negative
     _sofmin=0.e0 #salinity can't be negative
     _prs2pdb=1.e-4 

     def __init__(self,thflag) :
        if thflag==2 :
           self._pref=2000.e4      #ref. pressure in Pascals, sigma2
        else :
           self._pref=0.           #ref. pressure in Pascals, sigma2
        self._rpdb =self._pref*self._prs2pdb

#c --- reference pressure.
#      real, parameter :: prs2pdb=1.d-4     !Pascals to dbar
#csig0 real, parameter :: pref=   0.d0      !ref. pressure in Pascals, sigma0
#      real, parameter :: pref=2000.d4      !ref. pressure in Pascals, sigma2
#      real, parameter :: rpdb=pref*prs2pdb !ref. pressure in dbar
#c --- coefficients for 12-term rational function sig() at rpdb.
        self._c101=self._c001+self._rpdb*self._c007 #num. constant    coefficent
        self._c102=self._c002+self._rpdb*self._c008 #num.    T        coefficent
        self._c103=self._c003+self._rpdb*self._c009 #num.       S     coefficent
        self._c111=self._c011+self._rpdb*self._c017 #num. constant    coefficent
        self._c112=self._c012+self._rpdb*self._c018 #num.    T        coefficent
        self._c113=self._c013+self._rpdb*self._c019 #num.       S     coefficent


#c
#c --- -----------------
#c --- equation of state
#c --- -----------------
#c
#c --- sigma at rpdb (dbar) as a function of temp (deg c) and salinity (psu)
#c
#    sig_n(t,s) = c101+(c102+c004*t+c005*s)*t  +
#    &                  (c103+       c006*s)*s 
     def sig_n(self,t,s) :
        return self._c101+(self._c102+self._c004*t+self._c005*s)*t  + (self._c103+       self._c006*s)*s 

#    sig_d(t,s) = c111+(c112+c014*t+c015*s)*t  +
#    &                  (c113       +c016*s)*s 
     def sig_d(self,t,s) :
        return self._c111+(self._c112+self._c014*t+self._c015*s)*t  + (self._c113       +self._c016*s)*s 

#    sig_q(t,s) = aone/sig_d(t,s)
     def sig_q(self,t,s) :
        return 1.0/self.sig_d(t,s)

#    sig(  t,s) = sig_n(t,s)*sig_q(t,s)
     def sig(self,t,s) :
        return self.sig_n(t,s)*self.sig_q(t,s)


#c
#c --- d(sig)/dt
#      dsigdt_n(t,s) = c102+c004x2*t+c005*s
#      dsigdt_d(t,s) = c112+c014x2*t+c015*s
#      dsigdt(  t,s) = (dsigdt_n(t,s)-
#     &                 dsigdt_d(t,s)*sig_n(t,s)*sig_q(t,s))*sig_q(t,s)
#c
#c --- d(sig)/ds
#      dsigds_n(t,s) = c103+c005*t+c006x2*s
#      dsigds_d(t,s) = c113+c015*t+c016x2*s
#      dsigds(  t,s) = (dsigds_n(t,s)-
#     &                 dsigds_d(t,s)*sig_n(t,s)*sig_q(t,s))*sig_q(t,s)
#c
#c --- temp (deg c) as a function of sigma and salinity (psu)
#c --- find a quadratic polynominal root of a*t**2+b*t+c=0
#      tofsig_a(r,s)=(   c004 -
#     &               r* c014   )                  !quadratic coefficient
#      tofsig_b(r,s)=(  (c102+      c005*s) -
#     &               r*(c112+      c015*s)  )     !linear    coefficient
#      tofsig_c(r,s)=(  (c101+(c103+c006*s)*s) -
#     &               r*(c111+(c113+c016*s)*s)  )  !constant  coefficient
#      tofsig(r,s)=( -tofsig_b(r,s)
#     &              -sqrt(max(sqrmin,
#     &                            tofsig_b(r,s)**2 -
#     &                        4.0*tofsig_a(r,s)*tofsig_c(r,s))) ) /
#     &            (2.0*tofsig_a(r,s))
#c
#c --- salinity (psu) as a function of sigma and temperature (deg c)
#c --- find a quadratic polynominal root of a*s**2+b*s+c=0
#      sofsig_a(r,t)=(   c006 -
#     &               r* c016   )                  !quadratic coefficient
#      sofsig_b(r,t)=(  (c103+      c005*t) -
#     &               r*(c113+      c015*t)  )     !linear    coefficient
#      sofsig_c(r,t)=(  (c101+(c102+c004*t)*t) -
#     &               r*(c111+(c112+c014*t)*t)  )  !constant  coefficient
#      sofsig(r,s)=max(sofmin,
#     &                ( -sofsig_b(r,s)
#     &                  +sqrt(max(sqrmin,
#     &                                sofsig_b(r,s)**2 -
#     &                            4.0*sofsig_a(r,s)*sofsig_c(r,s))) ) /
#     &                (2.0*sofsig_a(r,s)) )
#c
#c --- locally referenced sigma, using the 18-term equation of state.
#c --- t: potential temperature; s: psu; prs: pressure
#c
#      sigloc_n(t,s,pdb) = c001+(c002+c004*t+c005*s)*t  +
#     &                         (c003+       c006*s)*s  +
#     &                         (c007+c008*t+c009*s)*pdb
#      sigloc_d(t,s,pdb) = c011+(c012+c014*t+c015*s)*t  +
#     &                         (c013       +c016*s)*s  +
#     &                         (c017+c018*t+c019*s)*pdb
#      sigloc_q(t,s,pdb) = aone/sigloc_d(t,s,pdb)
#      sigloc(  t,s,prs)=sigloc_n(t,s,prs*prs2pdb)*
#     &                  sigloc_q(t,s,prs*prs2pdb)
#c
#c --- d(sig)/dt
#      dsiglocdt_n(t,s,pdb) = c002+c004x2*t+c005*s+c008*pdb
#      dsiglocdt_d(t,s,pdb) = c012+c014x2*t+c015*s+c018*pdb
#      dsiglocdt(  t,s,prs)=(dsiglocdt_n(t,s,prs*prs2pdb)-
#     &                      dsiglocdt_d(t,s,prs*prs2pdb)*
#     &                         sigloc_n(t,s,prs*prs2pdb)*
#     &                         sigloc_q(t,s,prs*prs2pdb) ) *
#     &                         sigloc_q(t,s,prs*prs2pdb)
#c
#c --- d(sig)/ds
#      dsiglocds_n(t,s,pdb) = c003+c005*t+c006x2*s+c009*pdb
#      dsiglocds_d(t,s,pdb) = c013+c015*t+c016x2*s+c019*pdb
#      dsiglocds(  t,s,prs)=(dsiglocds_n(t,s,prs*prs2pdb)-
#     &                      dsiglocds_d(t,s,prs*prs2pdb)*
#     &                         sigloc_n(t,s,prs*prs2pdb)*
#     &                         sigloc_q(t,s,prs*prs2pdb) ) *
#     &                         sigloc_q(t,s,prs*prs2pdb)




class Sigma17Term(object) :
#
# --- Jackett, McDougall, Feistel, Wright and Griffies (2006), 
# --- Algorithms for Density, Potential Temperature, Conservative
# --- Temperature, and the Freezing Temperature of Seawater, JAOT
#
# --- coefficients for 25-term rational function sigloc().
   _c001= 9.9984085444849347e+02 #num. constant    coefficent
   _c002= 7.3471625860981584e+00 #num.    T        coefficent
   _c003=-5.3211231792841769e-02 #num.    T^2      coefficent
   _c004= 3.6492439109814549e-04 #num.    T^3      coefficent
   _c005= 2.5880571023991390e+00 #num.       S     coefficent
   _c006= 6.7168282786692355e-03 #num.    T  S     coefficent
   _c007= 1.9203202055760151e-03 #num.       S^2   coefficent
   _c008= 1.0000000000000000e+00 #den. constant    coefficent
   _c009= 7.2815210113327091e-03 #den.    T        coefficent
   _c010=-4.4787265461983921e-05 #den.    T^2      coefficent
   _c011= 3.3851002965802430e-07 #den.    T^3      coefficent
   _c012= 1.3651202389758572e-10 #den.    T^4      coefficent
   _c013= 1.7632126669040377e-03 #den.       S     coefficent
   _c014= 8.8066583251206474e-06 #den.    T  S     coefficent
   _c015= 1.8832689434804897e-10 #den.    T^3S     coefficent
   _c016= 5.7463776745432097e-06 #den.    T  S^1.5 coefficent
   _c017= 1.4716275472242334e-09 #den.    T^3S^1.5 coefficent
   #
   _c018= 1.1798263740430364e-02 #num. P           coefficent
   _c019= 9.8920219266399117e-08 #num. P  T^2      coefficent
   _c020= 4.6996642771754730e-06 #num. P     S     coefficent
   _c021= 2.5862187075154352e-08 #num. P^2         coefficent
   _c022= 3.2921414007960662e-12 #num. P^2T^2      coefficent
   _c023= 6.7103246285651894e-06 #den. P           coefficent
   _c024= 2.4461698007024582e-17 #den. P^2T^3      coefficent
   _c025= 9.1534417604289062e-18 #den. P^3T        coefficent
   #c --- additional coefficients for dsiglocdt().
   _c031= 7.3471625860981580e+00 #num. constant    coefficent
   _c032=-1.0642246358568354e-01 #num.    T        coefficent
   _c033= 1.0947731732944364e-03 #num.    T^2      coefficent
   _c034= 6.7168282786692355e-03 #num.       S     coefficent
   _c035= 7.2815210113327090e-03 #den. constant    coefficent
   _c036=-8.9574530923967840e-05 #den.    T        coefficent
   _c037= 1.0155300889740728e-06 #den.    T^2      coefficent
   _c038= 5.4604809559034290e-10 #den.    T^3      coefficent
   _c039=-8.8066583251206470e-06 #den.       S     coefficent
   _c040= 5.6498068304414700e-10 #den.    T^2S     coefficent
   _c041= 2.9432550944484670e-09 #den.    T  S^1.5 coefficent
   _c042= 1.9784043853279823e-07 #num. P  T        coefficent
   _c043= 6.5842828015921320e-12 #num. P^2T        coefficent
   _c044= 7.3385094021073750e-17 #den. P^2T^2      coefficent
   _c045= 9.1534417604289060e-18 #den. P^3         coefficent
   #c --- additional coefficients for dsiglocds().
   _c051= 2.5880571023991390e+00 #num. constant    coefficent
   _c052= 6.7168282786692355e-03 #num.    T        coefficent
   _c053= 3.8406404111520300e-03 #num.       S     coefficent
   _c054= 1.7632126669040377e-03 #den. constant    coefficent
   _c055=-8.8066583251206470e-06 #den.    T        coefficent
   _c056= 1.8832689434804897e-10 #den.    T^3      coefficent
   _c057= 8.6195665118148150e-06 #den.       S^0.5 coefficent
   _c058= 2.2074413208363504e-09 #den.    T^2S^0.5 coefficent
   _c059= 4.6996642771754730e-06 #num. P           coefficent
   
   _sqrmin=0.e0  #sqrt arg can't be negative
   # --- reference pressure.
   _prs2pdb=1.e-4     #Pascals t=o dbar
   _rhoref=1.e3       #rhoref=qthref kg/m^3



   def __init__(self,thflag) :
      if thflag==2 :
         self._pref=2000.e4      #ref. pressure in Pascals, sigma2
      elif thflag==0 :
         self._pref=0.           #ref. pressure in Pascals, sigma2
      else :
         raise ValueError,"thflag must be 0 or 2"
      self._rpdb =self._pref*self._prs2pdb

      # --- coefficients for 17-term rational function sig() at rpdb.
      self._c101=self._c001+(self._c018-self._c021*self._rpdb)*self._rpdb  #num. constant    coefficent
      self._c103=self._c003+(self._c019-self._c022*self._rpdb)*self._rpdb  #num.    T^2      coefficent
      self._c105=self._c005+self._c020*self._rpdb                          #num.       S     coefficent
      self._c108=self._c008+self._c023*self._rpdb                          #den. constant    coefficent
      self._c109=self._c009-self._c025*self._rpdb**3                       #den.    T        coefficent
      self._c111=self._c011-self._c024*self._rpdb**2                       #den.    T^3      coefficent
      # --- additional coefficients for dsigdt().
      self._c132=self._c032+(self._c042-self._c043*self._rpdb)*self._rpdb, #num.    T        coefficent
      self._c135=self._c035-self._c045*self._rpdb**3                       #den. constant    coefficent
      self._c137=self._c037-self._c044*self._rpdb**2                       #den.    T^2      coefficent
      # --- additional coefficients for dsigds().
      self._c151=self._c051+self._c059*self._rpdb                          #num. constant    coefficent


# --- -----------------
# --- equation of state
# --- -----------------
#
# --- sigma at rpdb (dbar) as a function of pot.temp (deg c) and salinity (psu)
#
#   sig_n(t,s) = c101 + t*(c002+t*(c103+t*c004)) +
#  &                    s*(c105-t*c006+s*c007)
   def sig_n(self,t,s) :
      return self._c101 + \
         t*(self._c002+t*(self._c103+t*self._c004)) + \
         s*(self._c105-t*self._c006+s*self._c007)

#  sig_d(t,s) = c108 + t*(c109+t*(c010+t*(c111+t*c012))) +
#  &                    s*(c013-t*(c014+t*t*c015) +
#  &  sqrt(max(sqrmin,s))*(c016+t*t*c017))
   def sig_d(self,t,s) :
      return self._c108 +  \
         t*(self._c109+t*(self._c010+t*(self._c111+t*self._c012))) + \
         s*(self._c013-t*(self._c014+t*t*self._c015) + \
         numpy.sqrt(numpy.maximum(self._sqrmin,s))*(self._c016+t*t*self._c017))

#  sig_q(t,s) = aone/sig_d(t,s)
   def sig_q(self,t,s) :
      return 1.0/self.sig_d(t,s)

#  sig(  t,s) = sig_n(t,s)*sig_q(t,s) - rhoref
   def sig(self,t,s) :
      return self.sig_n(t,s)*self.sig_q(t,s) - self._rhoref

#c
#c --- d(sig)/dt
#      dsigdt_n(t,s) = c031 + t*(c132+t*c033) -
#     &                       s* c034
#      dsigdt_d(t,s) = c135 + t*(c036+t*(c137+t*c038)) +
#     &                       s*(c039-t*t*c040+
#     &             sqrt(max(sqrmin,s))*t*c041 )
#      dsigdt(  t,s) = (dsigdt_n(t,s)-
#     &                 dsigdt_d(t,s)*sig_n(t,s)*sig_q(t,s))*sig_q(t,s)
#c
#c --- d(sig)/ds
#c --- additional coefficients for dsigds().
#      dsigds_n(t,s) = c151 - t*c052 + s*c053
#      dsigds_d(t,s) = c054 +       t*(c055-t*t*c056) +
#     &           sqrt(max(sqrmin,s))*(c057+t*t*c058)
#      dsigds(  t,s) = (dsigds_n(t,s)-
#     &                 dsigds_d(t,s)*sig_n(t,s)*sig_q(t,s))*sig_q(t,s)
#c
#c --- temp (deg c) as a function of sigma and salinity (psu)
#c --- NOT AVAILABLE AS AN EXPRESSION - DO NOT USE
#      tofsig(r,s)=99.0
#c
#c --- salinity (psu) as a function of sigma and temperature (deg c)
#c --- NOT AVAILABLE AS AN EXPRESSION - DO NOT USE
#      sofsig(r,t)=99.0
#c
#c --- locally referenced sigma, using the 25-term equation of state.
#c --- t: potential temperature (degC); s: salinity (psu); prs: pressure (dbar)
#      sigloc_n(t,s,pdb) =      c001 +
#     &                      t*(c002 +
#     &                      t*(c003 +
#     &                      t* c004  )) +
#     &                      s*(c005 -
#     &                      t* c006 +
#     &                      s* c007  ) +
#     &                    pdb*(c018 +
#     &                    t*t* c019 +
#     &                      s* c020 -
#     &                    pdb*(c021 +
#     &                    t*t* c022  ))
#      sigloc_d(t,s,pdb) =      c008 +
#     &                      t*(c009 +
#     &                      t*(c010 +
#     &                      t*(c011 +
#     &                      t* c012  ))) +
#     &                      s*(c013 -
#     &                      t*(c014 +
#     &                    t*t* c015  ) +
#     &    sqrt(max(sqrmin,s))*(c016 +
#     &                    t*t* c017  )) +
#     &                    pdb*(c023 -
#     &             pdb*t*(t*t* c024 +
#     &                    pdb* c025  ))
#      sigloc_q(t,s,pdb) = aone/sigloc_d(t,s,pdb)
#      sigloc(t,s,prs)=sigloc_n(t,s,prs*prs2pdb)*
#     &                sigloc_q(t,s,prs*prs2pdb) - rhoref
#c
#c --- d(sig)/dt
#      dsiglocdt_n(t,s,pdb) =   c031 +
#     &                      t*(c032 +
#     &                      t* c033  ) -
#     &                      s* c034 +
#     &                  pdb*t*(c042 -
#     &                    pdb* c043  )
#      dsiglocdt_d(t,s,pdb) =   c035 +
#     &                      t*(c036 +
#     &                      t*(c037 +
#     &                      t* c038  )) +
#     &                      s*(c039 -
#     &                    t*t* c040 +
#     &  sqrt(max(sqrmin,s))*t* c041  ) -
#     &           pdb*pdb*(t*t* c044 +
#     &                    pdb* c045  )
#      dsiglocdt(t,s,prs)=(dsiglocdt_n(t,s,prs*prs2pdb)-
#     &                    dsiglocdt_d(t,s,prs*prs2pdb)*
#     &                       sigloc_n(t,s,prs*prs2pdb)*
#     &                       sigloc_q(t,s,prs*prs2pdb) ) *
#     &                       sigloc_q(t,s,prs*prs2pdb)
#c
#c --- d(sig)/ds
#      dsiglocds_n(t,s,pdb) =   c051 -
#     &                      t* c052 +
#     &                      s* c053 +
#     &                    pdb* c059
#      dsiglocds_d(t,s,pdb) =   c054 +
#     &                      t*(c055 -
#     &                    t*t* c056  ) +
#     &    sqrt(max(sqrmin,s))*(c057 +
#     &                    t*t* c058  )
#      dsiglocds(t,s,prs)=(dsiglocds_n(t,s,prs*prs2pdb)-
#     &                    dsiglocds_d(t,s,prs*prs2pdb)*
#     &                       sigloc_n(t,s,prs*prs2pdb)*
#     &                       sigloc_q(t,s,prs*prs2pdb) ) *
#     &                       sigloc_q(t,s,prs*prs2pdb)





class Kappa(object) :

# --- coefficients for kappa^(theta)
# --- new values (w.r.t. t-toff,s-soff,prs) from Shan Sun, Sep.2004
# --- 1=Arctic/Antarctic; 2=Atlantic; 3=Mediterranean
   _rhoref=1.e3  #rhoref=qthref kg/m^3
   _sclkap=1.e-11
   _toff = [  0.0,             3.0,            13.0 ]
   _soff = [ 34.5,            35.0,            38.5 ]
   _qttt = [ -3.03869354E-05, -3.03869352E-05, -3.03869353E-05 ]
   _qtt  = [  4.56625601E-03,  4.29277358E-03,  3.38116552E-03 ]
   _qt   = [ -2.88801209E-01, -2.61828868E-01, -1.81335007E-01 ]
   _qs   = [ -1.08670290E-01, -1.05131061E-01, -9.33336309E-02 ]
   _qst  = [  7.90503772E-04,  7.71096940E-04,  1.07270585E-03 ]
   _qpt  = [  1.07813750E-09,  1.00638435E-09,  7.57239852E-10 ]
   _qpst = [  1.41541548E-11,  1.48598578E-11,  3.89226107E-12 ]
   _qptt = [ -1.31383708E-11, -1.31383707E-11, -1.31383708E-11 ]


   def __init__(self,kkf,pref) :

      self._kkf =kkf-1 #  Start from zero
      self._pref=pref  #  Pa


# --- thermobaric compressibility coefficient (integral from prs to pref)
# ---     Sun et.al. (1999) JPO 29 pp 2719-2729.
# --- kappaf1 used internally to simplify offsetting T and S,
# --- always invoke via kappaf.
# --- offset limits based on stability estimates from:
# ---     Hallberg (2005) Ocean Modelling 8 pp 279-300.
# --- t: potential temperature (degC); s: salinity (psu);
# --- r: potential density (sigma); prs: pressure; kkf: ref.state
# ---     example: kappaf(4.5,34.5,36.406,1.e7,1) = -0.12301201 
# ---     example: kappaf(4.5,34.5,36.406,1.e7,2) = -0.03356404 
# ---     example: kappaf(4.5,34.5,36.406,1.e7,3) =  0.05201003
#     kappaf1(t,s,r,prs,kkf)=(r+rhoref)*
#    &  (exp(sclkap*(prs-pref)*
#    &        ( s*( qs(kkf)+t* qst(kkf) ) +
#    &          t*( qt(kkf)+t*(qtt(kkf)+t*qttt(kkf))+
#    &              0.5*(prs+pref)*
#    &              (qpt(kkf)+s*qpst(kkf)+t*qptt(kkf)) ) ) )
#    &   -1.0)
#     kappaf(t,s,r,prs,kkf)=
#    &     kappaf1(max(-1.2,         t-toff(kkf) ),  !Hallberg,T-only: -1.8,0.9
#    &             max(-3.0,min(1.5, s-soff(kkf))),  !Hallberg,S-only: -4.2,2.1
#    &             r,prs,kkf)
   def kappaf1(self,t,s,r,prs):
      tmp1= 0.5*(prs+self._pref) *(self._qpt[self._kkf]+s*self._qpst[self._kkf]+t*self._qptt[self._kkf])
      tmp2= s*( self._qs[self._kkf]+t* self._qst[self._kkf] ) + \
         t*( self._qt[self._kkf]+t*(self._qtt[self._kkf]+t*self._qttt[self._kkf]) + \
         tmp1)
      return (r+self._rhoref)* (numpy.exp( self._sclkap*(prs-self._pref)* (tmp2))-1.)
#      return (r+self._rhoref)* (exp( self._sclkap*(prs-self._pref)* ( 
#               s*( self._qs[self._kkf]+t* self._qst[self._kkf] ) + t*( self._qt[self._kkf]+t*(self.qtt[self._kkf]+t*self._qttt[self._kkf])+ 
#               0.5*(sprs+self._pref)* (self._qpt[self._kkf]+s*self._qpst[self._kkf]+t*self._qptt[self._kkf]) ) ) )
#         -1.0)

   #Hallberg,T-only: -1.8,0.9
   #Hallberg,S-only: -4.2,2.1
   def kappaf(self,t,s,r,prs):
      return self.kappaf1(numpy.maximum(-1.2,                   t-self._toff[self._kkf] ),\
                          numpy.maximum(-3.0,numpy.minimum(1.5, s-self._soff[self._kkf])),\
                          r,prs)


# For unittest of kappa (pref=0)
# ---     example: kappaf(4.5,34.5,36.406,1.e7,1) = -0.12301201 
# ---     example: kappaf(4.5,34.5,36.406,1.e7,2) = -0.03356404 
# ---     example: kappaf(4.5,34.5,36.406,1.e7,3) =  0.05201003

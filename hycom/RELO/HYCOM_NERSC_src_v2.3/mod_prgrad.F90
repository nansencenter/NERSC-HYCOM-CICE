      module mod_prgrad
      implicit none

      contains

#if defined(PRGRAD_FROM_GEOPOT_7TERM)
! ***********************************************************************************
! --- alpha = 1/rho = 1. / (c_0 + c_1*p + c_2*p**2) 
! ***********************************************************************************
!
! --- Convenience functions
      real function myc0(myt,mys) 
      implicit none
      include 'stmt_fns.h'
      real, intent(in) :: myt,mys
         myc0=1000. + alphap(1)+  &
                    alphap(2)*myt      +  &
                    alphap(3)*mys      +  &
                    alphap(4)*myt**2   +  &
                    alphap(5)*myt*mys    +  &
                    alphap(6)*myt**3   +  &
                    alphap(7)*myt**2*mys
      end function myc0
!
      real function  myc1(myt,mys)
      implicit none 
      include 'stmt_fns.h'
      real, intent(in) :: myt,mys
      myc1=        betap(1)        +   &
                 betap(2)*myt      +   &
                 betap(3)*mys      +   &
                 betap(4)*myt**2   +   &
                 betap(5)*myt*mys    +   &
                 betap(6)*myt**3   +   &
                 betap(7)*myt**2*mys
      myc1=myc1*prs2pb
      end function myc1
!
      real function  myc2(myt,mys)
      implicit none 
      include 'stmt_fns.h'
      real, intent(in) :: myt,mys
      myc2 =       gammap(1)        +  &
                 gammap(2)*myt      +  &
                 gammap(3)*mys      +  &
                 gammap(4)*myt**2   +  &
                 gammap(5)*myt*mys    +  &
                 gammap(6)*myt**3   +  &
                 gammap(7)*myt**2*mys
      myc2= myc2*prs2pb**2
      end function myc2
!
      real function d1(t,s)
      implicit none 
      real, intent(in) :: t,s
      d1 = myc1(t,s) / (2.*myc2(t,s))
      end function d1
!
      real function d2(t,s)
      implicit none
      real, intent(in) :: t,s
      d2= - (4*myc0(t,s)*myc2(t,s) - myc1(t,s)**2) / (4*myc2(t,s)**2) 
      d2= sqrt(d2)
      end function d2
!
      real function rholoc(myt,mys,myp)
      implicit none
      include 'stmt_fns.h'
      real, intent(in) :: myt,mys,myp
      rholoc=1000.+sigloc(myt,mys,myp)
      end function rholoc
!
! --- Integrated specific density from pb to pt. T,S is 
! --- Assumed constant over layer.
      real function alphint(t,s,phib,pb,pt) 
      implicit none
      real, intent(in) :: t,s,phib,pb,pt
      real :: ub,ut, tmp0,tmp1,tmp2
      ub=pb+d1(t,s)
      ut=pt+d1(t,s)
      tmp0 = 1./(myc2(t,s) * d2(t,s) *2.)
      tmp1 = log(d2(t,s)-ub) - log(d2(t,s)+ub) 
      tmp1 = tmp0*tmp1
      tmp2 = log(d2(t,s)-ut) - log(d2(t,s)+ut) 
      tmp2 = tmp0*tmp2
      alphint = phib + tmp2 - tmp1
      end function alphint
!
! --- Integrated specific density from pb to pt. T,S is 
! --- Assumed constant over layer. Linearized version
      real function alphintl(t,s,phib,pb,pt) 
      implicit none
      real, intent(in) :: t,s,phib,pb,pt
      real :: tmp0, t0, t1, t2, t3, t4, t5
      tmp0 = myc1(t,s)+2*myc2(t,s)*pb
      t0 = 0.
      t1 = phib
      t2 = 1./rholoc(t,s,pb)
      t3 = - tmp0/rholoc(t,s,pb)**2
      t4 = 2* tmp0**2 /rholoc(t,s,pb)**3
      t4 = t4 - 2*myc2(t,s)/rholoc(t,s,pb)**2
      t5 =    - 6 * tmp0**3 /rholoc(t,s,pb)**4 
      t5 = t5 + 8 * myc2(t,s)*tmp0/rholoc(t,s,pb)**3
      t5 = t5 + 4 * myc2(t,s)*tmp0/rholoc(t,s,pb)**3
      alphintl = t1 + t2*(pt-pb) + t3*(pt-pb)**2/2.  &
                    + t4*(pt-pb)**3/6. + t5*(pt-pb)**4/24.
      end function alphintl
!
! --- Integrated geopotential density from pb to pt. T,S is 
! --- Assumed constant over layer.
      real function intgp(t,s,phib,pb,pt) 
      implicit none
      real, intent(in) :: t,s,phib,pb,pt
      real :: ub,ut,tmp1,tmp2,tmp3
      ub=pb+d1(t,s)
      ut=pt+d1(t,s)
      tmp1 = - d2(t,s) / (2.*myc2(t,s)*d2(t,s)) *  &
         log(  (d2(t,s)+ut) * (d2(t,s)-ut)  &
             /  &
            ( (d2(t,s)+ub) * (d2(t,s)-ub) )  &
            )
      tmp2 = - ut  / (2.*myc2(t,s)*d2(t,s)) *  &
         log(  (d2(t,s)+ut) * (d2(t,s)-ub)  &
             /  &
             ( (d2(t,s)+ub) * (d2(t,s)-ut) )  &
         )
      tmp3 = phib*(ut-ub)
      intgp = tmp1 + tmp2 + tmp3
      end function intgp
!
! --- Integrated geopotential density from pb to pt. T,S is 
! --- Assumed constant over layer. Linearized version
      real function intgpl(t,s,phib,pb,pt)
      implicit none
      real, intent(in) :: t,s,phib,pb,pt
      real :: tmp0, t0, t1, t2, t3, t4, t5
      tmp0 = myc1(t,s)+2*myc2(t,s)*pb
      t0 = 0.
      t1 = phib
      t2 = 1./rholoc(t,s,pb)
      t3 = - tmp0/rholoc(t,s,pb)**2
      t4 = 2* tmp0**2 /rholoc(t,s,pb)**3
      t4 = t4 - 2*myc2(t,s)/rholoc(t,s,pb)**2
      t5 =    - 6 * tmp0**3 /rholoc(t,s,pb)**4 
      t5 = t5 + 8 * myc2(t,s)*tmp0/rholoc(t,s,pb)**3
      t5 = t5 + 4 * myc2(t,s)*tmp0/rholoc(t,s,pb)**3
      intgpl = t0 + t1*(pt-pb) + t2*(pt-pb)**2/2.  &
         + t3*(pt-pb)**3/6. + t4*(pt-pb)**4/24. + t5*(pt-pb)**5/120.
      end function intgpl
!
! --- Get pressure at bottom from geopot and pt. Can be used to construct
! --- pressure at interfaces from geopotentials
      real function  pb_from_geop(t,s,phib,phit,pt)
      implicit none
      real, intent(in) :: t,s,phib,phit,pt
      real  :: f1, f2, ut, ub
      ut=pt+d1(t,s)
      f1=(phit-phib)*2*myc2(t,s)*d2(t,s)
      f1=exp(f1)
      f2=(d2(t,s)+ut)/(d2(t,s)-ut)
      ub=(f1*f2-1)*d2(t,s)/(1+f1*f2)
      pb_from_geop=ub-d1(t,s)
      end function pb_from_geop
!
#elif defined(PRGRAD_FROM_GEOPOT_12TERM)
! ***********************************************************************************
! --- alpha = 1/rho = ( n_0 + n_1*p) / (d_0 + d_1*p) 
! ***********************************************************************************
!
      real function n_0(t,s) 
      implicit none
      real, intent(in) :: t,s
      n_0 = c011+(c012+c014*t+c015*s)*t  +  (c013       +c016*s)*s  
      end function n_0
!
      real function n_1(t,s) 
      implicit none
      real, intent(in) :: t,s
      n_1 = (c017+c018*t+c019*s)*prs2pdb
      end function n_1
!
      real function d_0(t,s)
      implicit none
      real, intent(in) :: t,s
      d_0 = 1000.*n_0(t,s) + c001 + (c002+c004*t+c005*s)*t + &
            (c003+c006*s)*s
      end function d_0
!
      real function d_1(t,s)
      implicit none
      real, intent(in) :: t,s
      d_1 = 1000.*n_1(t,s) + (c007+c008*t+c009*s)*prs2pdb
      end function
!
      real function rholoc(t,s,p)
      implicit none
      real, intent(in) :: t,s,p
      rholoc=1000.+sigloc(t,s,p)
      end function rholoc
!
! --- Integrated specific density from pb to pt. T,S is 
! --- Assumed constant over layer.
      real function alphint(t,s,phib,pb,pt) 
      implicit none
      real, intent(in) :: t,s,phib,pt,pb
      real :: fac
      fac = (n_0(t,s)*d_1(t,s) - n_1(t,s)*d_0(t,s)) / d_1(t,s)**2
      t1 = fac * log( (d_0(t,s)+d_1(t,s)*pt)/(d_0(t,s)+d_1(t,s)*pb))
      t2 = n_1(t,s) * (pt-pb) / d_1(t,s) + phib
      alphint =  t1+t2
      end function alphint
!
! --- Integrated specific density from pb to pt. T,S is 
! --- Assumed constant over layer. Linearized version
      real function alphintl(t,s,phib,pb,pt) 
      implicit none
      real, intent(in) :: t,s,phib,pt,pb
      real :: fac, t0, t1, t2, t3, t4, t5
      t0  =  phib 
      t1  =  1./rholoc(t,s,pb)
      t2  =  ( n_1(t,s)*d_0(t,s) - n_0(t,s)*d_1(t,s)  )  &
         / (d_0(t,s) + d_1(t,s)*pb)**2
      t3  = -2. * d_1(t,s) * t2 /  (d_0(t,s) + d_1(t,s)*pb)
      t4  = -3. * d_1(t,s) * t3 /  (d_0(t,s) + d_1(t,s)*pb)
      t5  = -4. * d_1(t,s) * t4 /  (d_0(t,s) + d_1(t,s)*pb)
      alphintl=t0 + t1*(pt-pb) + t2*(pt-pb)**2/2. + t3*(pt-pb)**3/6.  &
         + t4*(pt-pb)**4/24. + t5*(pt-pb)**5/120.
      end function alphintl
!
! --- Helper function for intgp
      real function intgp_indefinite(t,s,phib,pb,p) 
      implicit none
      real, intent(in) :: t,s,phib,pt,pb
      real :: fac, t0, t1, t2, t3, t4, t5
!
      fac = (n_0(t,s)*d_1(t,s) - n_1(t,s)*d_0(t,s)) / d_1(t,s)**2
      t1 = fac * ((d_0(t,s)+d_1(t,s)*p) * &
         log(d_0(t,s)+d_1(t,s)*p) - d_1(t,s)*p) / d_1(t,s)
      t2 = n_1(t,s)*p**2/(2*d_1(t,s)) 
      t3 = phib* p
      t4 = - fac * log (d_0(t,s)+d_1(t,s)*pb) * p
      t5 = - n_1(t,s) * pb * p / d_1(t,s)
      intgp_indefinite = t1+t2+t3+t4+t5
      end function intgp_indefinite
!
! --- Integrated geopotential density from pb to pt. T,S is 
! --- Assumed constant over layer. 
      real function intgp(t,s,phib,pb,p) 
      implicit none
      real, intent(in) :: t,s,phib,pt,pb
      intgp = intgp_indf(t,s,phib,pb,pt) - intgp_indf(t,s,phib,pb,pb) 
      end function intgp
!
! --- Integrated geopotential density from pb to pt. T,S is 
! --- Assumed constant over layer. Linearized version
      real function intgpl(t,s,phib,pb,p) 
      implicit none
      real, intent(in) :: t,s,phib,pt,pb
      real :: fac, t0, t1, t2, t3, t4, t5
!
! --- Taylor series based on alphint expression. Below are the derivatives of alphint at pt=pb
      t0  =  phib 
      t1  =  1./rholoc(t,s,pb)
      t2  =  ( n_1(t,s)*d_0(t,s) - n_0(t,s)*d_1(t,s)  )  &
         / (d_0(t,s) + d_1(t,s)*pt)**2
      t3  = -2 * d_1(t,s) * t2 /  (d_0(t,s) + d_1(t,s)*pb)
      t4  = -3 * d_1(t,s) * t3 /  (d_0(t,s) + d_1(t,s)*pb)
      t5  = -4 * d_1(t,s) * t4 /  (d_0(t,s) + d_1(t,s)*pb)
      intgpl = 0. + t0*(pt-pb) + t1*(pt-pb)**2/2. + t2*(pt-pb)**3/6.  &
         + t3*(pt-pb)**4/24. + t4*(pt-pb)**5/120. + t5*(pt-pb)**6/720.
      end function intgpl
!
#elif defined(PRGRAD_FROM_GEOPOT_17TERM)
! ***********************************************************************************
! --- alpha = 1/rho = ( n_0 + n_1*p + n_2*p^2 + n_3*p^3) / ( d_0 + d_1*p + d_2*p^2) 
! ***********************************************************************************
! --- Convenience functions
      real function n_0(t,s) 
      implicit none
      real, intent(in) :: t,s
      n_0 =        c008 +  &
               t*(c009 + &
               t*(c010 + &
               t*(c011 + &
               t* c012  ))) + &
               s*(c013 - &
               t*(c014 + &
             t*t* c015  ) + &
          sqrt(max(0.,s))*(c016 + &
             t*t* c017  )) 
      end function n_0
!
      real function n_1(t,s) 
      implicit none
      real, intent(in) :: t,s
      n_1 = c023*prs2pdb
      end function n_2
!
      real function n_2(t,s) 
      implicit none
      real, intent(in) :: t,s
      n_2 = - t*t*t*c024*prs2pdb**2
      end function
!
      real function n_3(t,s) 
      implicit none
      real, intent(in) :: t,s
      n_3 = -  t*c025*prs2pdb**3
      end function n_3
!
      real function d_0(t,s) 
      implicit none
      real, intent(in) :: t,s
      d_0  =         c001 + &
             t*(c002 + &
             t*(c003 + &
             t* c004  )) + &
             s*(c005 - &
             t* c006 + &
             s* c007  )
      end function d_0

      real function d_1(t,s) 
      implicit none
      real, intent(in) :: t,s
      d_1 =  (c018 + t*t* c019 + s* c020)*prs2pdb
      end function

      real function d_2(t,s) 
      implicit none
      real, intent(in) :: t,s
      d_2 = -(c021 +t*t*c022)*prs2pdb**2
      end function
!
      real function d_discriminant(t,s) 
      implicit none
      real, intent(in) :: t,s
      d_discriminant=d_1(t,s)**2  - 4*d_0(t,s)*d_2(t,s)
      end function d_discriminant

      function d_roots(t,s) 
      implicit none
      real, intent(in) :: t,s
      real d_roots(2), rootsq
      rootsq =  d_1(t,s)**2  - 4*d_0(t,s)*d_2(t,s)
      d_roots(1)   = (-d_1(t,s) + sqrt(rootsq)) /(2*d_2(t,s))
      d_roots(2)   = (-d_1(t,s) - sqrt(rootsq)) /(2*d_2(t,s))
      end function d_roots
!
      real function  alphastar(t,s) 
      implicit none
      real, intent(in) :: t,s
      alphastar = n_3(t,s)/d_2(t,s)
      end function alphastar
!
      real function betastar(t,s) 
      implicit none
      real, intent(in) :: t,s
      betastar = (n_2(t,s) - n_3(t,s)*d_1(t,s) / d_2(t,s)) / d_2(t,s)
      end function betastar
!
      real function mstar(t,s) 
      implicit none
      real, intent(in) :: t,s
      mstar = n_1(t,s) - n_3(t,s)*d_0(t,s) / d_2(t,s)  &
         - ( n_2(t,s) - n_3(t,s)*d_1(t,s)/d_2(t,s)) * d_1(t,s)/d_2(t,s)
      end function mstar
!
      real function nstar(t,s) 
      implicit none
      real, intent(in) :: t,s
      nstar = n_0(t,s)  &
         - (n_2(t,s) - n_3(t,s)*d_1(t,s)/d_2(t,s)) * d_0(t,s)/d_2(t,s)
      end function nstar
!
      real function astar(t,s) 
      implicit none
      real, intent(in) :: t,s
      astar = d_2(t,s)
      end function astar
!
      real function bstar(t,s) 
      implicit none
      real, intent(in) :: t,s
      bstar = d_1(t,s)
      end function bstar
!
      real function cstar(t,s) 
      implicit none
      real, intent(in) :: t,s
      cstar = d_0(t,s)
      end function cstar

      real function dstar(t,s)
      implicit none
      real, intent(in) :: t,s
      real tmp1, tmp2
      tmp1 = sqrt(d_discriminant(t,s))
      tmp2 = 2*astar(t,s)
      dstar = tmp1/tmp2
      end function dstar
!
      real function estar(t,s)
      implicit none
      real, intent(in) :: t,s
      real :: tmp1, tmp2
      tmp1 = 2*astar(t,s)*nstar(t,s) - bstar(t,s)*mstar(t,s)
      tmp2 = astar(t,s)*numpy.sqrt(d_discriminant(t,s))
      estar = tmp1/tmp2
      end function estar
!
! --- Helper function for exact alphint 
      real function alphint_indf(t,s,p) :
      implicit none
      real, intent(in) :: t,s,p
      real :: t1,t2,t3,t4
!
      t1 = 0.5 * alphastar(t,s) * p**2 
      t2 = betastar(t,s) * p
      t3 = mstar(t,s) * log(astar(t,s)*p**2 + bstar(t,s)*p+cstar(t,s)) &
         / (2*astar(t,s))
      t4n= 2*astar(t,s)*p + bstar(t,s)
      t4d= sqrt(d_discriminant(t,s))
      t4=t4n/t4d
      t4 =-estar(t,s)*numpy.arctanh(t4)
      alphint_indf = t1 + t2 + t3 + t4
      end function alphint_indf
!
! --- Integrated specific density from pb to pt. T,S is 
! --- Assumed constant over layer.
      real function alphint(t,s,phib,pb,pt) 
      implicit none
      real, intent(in) :: t,s,phib,pt,pb
      alphint = alphint_indf(t,s,pt) - alphint_indf(t,s,pb) +phib
      end function alphint
!
! --- Integrated specific density from pb to pt. T,S is 
! --- Assumed constant over layer. Linearized version
      real function alphintl(t,s,phib,pb,pt)
      implicit none
      real, intent(in) :: t,s,phib,pt,pb
      real t0, t1, t2, t3, t4, t5
      ! TODO taylor series expansion for alphint. Not enough terms atm  ...
      t0  =  phib 
      t1  =  1./rholoc(t,s,pb)
      t2  =  ( n_1(t,s) + 2*n_2(t,s)*pb + 3.*n_3(t,s)*pb**2)  &
         *(d_0(t,s)+d_1(t,s)*pb+d_2(t,s)*pb)
      t2  = t2 - (d_1(t,s)+2*d_2(t,s)*pb)  &
         * ( n_0(t,s) + n_1(t,s)*pb + n_2(t,s)*pb**2 + n_3(t,s)*pb**3)
      t2  = t2 / (d_0(t,s) + d_1(t,s)*pb + d_2(t,s)*pb**2)**2
      !t3  =  ?
      !t4  =  ?
      !t5  =  ?
      alphintl = t0 + t1*(pt-pb) + t2*(pt-pb)**2/2. !+ t3*(pt-pb)**3/6. + t4*(pt-pb)**4/24. + t5*(pt-pb)**5/120.
      end function alphintl
!
! --- Helper function for intgp
      real function int_C(t,s,pb,phib) :
      implicit none
      real, intent(in) :: t,s,phib,pb
      int_C  = - alphint_indf(t,s,pb) +phib
      end function

      real function u_from_p(t,s,p) 
      implicit none
      real, intent(in) :: t,s,p
      u_from_p  = ( p + bstar(t,s)/(2.*astar(t,s)) ) / dstar(t,s)
      end function u_from_p
!
! --- Helper function for intgp
      real function intgp_indf(t,s,p,c) :
      implicit none
      real, intent(in) :: t,s,phib,pb
      real :: u,p0, p1,tmpv(2), t1, t2, t3, t4, t5, t6
      u = u_from_p(t,s,p)
      tmpv  = d_roots(t,s)
      p0=tmpv(1)
      p1=tmpv(2)
      ! KAL - double check !
      !f0=numpy.sign(p-p0)
      !f1=numpy.sign(p-p1)
      f0=sign(1.,p-p0)
      f1=sign(1.,p-p1)
      t1=alphastar(t,s) * p**3 / 6.
      t2=betastar(t,s)  * p**2 / 2.
      t3=C * p
      t4=mstar(t,s) * ( &
           p * log(abs(astar(t,s))) &
         + (abs(p-p0)*log(abs(p-p0))-f0*p)/f0 &
         + (abs(p-p1)*log(abs(p-p1))-f1*p)/f1 &
         ) / (2*astar(t,s))
      t5= - estar(t,s) * dstar(t,s) * u * atanh(u) 
      t6= - estar(t,s) * dstar(t,s) * 0.5*log(1.-u**2)
      intgp_indf = t1+t2+t3+t4+t5+t6
      end function intgp_indf
!
! --- Integrated geopotential density from pb to pt. T,S is 
! --- Assumed constant over layer. 
      real function intgp(t,s,phib,pb,p) 
      implicit none
      real, intent(in) :: t,s,phib,pt,pb
      real :: c
      c=int_c(t,s,pb,phib)
      intgp = intgp_indf(t,s,pt,c)  - self.intgp_indf(t,s,pb,c)
      end function intgp
!
! --- Integrated geopotential density from pb to pt. T,S is 
! --- Assumed constant over layer. Linearized version
      real function intgpl(t,s,phib,pb,p) 
      implicit none
      real, intent(in) :: t,s,phib,pt,pb
      real :: fac, t0, t1, t2, t3, t4, t5
!
      ! Integration of taylor series expansion for alphint. TODO: Not enough terms at the moment
      t0  =  phib 
      t1  =  1./rholoc(t,s,pb)
      t2  =  ( n_1(t,s) + 2*n_2(t,s)*pb + 3.*n_3(t,s)*pb**2)  &
         *(d_0(t,s)+d_1(t,s)*pb+d_2(t,s)*pb)
      t2  = t2 - (d_1(t,s)+2*d_2(t,s)*pb)  &
         * ( n_0(t,s) + n_1(t,s)*pb + n_2(t,s)*pb**2 + n_3(t,s)*pb**3)
      t2  = t2 / (d_0(t,s) + d_1(t,s)*pb + d_2(t,s)*pb**2)**2
      !t3  =  ?
      !t4  =  ?
      !t5  =  ?
      intgpl = 0+ t0*(pt-pb) + t1*(pt-pb)**2/2. + t2*(pt-pb)**3/6. ! ...
      end function intgpl
#endif /*defined(PRGRAD_FROM_GEOPOT_XXTERM)*/


#if defined(PRGRAD_FROM_GEOPOT_7TERM) && defined(PRGRAD_FROM_GEOPOT_12TERM) && defined(PRGRAD_FROM_GEOPOT_17TERM)
!
! --- Calculate geopotential at top of layer from t, s
! --- assumed constant over layer. Exact version
      real function phitop(t,s,phib,pb,pt)
      implicit none
      real, intent(in) ::  t,s,phib,pb,pt   ! Temperature, salinity and pressure in bottom 
      phitop = alphint(t,s,phib,pb,pt) 
      end function 
!
! --- Calculate geopotential at top of layer from t, s
! --- assumed constant over layer. Linearized version
      real function phitopl(t,s,phib,pb,pt)
      implicit none
      real, intent(in) ::  t,s,phib,pb,pt   ! Temperature, salinity and pressure in bottom 
      phitopl = alphintl(t,s,phib,pb,pt) 
      end function 
!
! --- Calculate geopotential at top of layer from t, s
! --- assumed constant over staircase profile. "Exact" version
      real function phitop_sc(tb,sb,pb,tt,st,pt,phib,nump)
      implicit none
      real, intent(in) ::  tb,sb,pb      ! t, s and p at bottom of cell
      real, intent(in) ::  tt,st,pt      ! t, s and p at top of cell
      real, intent(in) ::  phib          ! phi at bottom
      integer, intent(in) :: nump        ! number of points in staircase profile
      real :: t,s,w,p, pb0
      integer :: k
!
      phitop_sc=phib
      pb0      =pb   ! pressure at bottom of staircase
      do k=1,nump 
         w   = float(k)/nump
         s   = w*st   + (1.-w)*sb
         t   = w*tt   + (1.-w)*tb
         p   = w*pt   + (1.-w)*pb
         phitop_sc = alphint(t,s,phitop_sc,pb0,p) 
         pb0=p
      end do
      end function  phitop_sc
!
! --- Calculate geopotential at top of layer from t, s
! --- assumed constant over staircase profile. "Linearized" version
      real function phitop_scl(tb,sb,pb,tt,st,pt,phib,nump)
      implicit none
      real, intent(in) ::  tb,sb,pb      ! t, s and p at bottom of cell
      real, intent(in) ::  tt,st,pt      ! t, s and p at top of cell
      real, intent(in) ::  phib          ! phi at bottom
      integer, intent(in) :: nump        ! number of points in staircase profile
      real :: t,s,w,p, pb0
      integer k
!
      phitop_scl=phib
      pb0      =pb   ! pressure at bottom of staircase
      do k=1,nump 
         w   = float(k)/nump
         s   = w*st   + (1.-w)*sb
         t   = w*tt   + (1.-w)*tb
         p   = w*pt   + (1.-w)*pb
         phitop_scl = alphint(t,s,phitop_scl,pb0,p) 
         pb0=p
      end do
      end function  phitop_scl
!
! --- Calculate geopotential at top of layer from t, s
! --- assumed constant over layer. "Exact" version 
      function phitopedge(tbl,sbl,pbl,tbr,sbr,pbr,ptl,ptr,phib,nump)
      implicit none
      integer, intent(in) :: nump
      real, intent(in) ::  tbl,sbl,pbl          ! t,s,p at bottom left
      real, intent(in) ::  tbr,sbr,pbr          ! t,s,p at bottom right
      real, intent(in) ::  ptl,ptr              ! p at top left and top right
      real, intent(in), dimension(nump) :: phib ! geopotential along bottom edge
      real,             dimension(nump) :: phitopedge
      real :: w,s,t,pb,pt
      integer :: i
      ! Linear interpolation of s, t from pb to pt
      do i=1,nump 
         w=float(i-1)/(nump-1)
         s = w*sbr + (1.-w)*sbl
         t = w*tbr + (1.-w)*tbl
         pb= w*pbr + (1.-w)*pbl
         pt= w*ptr + (1.-w)*ptl
         phitopedge(i) = phitop(t,s,phib(i),pb,pt) 
      end do
      end function phitopedge
!
! --- Calculate geopotential at top of layer from t, s
! --- assumed constant over layer. "Linearized" version 
      function phitopedgel(tbl,sbl,pbl,tbr,sbr,pbr,ptl,ptr,phib,nump)
      implicit none
      integer, intent(in) :: nump
      real, intent(in) ::  tbl,sbl,pbl          ! t,s,p at bottom left
      real, intent(in) ::  tbr,sbr,pbr          ! t,s,p at bottom right
      real, intent(in) ::  ptl,ptr              ! p at top left and top right
      real, intent(in), dimension(nump) :: phib ! geopotential along bottom edge
      real,             dimension(nump) :: phitopedgel
      real :: w,s,t,pb,pt
      integer :: i
      ! Linear interpolation of s, t from pb to pt
      do i=1,nump 
         w=float(i-1)/(nump-1)
         s = w*sbr + (1.-w)*sbl
         t = w*tbr + (1.-w)*tbl
         pb= w*pbr + (1.-w)*pbl
         pt= w*ptr + (1.-w)*ptl
         phitopedgel(i) = phitopl(t,s,phib(i),pb,pt) 
      end do
      end function phitopedgel
!
! --- Integration rules, Boole's or Simpson's rules
      real function intgp_numq(pl,pr,phi,nump)
      implicit none
      integer, intent(in) :: nump
      real,    intent(in) :: pl,pr        ! Pressure at left and right sides (bottom/top?)
      real,    intent(in) :: phi(nump)    ! Phi at nump points along edge
!
! --- Simpsons rule
      if (nump.eq.3) then
         intgp_numq = (pr-pl) * ( phi(1) + 4*phi(2) + phi(3)) / 6.
      elseif (nump.eq.5) then
         intgp_numq = (pr-pl) * ( 7*phi(1) + 32*phi(2) + 12*phi(3)  &
            + 32*phi(4) + 7*phi(5)) / 90.
      else
         ! brutalt !
         stop ('intgp_numq:nump not 3 or 5')
      end if
      end function intgp_numq
#endif /* defined(PRGRAD_FROM_GEOPOT_7TERM) && defined(PRGRAD_FROM_GEOPOT_12TERM) && defined(PRGRAD_FROM_GEOPOT_17TERM)*/
!
      end module mod_prgrad

!undef PARTICLE   # Include code for particle tracers
!undef ICE        # Save and read ice restart files
!undef HIBLER     # Include ice dynamics from Hardner
!undef ENSEMBLE   # Include code for ensemble stuff
!undef RIVER      # Include River fluxes (edit the ./Data/river.dat file)
!undef MIXDIAG    # computes and accumulates vertical transports
!undef NEST_OUTER # Save nesting conditions to be used in local models
!undef NEST_INNER # read nesting conditions saved by a global model
!undef OPEN_MP    # Parallelization using OPEN_MP directives.
!undef CLIMPRT    # Dont add just print current field and stop
!undef DEBUG_FRC
!undef DEBUG_RELAX
!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#define 
!#define ICESTATE_TEST
module m_insitu2theta
contains
   real function atg(s,t,p)
   ! ****************************
   ! adiabatic temperature gradient deg c per decibar
   ! ref: bryden,h.,1973,deep-sea res.,20,401-408
   ! units:
   !       pressure        p        decibars
   !       temperature     t        deg celsius (ipts-68)
   !       salinity        s        (ipss-78)
   !       adiabatic       atg      deg. c/decibar
   ! checkvalue: atg=3.255976e-4 c/dbar for s=40 (ipss-78),
   ! t=40 deg c,p0=10000 decibars
   implicit none
   real, intent(in) :: s,t,p
   real :: ds
   ds = s - 35.0
   atg = (((-2.1687e-16*t+1.8676e-14)*t-4.6206e-13)*p     &
       +((2.7759e-12*t-1.1351e-10)*ds+((-5.4481e-14*t     &
       +8.733e-12)*t-6.7795e-10)*t+1.8741e-8))*p          &
       +(-4.2393e-8*t+1.8932e-6)*ds                       &
       +((6.6228e-10*t-6.836e-8)*t+8.5258e-6)*t+3.5803e-5
   return
   end function atg
   real function insitu2theta(s,t0,p0,pr)
   ! ***********************************
   ! to compute local potential temperature at pr
   ! using bryden 1973 polynomial for adiabatic lapse rate
   ! and runge-kutta 4-th order integration algorithm.
   ! ref: bryden,h.,1973,deep-sea res.,20,401-408
   ! fofonoff,n.,1977,deep-sea res.,24,489-491
   ! units:
   !       pressure        p0       decibars
   !       temperature     t0       deg celsius (ipts-68)
   !       salinity        s        (ipss-78)
   !       reference prs   pr       decibars
   !       potential tmp.  theta    deg celsius
   ! checkvalue: theta= 36.89073 c,s=40 (ipss-78),t0=40 deg c,
   ! p0=10000 decibars,pr=0 decibars
   !
   !      set-up intermediate temperature and pressure variables
   implicit none
      
   real, intent(in) :: s, t0, p0, pr
   real :: p,t,h,xk,q
   p=p0
   t=t0
   !**************
   h = pr - p
   xk = h*atg(s,t,p)
   t = t + 0.5*xk
   q = xk
   p = p + 0.5*h
   xk = h*atg(s,t,p)
   t = t + 0.29289322*(xk-q)
   q = 0.58578644*xk + 0.121320344*q
   xk = h*atg(s,t,p)
   t = t + 1.707106781*(xk-q)
   q = 3.414213562*xk - 4.121320344*q
   p = p + 0.5*h
   xk = h*atg(s,t,p)
   insitu2theta = t + (xk-2.0*q)/6.0
   return
   end function insitu2theta
end module m_insitu2theta

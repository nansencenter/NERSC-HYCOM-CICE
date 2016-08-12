! CVS: $Id: mod_evp.F,v 1.1 2007/01/08 08:04:43 knutali Exp $
! CVS: $Source: /home/fimm/nersc/knutali/cvsroot/HYCOM_2.1.03_MPI/EVP_MPI/mod_evp.F,v $
! CVS: $Name:  $
!=======================================================================
!BOP
!
! !MODULE: ice_dyn_evp - elastic-viscous-plastic sea ice dynamics model 
!
! !DESCRIPTION:
!
! Elastic-viscous-plastic sea ice dynamics model \\
! Computes ice velocity and deformation \\
!
! See: \\
!
! Hunke, E. C., and J. K. Dukowicz (1997). An elastic-viscous-plastic model
! for sea ice dynamics. {\em J. Phys. Oceanogr.}, {\bf 27}, 1849--1867. \\
!
! Hunke, E. C. (2001).  Viscous-Plastic Sea Ice Dynamics with the EVP Model:
! Linearization Issues. {\em Journal of Computational Physics}, {\bf 170},
! 18--38. \\
!
! Hunke, E. C., and J. K. Dukowicz (2002).  The Elastic-Viscous-Plastic
! Sea Ice Dynamics Model in General Orthogonal Curvilinear Coordinates
! on a Sphere---Incorporation of Metric Terms. {\em Monthly Weather Review},
! {\bf 130}, 1848--1865.\\
!
! Hunke, E. C., and J. K. Dukowicz (2003).  The sea ice momentum
! equation in the free drift regime.  Los Alamos Tech. Rep. LA-UR-03-2219.\\
!
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! Summer 2003:  Vectorized by Clifford Chen (Fujitsu) and
!               William H. Lipscomb (LANL)
!
! !INTERFACE:
!
      module mod_evp
!
! !USES:
!
!
!EOP
!
      use mod_xc
      implicit none

      ! KAL - added here
      real, parameter :: 
     &   rhow    = 1026   ! Density of seawater
     &,  rhoa    = 1.3    ! Density of air

      ! KAL - time steps etc - added here - TODO: hardcoded now, modify
      ! logic a bit ...
      integer, parameter :: 
     &   ndyn_dt=1    ! dynamic time steps per thermo time step dyn_dt=dt/ndyn_dt


      integer ::
     &   kdyn         ! type of dynamics ( 1 = evp )
     &,  ndte         ! number of subcycles:  ndte=dyn_dt/dte

      ! KAL - time steps etc - added here - TODO: hardcoded now, modify
      ! logic a bit ...
      real ::
     &   dyn_dt       ! dynamic time step
     &,      dt       ! ice thermodynamic time step
     &,      dtei     !inverse subcycling time step (1/dte)

      logical ::
     &   evp_damping  ! if true, use evp damping procedure

      real, parameter ::
     &   dragw = 0.00536 * rhow 
                      ! drag coefficient for water on ice *rhow (kg/m^3)
     &,  draga = 0.00536 * 0.29 * rhoa
                      ! drag coefficient for water on ice *rhow (kg/m^3)

      real, parameter ::
     &   eyc = 0.36
                      ! coefficient for calculating the parameter E
     &,  cosw = 1.0   ! cos(ocean turning angle)  ! turning angle = 0
     &,  sinw = 0.0   ! sin(ocean turning angle)  ! turning angle = 0
     &,  cosa = 1.0   ! cos(atmos turning angle)  ! turning angle = 0
     &,  sina = 0.0   ! sin(atmos turning angle)  ! turning angle = 0

      real            ::  Pstar0
      real            ::  Cstar0

      real ::
     &   ecci         ! 1/e^2 
     &,  dte2T        ! dte/2T
     &,  denom        ! constants for stress equation
     &,  denom2       ! denom2 now only used in MIZ stress 

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &   shear        ! strain rate II component (1/s)
     &,  divu         ! strain rate I component, velocity divergence (1/s)
     &,  Delta        ! function of strain rates (1/s)

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &   waterx       ! for ocean stress calculation, x (m/s)
     &,  watery       ! for ocean stress calculation, y (m/s)
     &,  forcex       ! work array: combined atm stress and ocn tilt, x
     &,  forcey       ! work array: combined atm stress and ocn tilt, y
     &,  umassdtei    ! mass of U-cell/dte (kg/m^2 s)
     &,  fcor         ! Coriolis parameter (1/s)
     &,  fm           ! Coriolis param. * mass in U-cell (kg/s)
     &,  rcon_evp     ! for damping criterion (kg/s)
     &,  rcon_miz     ! for damping criterion in the MIZ (kg/s)
     &,  prss         ! pressure P (centered in T-cell) (kg/s)
     &,  prs_sig      ! replacement pressure, for stress calc
     &,  strintx      ! divergence of internal ice stress, x (N/m^2)
     &,  strinty      ! divergence of internal ice stress, y (N/m^2)
     &,  sig1         ! principal stress component (diagnostic)
     &,  sig2         ! principal stress component (diagnostic)
     &,  aiu          ! ice fraction on u-grid

      ! ice stress tensor in each corner of T cell (kg/s^2)
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &   stressp_1, stressp_2, stressp_3, stressp_4   ! sigma11+sigma22
     &,  stressm_1, stressm_2, stressm_3, stressm_4   ! sigma11-sigma22
     &,  stress12_1,stress12_2,stress12_3,stress12_4  ! sigma12

CKAL  ! -- For vectorized version. Not used yet
CKAL  !KAL integer (kind=int_kind), dimension (1:(ihi-ilo+1)*(jhi-jlo+1)) ::
CKAL  integer, dimension (1:(idm+2*nbdy)*(jdm+2*nbdy)) ::
CKAL &   indxti       ! compressed index in i-direction
CKAL &,  indxtj       ! compressed index in j-direction
CKAL &,  indxui       ! compressed index in i-direction
CKAL &,  indxuj       ! compressed index in j-direction

CKAL  ! -- For vectorized version. Not used yet
CKAL  integer ::
CKAL &   icellt       ! total count when icetmask is true
CKAL &,  icellu       ! total count when iceumask is true

      ! These are from the grid spec of cice
      logical, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &   tmask, ! land/boundary mask, thickness (T-cell)
     &   umask  ! land/boundary mask, velocity (U-cell)

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &   dxt        ! width of T-cell through the middle (m)
     &,  dyt        ! height of T-cell through the middle (m)
     &,  HTE        ! length of eastern edge of T-cell (m)
     &,  HTN        ! length of northern edge of T-cell (m)
     &,  dxhy       ! 0.5*(HTE - HTE)
     &,  dyhx       ! 0.5*(HTN - HTN)
     &,  dxt2       ! 0.5*dxt
     &,  dyt2       ! 0.5*dyt
     &,  dxt4       ! 0.25*dxt
     &,  dyt4       ! 0.25*dyt
     &,  dxu        ! width of U-cell through the middle (m)
     &,  dyu        ! height of U-cell through the middle (m)
     &,  ULAT       ! latitude of U-cell center
     &,  cyp        ! 1.5*HTE - 0.5*HTE
     &,  cxp        ! 1.5*HTN - 0.5*HTN
     &,  cym        ! 0.5*HTE - 1.5*HTE
     &,  cxm        ! 0.5*HTN - 1.5*HTN
     &,  tarea      ! area of T-cell (m^2)
     &,  uarea      ! area of U-cell (m^2)
     &,  tarear     ! 1/tarea
     &,  uarear     ! 1/uarea
     &,  tinyarea   ! puny*tarea





      logical, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &   icetmask     ! ice mask for t-cells
     &,  iceumask     ! ice mask for u-cells

      ! KAL From CSIM ice_state module
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &   uvel         ! x-component of velocity (m/s)
     &,  vvel         ! y-component of velocity (m/s)
     &,  uocn         ! x-component of velocity (m/s) - ocean
     &,  vocn         ! y-component of velocity (m/s) - ocean
     &,  uair         ! x-component of velocity (m/s) - air
     &,  vair         ! y-component of velocity (m/s) - air
     &,  tmass        ! Ice mass in t-cells
     &,  aice         ! Ice conc in t-cells
     &,  vice         ! Ice  volume in t-cells
     &,  vsno         ! Snow volume in t-cells
     &,  strength     ! ice strength (N/m)


      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &   strocnx      ! Ocean stress x - u-cell
     &,  strocny      ! Ocean stress y - u-cell
     &,  strocnxt     ! Ocean stress x - t-cell
     &,  strocnyt     ! Ocean stress y - t-cell
     &,  strairx      ! Ocean stress x - u-cell
     &,  strairy      ! Ocean stress y - u-cell
     &,  strtltx      ! Ocean tilt + stress  x - u-cell
     &,  strtlty      ! Ocean tilt + stress  y - u-cell
     &,  ss_tltx      ! Ocean tilt x - u-cell
     &,  ss_tlty      ! Ocean tilt y - u-cell
!    &,  strairxt     ! Ocean stress x - t-cell
!    &,  strairyt     ! Ocean stress y - t-cell
      integer, save :: imargin
      integer, save :: e_itst,e_jtst

      real, parameter ::
     &   gravit    = 9.80616  ! gravitational acceleration (m/s^2)
     &,  omega     = 7.292e-5 ! angular velocity of earth (rad/sec)

      real, parameter ::
     &  eps04  = 1.0e-4,
     &  eps11  = 1.0e-11,
     &  eps12  = 1.0e-12,
     &  eps13  = 1.0e-13,
     &  eps15  = 1.0e-15,
     &  eps16  = 1.0e-16,
     &  puny = eps13

      end module mod_evp

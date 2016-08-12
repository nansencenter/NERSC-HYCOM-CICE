   ! ======================================================================
   !  This module contains definitions of ice classes, grid variables (taken
   !  from HYCOM), the definition of types and the icestate variable containing
   !  all state variables in ICESTATE.
   !
   !  Also contains the definition of model parameters
   ! ======================================================================
   !Written by: David Salas y Melia  and  Ragnhild Hansen
   !KAL changed 18.01.2009 - Combined several modules

MODULE mod_icestate
   use mod_xc 
   implicit none

   ! Test parameters
   integer, parameter :: itst=13, jtst=30

   ! Albedo values
   REAL, save      :: albi_max = .80      ! Max albedo over dry ice []
   real, save      :: albi_mlt = .60      ! Albedo over melting ice []
   REAL, save      :: albsmin  = .75      ! Minimum albedo over snow []
   REAL, save      :: albsmax  = .85      ! Maximum albedo over snow []
   REAL, PARAMETER :: albw     = .065     ! Albedo over water []
   REAL, PARAMETER :: coalbw= 1.-albw     ! Water SW absorbtion

   !  Specific heats [J/(kg K)]
   REAL, PARAMETER :: cpice    = 2080.    ! ice 
   REAL, PARAMETER :: cpsw     = 3987.    ! water 
   REAL, PARAMETER :: cpw      = 4180.    ! fresh water 

   ! Phase transition heat parameters
   REAL, PARAMETER :: hofusn0  = 3.02E+08  ! Heat of fusion of ice [J/m^3]
   REAL, PARAMETER :: hofusni0= 1./hofusn0 ! Inverse of heat of fusion [m^3 / J]

   ! Ice and snow density [kg/m^3]
   real, parameter :: rhoice   = 910      ! Ice
   REAL, PARAMETER :: rhosnwmax= 330.     ! Maximum density for old snow Maykut (1971)
   !REAL, PARAMETER :: rhosnwmax_sum= 450.! Maximum density for old snow (summer) Maykut (1971)
   REAL, PARAMETER :: rhosnwmin= 150.     ! fresh snow 
   REAL, PARAMETER :: rhow     = 1000.    ! fresh water
   REAL, save      :: rhoref              ! Ocean reference specific density -- from OGCM
   REAL, save      :: thref               ! 1/rhoref [m^3/kg]

   ! Phase transition temperatures
   REAL, PARAMETER :: tice_m   = 273.05   ! Melting point of ice [K]
   REAL, PARAMETER :: t0deg    = 273.15   ! 0 C in [K]

   ! Decay coeffs for snow albedo
   REAL, PARAMETER :: taua     = .008     ! Linear decay coef. for snow albedo []
   REAL, PARAMETER :: tauf     = .24      ! Exponential decay coef. for snow albedo []

   ! Various constants
   real, parameter :: g=9.82              ! g ...
   REAL, PARAMETER :: rkice    = 2.04     ! Ice conductivity [W/(mK)]
   REAL, PARAMETER :: sice     = 6.       ! Salt concentrtation per mil for sea ice []


   ! "Tunable" parameters
   REAL, PARAMETER :: wnew     = 2.E-3    ! Minimum snowfall that allows snow 
                                          ! albedo to be refreshed to albsmax [m]
   real, save      :: snwlim  = .3        ! Limiting snow thickness
   real, save      :: qst_frac = 0.3      ! Max brine storage relative to total latent heat of ice
   real, save      :: fice_max = 0.999    ! Max allowed ice conc.
   real, parameter :: fice_min = 0.02     ! Min allowed ice conc.


   !Some useful parameters :
   real,    parameter :: pi2 = 6.283185307178 ! 2 times pi
   real,    parameter :: pi = 3.141592653589  !  pi
   REAL,    PARAMETER :: rmnthi   = 3.86E-07  ! Precipitations: 1m/month in m/s [m/s]
   real,    parameter :: radian=57.2957795


   REAL, PARAMETER :: aice     = 9.5      ! Vapor pressure first parameter over ice []
   REAL, PARAMETER :: bice     = 7.66     ! Vapor pressure second parameter over ice [K]
   REAL, PARAMETER :: awater   = 7.5      ! Vapor pressure first parameter over water []
   REAL, PARAMETER :: bwater   = 35.86    ! Vapor pressure second parameter over water [K]
   REAL, PARAMETER :: gasconst = .287E+03 ! Perfect gas constant [Pa.m^3/(K.kg)]

   REAL, PARAMETER :: emiss    = .97      ! Emissivity of water []
   REAL, PARAMETER :: stefanb  = 5.67E-08 ! Stefan-boltzman constant [W/(m^2 K^4)]
   REAL, PARAMETER :: cpair    = 1004.    ! Density of dry air
   REAL, PARAMETER :: hosubl   = 2.834E+06 ! Heat of sublimation [J/kg]

   REAL, PARAMETER :: hocond   = 2.5E+06   ! Heat of condensation [J/kg]
   REAL, PARAMETER :: hocondi = 1./hocond  ! Inverse of ice heat of condensation [kg/J]

   ! Array holding numerical values of lat. transf. coefficient
   real, save, dimension(16,29) :: clat


   ! Time specifics
   real, save :: dtt       ! Timestep in seconds for thermo routine

   ! Horizontal melt parameterization flags
   logical, save :: hmelt_hak,hmelt_solar


   ! Tolerance parameters moved from mod_icestate_thermo in here
   REAL, PARAMETER :: epsil0   = 1.E-6
   REAL, PARAMETER :: epsil1   = 1.E-10
   REAL, PARAMETER :: epsil2   = 1.E-20  


   ! ====================================================================
   ! ====================definition of ice classes ======================
   ! ====================================================================
   ! Vector thickl represents the boundaries of each ice class. For 
   ! example, ice class number i consists of ice whose thickness is 
   ! between thickl(i) and thickl(i+1). The thickest ice has boundaries
   ! thickl(nlay) and +infinity.
   ! Vector nlay defines vertical layers of each ice class

#if defined(ICESTATE_1CLASS)
   ! - 1 class
   integer, parameter              :: nthick=1   
   INTEGER, PARAMETER              :: hklim =1   ! Number of "thin" ice categories
   REAL, save, DIMENSION(nthick)   :: thickl    
   DATA thickl/ epsil1  /
   integer, parameter ::nlaymax = 1
   integer, parameter, dimension(nthick)   :: nlay=(/ 1  /)       

#elif defined (ICESTATE_2CLASS)
   ! - 2 classes
   integer, parameter              :: nthick=2
   INTEGER, PARAMETER              :: hklim =2   ! Number "thin" ice categories
   REAL, save, DIMENSION(nthick)   :: thickl     ! thickl describe the the lower vert.boundaries 
   DATA thickl/ epsil1 ,  1.0 /
   integer, parameter ::nlaymax = 3
   integer, parameter, dimension(nthick)   :: nlay=(/ 1 , 3  /)       

#elif defined (ICESTATE_5CLASS)
   ! - 5 classes
   integer, parameter              :: nthick = 5
   INTEGER, PARAMETER              :: hklim  = 2   ! Number of "thin" ice categories
   REAL, save, DIMENSION(nthick)   :: thickl      
   DATA thickl/ epsil1 , .5 , 1.0 , 2.0 , 5.0 /
   integer, parameter ::nlaymax = 5
   integer, parameter, dimension(nthick)   :: nlay=(/ 1 , 1 , 3 , 4 , 5 /)    
#if defined(SSNOWD)
    !Coefficient of variation (CV) for snow depth distribution.
    ! Each ice category has a given CV value based on observation
   real, parameter, dimension(nthick)     :: cvsnw=(/ 0.4 , 0.4 , 0.53 , 0.53 , 0.7 /) 
#endif
#else
#error - No CPP set for ICESTATE class selection
#endif

   !  - Type describing the state of the ice class
   type t_ice
      integer :: nlay                   ! Number of vtp layers
      real    :: qstore                 ! Energy stored in brine cells
      real    :: albs                   ! Albedo of snow
      real    :: fice                   ! Cell fraction with ice 
      real    :: hice                   ! Average thickness of ice 
      REAL    :: hsnw                   ! Average thickness of snow
      REAL    :: rhosnw                 ! Density of snow
      REAL    :: tsrf                   ! Temperature of surface 
      real, dimension(1:nlaymax) :: vtp ! Vert. temp profile of ice
#if defined(SSNOWD)
      real    :: hprcp                 ! Accumulated snow depth in the absence
                                       ! of any melt
      real    :: hmelt                 ! Accumulated snow melt depth
      real    :: cv                    ! Coefficient of variation (depends on the ice category)  
#endif     
#if defined(ICEAGE)
      real    :: age                     !Ice age (days)
#endif
   end type t_ice

   !  - Type describing the gridcell (ice and ml)
   type t_istate_cell
      type(t_ice),    dimension(nthick) :: ice
      real                              :: hml
      real                              :: sml
      real                              :: tml
   end type t_istate_cell

   ! The definition of "thin" and "thick" ice is given in module mod_icestate_init
   logical, save, dimension(nthick)                :: thin       ! True if ice is 'thin'

   !  - The array storing info on the ice:
   type(t_istate_cell), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), save :: icestate

   ! Latitude. (Used in solar flux and to deduce coriolis force). 
   REAL,    save, DIMENSION(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: rlatm      ! Latitude in radians

#if defined (ICEAGE)
    ! Average ice age for a given grid cell
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::               &
     &   ave_age
      real age_restart,age_max
#endif

#if defined (SSNOWD)
    ! Average snow cover fraction 
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::               &
     &   snow_cov_frac
#endif
end module mod_icestate

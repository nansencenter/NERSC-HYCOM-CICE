module m_icestate_exchange
!The theory here is presented in McPhee et al(1987) and Holland&Jenkins(1999).

implicit none
private
   real,parameter :: psi_n    = 0.052   ! An empirical constant
   real,parameter :: vonk     = 0.400   ! Von Karman constant
   real,parameter :: visc_mol = 1.95e-6 ! Molecular viscosity of water
   real,parameter :: minustar = 2.00e-3 ! corresponds to relative ice-ocean speed of 3 cm/s
   real,parameter :: prandtl  = 13.8    ! Prantdl Number (visc / thermal diff)
   real,parameter :: etastar  = 1.      ! Stability parameter(McPhee et al. 1987), set to one here.

   public :: heatexch_coeff
contains



   real function heatexch_coeff(ustar,fcor) 
      implicit none
      real,intent(in) :: ustar,fcor
      real :: den

      ! Exchange coeff From Holland and Jenkins(1999). 
      ! NB: coriolis sign unimportant
      den = ( log(max(ustar,minustar)*psi_n*etastar**2)   -   &
              log(abs(fcor)*sublayer_thickness(ustar)      )   +   &
              vonk/(2*psi_n*etastar)   -    1. ) /vonk
      !print *,'Denominator= ',den,ustar,fcor
      den = den + 12.5*prandtl**(.67) - 6
      den = max(den,.1)
      heatexch_coeff = max(ustar,minustar)/den
   end function heatexch_coeff



   real function sublayer_thickness(ustar)
      implicit none
      real,intent(in) :: ustar

      ! Viscous sublayer thickness from Holland and Jenkins(1999)
      sublayer_thickness = 5.*visc_mol/max(minustar,ustar)
   end function sublayer_thickness


end module  m_icestate_exchange








      















   






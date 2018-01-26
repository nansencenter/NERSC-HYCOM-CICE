! Copyright 2017 Helmholtz-Zentrum Geesthacht
! Licensed under the Apache License, Version 2.0
!  (http://www.apache.org/licenses/LICENSE-2.0)
! Author(s): Richard Hofmeister, Ute Daewel

#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_nersc_ecosmo --- ECOSMO biogeochemical model developed by NERSC
!
! !INTERFACE:
   module ecosmo_bg
!
! !DESCRIPTION:
!
! The ECOSMO model is based on Daewel & Schrum (JMS,2013)
!
! !USES:
   use fabm_types
   use fabm_expressions

   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_nersc_ecosmo_bg
!
! !PRIVATE DATA MEMBERS:
   real(rk), parameter :: secs_pr_day = 86400.0_rk
   real(rk), parameter :: sedy0 = 86400.0_rk
   real(rk), parameter :: mmolm3_in_mll = 44.6608009_rk
   real(rk)            :: redf(20)=0.0_rk
   real(rk)            :: PrmBioC(45)=0.0_rk,BioC(45)=0.0_rk
   real(rk)            :: PrmGI(2,6)=0.0_rk,GI(2,6)=0.0_rk
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model) :: type_nersc_ecosmo_bg
!     Variable identifiers
      type (type_state_variable_id)         :: id_no3, id_nh4, id_pho, id_sil
      type (type_state_variable_id)         :: id_det, id_dom, id_oxy
      type (type_state_variable_id)         :: id_bg
      type (type_dependency_id)             :: id_temp, id_salt, id_par
      type (type_dependency_id)             :: id_parmean
      type (type_horizontal_dependency_id)  :: id_sfpar, id_meansfpar
      type (type_diagnostic_variable_id)    :: id_primprod
      type (type_diagnostic_variable_id)    :: id_parmean_diag

!     Model parameters
      real(rk) :: BioC(45)
      real(rk) :: nfixation_minimum_daily_par
      real(rk) :: bg_growth_minimum_daily_rad
      real(rk) :: frr
      contains

!     Model procedures
      procedure :: initialize
      procedure :: do

   end type type_nersc_ecosmo_bg
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the ECOSMO model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the ecosmo namelist is read and the variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_nersc_ecosmo_bg),intent(inout),target  :: self
   integer,                intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister
!
! !LOCAL VARIABLES:
   real(rk) :: no3_init, nh4_init,pho_init,oxy_init,&
               bg_init, det_init, dom_init
   real(rk) :: frr
   real(rk) :: nfixation_minimum_daily_par=40.0
   real(rk) :: bg_growth_minimum_daily_rad=120.0

   character(len=attribute_length) :: no3_name
! add for ph
   integer  :: i

   namelist /nersc_ecosmo_bg/  PrmBioC, &
                          bg_init,no3_name, &
                          nfixation_minimum_daily_par, bg_growth_minimum_daily_rad

!EOP
!-----------------------------------------------------------------------
!BOC

!include 'ECOSMparamNSBS.f'

frr = 0.4_rk

   ! Read the namelist
   if (configunit>0) read(configunit,nml=nersc_ecosmo_bg,err=99,end=100)

   self%frr = frr
   self%nfixation_minimum_daily_par = nfixation_minimum_daily_par
   self%bg_growth_minimum_daily_rad = bg_growth_minimum_daily_rad

   ! set Redfield ratios:
   redf(1) = 6.625_rk      !C_N
   redf(2) = 106.0_rk      !C_P
   redf(3) = 6.625_rk      !C_SiO
   redf(4) = 16.0_rk       !N_P
   redf(5) = 1.0_rk        !N_SiO
   redf(6) = 12.01_rk      !C_Cmg
   redf(7) = 44.6608009_rk !O2mm_ml
   redf(8) = 14.007_rk     !N_Nmg
   redf(9) = 30.97_rk      !P_Pmg
   redf(10) = 28.09_rk     !Si_Simg
   do i=1,10
     redf(i+10) = 1._rk/redf(i)
   end do

   !  change units 1/day to 1/sec and mmolN,P,Si to mmolC
           BioC(1) =   PrmBioC(1) /sedy0                    !  1/day
           BioC(2) =   PrmBioC(2) /sedy0                    !  1/day
           BioC(3) =   PrmBioC(3)                           !  m**2/W
           BioC(4) =   PrmBioC(4)                           !  1/m
           BioC(5) =   PrmBioC(5)/(REDF(1)*REDF(6)) !  m**2/mmolN
           BioC(6) =   PrmBioC(6)*REDF(1)*REDF(6)   !  mmolN/m**3
           BioC(7) =   PrmBioC(7)*REDF(1)*REDF(6)   !  mmolN/m**3
           BioC(8) =   PrmBioC(8)/(REDF(1)*REDF(6)) !  m**3/mmolN
           BioC(9) =   PrmBioC(9) /sedy0                    !  1/day
           BioC(10)=   PrmBioC(10)/sedy0                    !  1/day
           BioC(11)=   PrmBioC(11)/sedy0                    !  1/day
           BioC(12)=   PrmBioC(12)/sedy0                    !  1/day
           BioC(13)=   PrmBioC(13)/sedy0                    !  1/day
           BioC(14)=   PrmBioC(14)*REDF(1)*REDF(6)   !  mmolN/m**3
           BioC(15)=   PrmBioC(15)/sedy0                    !  1/day
           BioC(16)=   PrmBioC(16)/sedy0                    !  1/day
           BioC(17)=   PrmBioC(17)/sedy0                    !  1/day
           BioC(18)=   PrmBioC(18)/sedy0                    !  1/day
           BioC(19)=   PrmBioC(19)                          !  1
           BioC(20)=   PrmBioC(20)                          !  1
           BioC(21)=   PrmBioC(21)                          !  1

           BioC(22)=   PrmBioC(22)/sedy0                    !   1/day
           BioC(23)=   PrmBioC(23)/sedy0                    !   m/day
           BioC(24)=   PrmBioC(24)/sedy0                    !   m/day
! set as mmolN/m**3 =0.25  !!!! need a correction
           BioC(25)=   PrmBioC(25)*REDF(2)*REDF(6)   ! mmolP/m**3
           BioC(26)=   PrmBioC(26)*REDF(3)*REDF(6)   ! mmolSi/m**3
           BioC(27)=   PrmBioC(27)/sedy0                    !   1/day

           BioC(28)=   PrmBioC(28)/sedy0                    !  1/day

           BioC(29)=   PrmBioC(29)                          ! 1/degC

           BioC(30)=   PrmBioC(30)                         ! degC

           BioC(31)=   PrmBioC(31)/sedy0                   ! 1/d

           BioC(32)=   PrmBioC(32)/sedy0                   ! 1/d

           BioC(33)=   PrmBioC(33)/sedy0                   ! m/d

           BioC(34)=   PrmBioC(34)                         !  N/m**2

           BioC(35)=   PrmBioC(35)/sedy0                   ! 1/day

           BioC(36)=   PrmBioC(36)/sedy0                   ! m/day

           BioC(37)=   PrmBioC(37)/sedy0                   ! 1/day

           BioC(38)=   PrmBioC(38)/sedy0                   ! 1/day

           BioC(39)=   PrmBioC(39)                         ! 1/degC

           BioC(40)=   PrmBioC(40)

           BioC(41)=   PrmBioC(41)

           BioC(42)=   PrmBioC(42)/sedy0

           BioC(43)=   PrmBioC(43)/sedy0                    !   m/day
           BioC(44)=   PrmBioC(44)/sedy0                    !   m/day
           BioC(45)=   PrmBioC(45)/sedy0                    !   m/day

   ! Store parameter values in our own derived type
   do i=1,size(BioC)
     self%BioC(i) = BioC(i)
   end do
   !  growth fractions

   ! Register state variables
   call self%register_state_dependency(self%id_no3,'no3','mgC/m3','nitrate')
       call self%request_coupling(self%id_no3,no3_name)
   call self%register_state_dependency(self%id_nh4,'nh4','mgC/m3','ammonium')     

   call self%register_state_dependency(self%id_pho,'pho','mgC/m3','phosphate')     

   call self%register_state_variable(self%id_bg,'bg','mgC/m3','cyanobacteria',     &
                                     initial_value=bg_init*redf(1)*redf(6), &
                                     vertical_movement=-BioC(44), &
                                     minimum=1.0e-07_rk,specific_light_extinction=self%BioC(5))
                                     ! CAGLAR: Jorn (FABM developer) assured me that if we do this way
                                     ! BG light extinction will be taken into account

   call self%register_state_dependency(self%id_det,'det','mgC/m3','detritus')     

   call self%register_state_dependency(self%id_dom,'dom','mgC/m3','labile dissolved om')
   
   call self%register_state_dependency(self%id_oxy,'oxy','mmolO2/m3','oxygen')

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_primprod,'primprod','mgC/m**3/s', &
         'primary production rate')
   call self%register_diagnostic_variable(self%id_parmean_diag,'parmean','W/m**2', &
         'daily-mean photosynthetically active radiation')

   ! Register dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_salt,standard_variables%practical_salinity)
   call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_sfpar,standard_variables%surface_downwelling_photosynthetic_radiative_flux)
   ! use temporal mean of light for the last 24 hours
   call self%register_dependency(self%id_parmean,temporal_mean(self%id_par,period=86400._rk,resolution=3600._rk))
   call self%register_dependency(self%id_meansfpar,temporal_mean(self%id_sfpar,period=86400._rk,resolution=3600._rk))

   return

99 call self%fatal_error('nersc_ecosmo_initialize','Error reading namelist nersc_ecosmo.')

100 call self%fatal_error('nersc_ecosmo_initialize','Namelist nersc_ecosmo was not found.')

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of ECOSMO model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_nersc_ecosmo_bg),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister
!
! !LOCAL VARIABLES:
   real(rk) :: no3,nh4,pho,oxy
   real(rk) :: bg
   real(rk) :: det,dom
   real(rk) :: temp,salt,par
   real(rk) :: blight
   real(rk) :: Tbg
   real(rk) :: Prod,Bg_prod
   real(rk) :: Fs,Fl,ZlonBg,ZsonBg
   real(rk) :: up_no3,up_nh4,up_n,up_pho
   real(rk) :: rhs,dxxdet ! what is rhs?
   real(rk) :: mean_par, mean_surface_par, Bg_fix
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_temp,temp)
   _GET_(self%id_salt,salt)
   _GET_(self%id_par,par)
   _GET_(self%id_no3,no3)
   _GET_(self%id_nh4,nh4)
   _GET_(self%id_pho,pho)
   _GET_(self%id_bg,bg)
!   _GET_(self%id_det,det)
!   _GET_(self%id_dom,dom)
!   _GET_(self%id_oxy,oxy)
   _GET_(self%id_parmean,mean_par)
   _GET_HORIZONTAL_(self%id_meansfpar,mean_surface_par)

   ! nutrient limitation factors
   ! k denotes half-saturation values
   up_nh4 = nh4/(BioC(6)+nh4)
   up_no3 = no3/(BioC(7)+no3)*exp(-BioC(8)*nh4)
   up_n = up_nh4+up_no3
   up_pho = pho/(BioC(25)+pho)

   ! light limitation
   blight = max(tanh(BioC(3)*par),0.0_rk)

   ! temperature dependence
   if ((salt<=10.0) .and. (mean_surface_par > self%bg_growth_minimum_daily_rad)) then
     Tbg = 1.0_rk/(1.0_rk + exp(BioC(29)*(BioC(30)-temp)))
   else
     Tbg = 0.0_rk
   end if

   ! production and nutrient uptake
   !Bg_prod = 0.0_rk ! the light criterium restricts growth to surface waters
   Bg_prod = Tbg * min(blight, up_n, up_pho)
   Bg_fix=0.0_rk
   if (mean_par > self%nfixation_minimum_daily_par) then
     Bg_fix = Tbg * min(blight, up_pho) - Bg_prod
   end if
   Prod = BioC(28)*Bg_prod*bg ! cyanobacteria production
        ! production due to flagellates and diatoms are
        ! calculated in the main main ECOSMO code



! reaction rates
   _SET_ODE_( self%id_bg,  (BioC(28)*(Bg_prod + Bg_fix) - BioC(32))*bg )

   ! detritus
   dxxdet = ( BioC(32) * bg )

   rhs = (1.0_rk-self%frr) * dxxdet
   _SET_ODE_(self%id_det, rhs)

   ! labile dissolved organic matter
   _SET_ODE_(self%id_dom, self%frr*dxxdet)

   ! nitrate
   rhs = -(up_no3+0.5d-10)/(up_n+1.0d-10)*Prod
   _SET_ODE_(self%id_no3, rhs)

   ! ammonium
   rhs = -(up_nh4+0.5d-10)/(up_n+1.0d-10)*Prod
   _SET_ODE_(self%id_nh4, rhs)

   ! phosphate
   _SET_ODE_(self%id_pho, -Prod -BioC(28)*bg*Bg_fix )


! CAGLAR - FIX THIS
   ! oxygen
   rhs = ((6.625*up_nh4 + 8.125*up_no3+1.d-10))/(up_n+1.d-10)*Prod
   _SET_ODE_(self%id_oxy, rhs)

   ! Export diagnostic variables

   _SET_DIAGNOSTIC_(self%id_primprod, Prod + BioC(28)*bg*Bg_fix )

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC



 end module ecosmo_bg

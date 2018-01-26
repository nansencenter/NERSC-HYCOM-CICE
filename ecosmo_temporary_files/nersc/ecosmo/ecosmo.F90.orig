! Copyright 2017 Helmholtz-Zentrum Geesthacht
! Licensed under the Apache License, Version 2.0
!  (http://www.apache.org/licenses/LICENSE-2.0)
! Author(s): Richard Hofmeister, Ute Daewel

#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_hzg_ecosmo --- ECOSMO biogeochemical model
!
! !INTERFACE:
   module fabm_hzg_ecosmo
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
   public type_hzg_ecosmo
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
   type,extends(type_base_model) :: type_hzg_ecosmo
!     Variable identifiers
      type (type_state_variable_id)         :: id_no3, id_nh4, id_pho, id_sil
      type (type_state_variable_id)         :: id_opa, id_det, id_dia, id_fla
      type (type_state_variable_id)         :: id_mesozoo, id_microzoo, id_bg, id_dom, id_oxy
      type (type_bottom_state_variable_id)  :: id_sed1, id_sed2, id_sed3
      type (type_dependency_id)             :: id_temp, id_salt, id_par
      type (type_dependency_id)             :: id_parmean
      type (type_horizontal_dependency_id)  :: id_tbs
      type (type_horizontal_dependency_id)  :: id_sfpar, id_meansfpar
      type (type_diagnostic_variable_id)    :: id_denit, id_primprod, id_secprod
      type (type_diagnostic_variable_id)    :: id_parmean_diag

!     Model parameters
      real(rk) :: BioC(45)
      real(rk) :: zpr, frr
      real(rk) :: gf_zs_ps
      real(rk) :: gf_zs_pl
      real(rk) :: gf_zs_bg
      real(rk) :: gf_zs_det
      real(rk) :: gf_zl_ps
      real(rk) :: gf_zl_pl
      real(rk) :: gf_zl_bg
      real(rk) :: gf_zl_det
      real(rk) :: gf_zl_zs
      real(rk) :: surface_deposition_no3
      real(rk) :: surface_deposition_nh4
      real(rk) :: surface_deposition_pho
      real(rk) :: surface_deposition_sil
      real(rk) :: nfixation_minimum_daily_par
      real(rk) :: bg_growth_minimum_daily_rad

      contains

!     Model procedures
      procedure :: initialize
      procedure :: do
      procedure :: do_surface
      procedure :: do_bottom
      procedure :: get_light_extinction

   end type type_hzg_ecosmo
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
   class (type_hzg_ecosmo),intent(inout),target  :: self
   integer,                intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister
!
! !LOCAL VARIABLES:
   real(rk) :: no3_init, nh4_init,pho_init,sil_init,oxy_init, &
               dia_init,fla_init,bg_init, &
               mesozoo_init,microzoo_init, &
               det_init, dom_init, opa_init, &
               sed1_init, sed2_init, sed3_init
   real(rk) :: zpr, frr
   real(rk) :: surface_deposition_no3=0.0
   real(rk) :: surface_deposition_nh4=0.0
   real(rk) :: surface_deposition_pho=0.0
   real(rk) :: surface_deposition_sil=0.0
   real(rk) :: nfixation_minimum_daily_par=40.0
   real(rk) :: bg_growth_minimum_daily_rad=120.0


   integer  :: i

   namelist /hzg_ecosmo/  PrmBioC,GI,zpr,frr, &
                          no3_init,nh4_init,pho_init, &
                          sil_init,oxy_init, &
                          dia_init,fla_init,bg_init, &
                          mesozoo_init,microzoo_init, &
                          det_init, dom_init, opa_init, &
                          sed1_init, sed2_init, sed3_init, &
                          surface_deposition_no3, surface_deposition_nh4, &
                          surface_deposition_pho, surface_deposition_sil, &
                          nfixation_minimum_daily_par, bg_growth_minimum_daily_rad

!EOP
!-----------------------------------------------------------------------
!BOC

   !include 'ECOSMparamNSBS.f'

   zpr = 0.001_rk
   frr = 0.4_rk

   ! Read the namelist
   if (configunit>0) read(configunit,nml=hzg_ecosmo,err=99,end=100)

   self%zpr = zpr/secs_pr_day
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

   ! set surface fluxes in [mgC/m2/s]
   self%surface_deposition_no3 = surface_deposition_no3*redf(1)*redf(6)/secs_pr_day
   self%surface_deposition_nh4 = surface_deposition_nh4*redf(1)*redf(6)/secs_pr_day
   self%surface_deposition_pho = surface_deposition_pho*redf(2)*redf(6)/secs_pr_day
   self%surface_deposition_sil = surface_deposition_sil*redf(3)*redf(6)/secs_pr_day

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
   self%gf_zs_ps  = GI(1,1)
   self%gf_zs_pl  = GI(1,2)
   self%gf_zs_bg  = GI(1,6)
   self%gf_zs_det = GI(1,5)
   self%gf_zl_ps  = GI(2,1)
   self%gf_zl_pl  = GI(2,2)
   self%gf_zl_bg  = GI(2,6)
   self%gf_zl_det = GI(2,5)
   self%gf_zl_zs  = GI(2,3)

   ! Register state variables
   call self%register_state_variable(self%id_no3,'no3','mgC/m3','nitrate',     &
                                     initial_value=no3_init*redf(1)*redf(6), &
                                     vertical_movement=0.0_rk, &
                                     minimum=0.0_rk)

   call self%register_state_variable(self%id_nh4,'nh4','mgC/m3','ammonium',     &
                                     initial_value=nh4_init*redf(1)*redf(6), &
                                     vertical_movement=0.0_rk, &
                                     minimum=0.0_rk)

   call self%register_state_variable(self%id_pho,'pho','mgC/m3','phosphate',     &
                                     initial_value=pho_init*redf(2)*redf(6), &
                                     vertical_movement=0.0_rk, &
                                     minimum=0.0_rk)

   call self%register_state_variable(self%id_sil,'sil','mgC/m3','silicate',     &
                                     initial_value=sil_init*redf(3)*redf(6), &
                                     vertical_movement=0.0_rk, &
                                     minimum=0.0_rk)

   call self%register_state_variable(self%id_oxy,'oxy','mmolO2/m3','oxygen',     &
                                     initial_value=oxy_init, &
                                     vertical_movement=0.0_rk, &
                                     minimum=-100000._rk)

   call self%register_state_variable(self%id_dia,'dia','mgC/m3','large phytoplankton',     &
                                     initial_value=dia_init*redf(1)*redf(6), &
                                     vertical_movement=-BioC(45), &
                                     minimum=1.0e-07_rk)

   call self%register_state_variable(self%id_fla,'fla','mgC/m3','small phytoplankton',     &
                                     initial_value=fla_init*redf(1)*redf(6), &
                                     vertical_movement=0.0_rk, &
                                     minimum=1.0e-07_rk)

   call self%register_state_variable(self%id_bg,'bg','mgC/m3','cyanobacteria',     &
                                     initial_value=bg_init*redf(1)*redf(6), &
                                     vertical_movement=-BioC(44), &
                                     minimum=1.0e-07_rk)

   call self%register_state_variable(self%id_microzoo,'microzoo','mgC/m3','microzooplankton',     &
                                     initial_value=microzoo_init*redf(1)*redf(6), &
                                     vertical_movement=0.0_rk, &
                                     minimum=1.0e-7_rk)

   call self%register_state_variable(self%id_mesozoo,'mesozoo','mgC/m3','mesozooplankton',     &
                                     initial_value=mesozoo_init*redf(1)*redf(6), &
                                     vertical_movement=0.0_rk, &
                                     minimum=1.0e-7_rk)

   call self%register_state_variable(self%id_det,'det','mgC/m3','detritus',     &
                                     initial_value=det_init*redf(1)*redf(6), &
                                     vertical_movement=-BioC(23), &
                                     minimum=0.0_rk)

   call self%register_state_variable(self%id_opa,'opa','mgC/m3','opal',     &
                                     initial_value=opa_init*redf(3)*redf(6), &
                                     vertical_movement=-BioC(43), &
                                     minimum=0.0_rk)

   call self%register_state_variable(self%id_dom,'dom','mgC/m3','labile dissolved om',     &
                                     initial_value=dom_init*redf(1)*redf(6), &
                                     vertical_movement=0.0_rk, &
                                     minimum=0.0_rk)

   call self%register_state_variable(self%id_sed1,'sed1','mgC/m2','sediment detritus',     &
                                     initial_value=sed1_init*redf(1)*redf(6), &
                                     minimum=0.0_rk)

   call self%register_state_variable(self%id_sed2,'sed2','mgC/m2','sediment opal',     &
                                     initial_value=sed2_init*redf(3)*redf(6), &
                                     minimum=0.0_rk)

   call self%register_state_variable(self%id_sed3,'sed3','mgC/m2','sediment adsorbed phosporus',     &
                                     initial_value=sed3_init*redf(2)*redf(6), &
                                     minimum=0.0_rk)


   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_denit,'denit','mmolN/m**3/s', &
         'denitrification rate', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_primprod,'primprod','mgC/m**3/s', &
         'primary production rate', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_secprod,'secprod','mgC/m**3/s', &
         'secondary production rate', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_parmean_diag,'parmean','W/m**2', &
         'daily-mean photosynthetically active radiation', output=output_time_step_averaged)

   ! Register dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_salt,standard_variables%practical_salinity)
   call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_tbs,standard_variables%bottom_stress)
   call self%register_dependency(self%id_sfpar,standard_variables%surface_downwelling_photosynthetic_radiative_flux)
   ! use temporal mean of light for the last 24 hours
   call self%register_dependency(self%id_parmean,temporal_mean(self%id_par,period=86400._rk,resolution=3600._rk))
   call self%register_dependency(self%id_meansfpar,temporal_mean(self%id_sfpar,period=86400._rk,resolution=3600._rk))

   return

99 call self%fatal_error('hzg_ecosmo_initialize','Error reading namelist hzg_ecosmo.')

100 call self%fatal_error('hzg_ecosmo_initialize','Namelist hzg_ecosmo was not found.')

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
   class (type_hzg_ecosmo),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister
!
! !LOCAL VARIABLES:
   real(rk) :: no3,nh4,pho,sil,t_sil,oxy,fla,bg,dia
   real(rk) :: microzoo,mesozoo,opa,det,dom
   real(rk) :: temp,salt,par
   real(rk) :: frem, fremDOM, blight
   real(rk) :: Ts,Tl,Tbg
   real(rk) :: Prod,Ps_prod,Pl_prod,Bg_prod
   real(rk) :: Fs,Fl,ZlonPs,ZlonPl,ZsonD,ZlonD,ZlonBg,ZsonBg,ZsonPs,ZsonPl,ZlonZs
   real(rk) :: up_no3,up_nh4,up_n,up_pho,up_sil
   real(rk) :: bioom1,bioom2,bioom3,bioom4,bioom5,bioom6,bioom7,bioom8,Onitr
   real(rk) :: rhs,dxxdet
   real(rk) :: Zl_prod, Zs_prod
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
   _GET_(self%id_sil,sil)
   _GET_(self%id_dia,dia)
   _GET_(self%id_fla,fla)
   _GET_(self%id_bg,bg)
   _GET_(self%id_microzoo,microzoo)
   _GET_(self%id_mesozoo,mesozoo)
   _GET_(self%id_det,det)
   _GET_(self%id_dom,dom)
   _GET_(self%id_opa,opa)
   _GET_(self%id_oxy,oxy)
   _GET_(self%id_parmean,mean_par)
   _GET_HORIZONTAL_(self%id_meansfpar,mean_surface_par)

   ! remineralisation rate
   frem = BioC(22) * (1._rk+20._rk*(temp**2/(13._rk**2+temp**2)))
   fremDOM = 10._rk * frem

   ! nutrient limitation factors
   ! k denotes half-saturation values
   up_nh4 = nh4/(BioC(6)+nh4)
   up_no3 = no3/(BioC(7)+no3)*exp(-BioC(8)*nh4)
   up_n = up_nh4+up_no3
   up_pho = pho/(BioC(25)+pho)
   t_sil = max(sil-80._rk,0.0_rk)
   up_sil = t_sil/(BioC(26)+t_sil)

   ! light limitation
   blight = max(tanh(BioC(3)*par),0.0_rk)

   ! temperature dependence
   Ts = 1.0_rk
   Tl = 1.0_rk
   if ((salt<=10.0) .and. (mean_surface_par > self%bg_growth_minimum_daily_rad)) then
     Tbg = 1.0_rk/(1.0_rk + exp(BioC(29)*(BioC(30)-temp)))
   else
     Tbg = 0.0_rk
   end if

   ! production and nutrient uptake
   Ps_prod = Ts * min(blight, up_n, up_pho)
   Pl_prod = Tl * min(blight, up_n, up_pho, up_sil)
   !Bg_prod = 0.0_rk ! the light criterium restricts growth to surface waters
   Bg_prod = Tbg * min(blight, up_n, up_pho)
   Bg_fix=0.0_rk
   if (mean_par > self%nfixation_minimum_daily_par) then
     Bg_fix = Tbg * min(blight, up_pho) - Bg_prod
   end if
   Prod = BioC(1)*Pl_prod*dia + & ! diatoms production
          BioC(2)*Ps_prod*fla + & ! flagellates production
          BioC(28)*Bg_prod*bg ! cyanobacteria production

   ! grazing
   ! gf denotes grazing fraction

   Fs = self%gf_zs_ps*fla + self%gf_zs_pl*dia + self%gf_zs_det*det + self%gf_zs_bg*bg
   Fl = self%gf_zl_ps*fla + self%gf_zl_pl*dia + self%gf_zl_zs*microzoo + &
          self%gf_zl_det*det + self%gf_zl_bg*bg

   ZsonPs = BioC(12) * self%gf_zs_ps * fla/(BioC(14) + Fs)
   ZsonPl = BioC(12) * self%gf_zs_pl * dia/(BioC(14) + Fs)
   ZsonD = BioC(12) * self%gf_zs_det * det/(BioC(14) + Fs)
   ZsonBg = BioC(31) * self%gf_zs_bg * bg/(BioC(14) + Fs)

   ZlonPs = BioC(11) * self%gf_zl_ps * fla/(BioC(14) + Fl)
   ZlonPl = BioC(11) * self%gf_zl_pl * dia/(BioC(14) + Fl)
   ZlonD = BioC(11) * self%gf_zl_det * det/(BioC(14) + Fl)
   ZlonZs = BioC(13) * self%gf_zl_zs * microzoo/(BioC(14) + Fl)
   ZlonBg = BioC(31) * self%gf_zl_bg * bg/(BioC(14) + Fl)

   ! nitrification
   Onitr = 0.01_rk * redf(7) !according to Neumann  (Onitr in mlO2/l see also Stigebrand and Wulff)
   bioom1 = 0.0_rk
   bioom2 = 0.0_rk
   bioom3 = 0.0_rk
   bioom4 = 0.0_rk
   bioom5 = 0.0_rk
   bioom6 = 0.0_rk
   bioom7 = 0.0_rk
   bioom8 = 0.0_rk
   if (oxy > 0) then
     bioom1 = 0.1_rk/secs_pr_day * exp(temp*0.11_rk) * oxy/(Onitr+oxy)
     bioom2 = bioom1
     bioom6 = 1.0_rk
   else
     if (no3>0) then
       bioom5 = 5.0_rk
       bioom8 = 1.0_rk
     else
       bioom7 = 1.0_rk
     end if
   end if

! reaction rates
   _SET_ODE_(self%id_fla, (BioC(2)*Ps_prod - BioC(10))*fla - ZsonPs*microzoo - ZlonPs*mesozoo)
   _SET_ODE_(self%id_dia, (BioC(1)*Pl_prod - BioC(9))*dia - ZsonPl*microzoo - ZlonPl*mesozoo)
   _SET_ODE_(self%id_bg,  (BioC(28)*(Bg_prod + Bg_fix) - BioC(32))*bg - ZsonBg*microzoo - ZlonBg*mesozoo)

   ! microzooplankton
   Zs_prod = BioC(20)*(ZsonPs + ZsonPl + ZsonBg) + BioC(21)*ZsonD
   rhs = (Zs_prod - BioC(16) - BioC(18) - self%zpr) * microzoo &
         - ZlonZs * mesozoo
   _SET_ODE_(self%id_microzoo, rhs)

   ! mesozooplankton
   Zl_prod = BioC(19)*(ZlonPs + ZlonPl + ZlonBg + ZlonZs) + BioC(21)*ZlonD 
   rhs = (Zl_prod - BioC(15) - BioC(17) - self%zpr) * mesozoo
   _SET_ODE_(self%id_mesozoo, rhs)

   ! detritus
   dxxdet = (  ((1.0_rk-BioC(20))*(ZsonPs + ZsonPl + ZsonBg) & 
              + (1.0_rk-BioC(21))*ZsonD) * microzoo &
              + ((1.0_rk-BioC(19))*(ZlonPs + ZlonPl + ZlonBg + ZlonZs) &
              + (1.0_rk-BioC(21))*ZlonD) * mesozoo &
              + BioC(16) * microzoo &
              + BioC(15) * mesozoo &
              + BioC(10) * fla &
              + BioC(9)  * dia &
              + BioC(32) * bg )
   
   rhs = (1.0_rk-self%frr) * dxxdet &
         - ZsonD * microzoo &
         - ZlonD * mesozoo &
         - frem * det
   _SET_ODE_(self%id_det, rhs)

   ! labile dissolved organic matter
   _SET_ODE_(self%id_dom, self%frr*dxxdet - fremdom * dom)

   ! nitrate
   rhs = -(up_no3+0.5d-10)/(up_n+1.0d-10)*Prod &
         + bioom1 * nh4 &
         - bioom3 * no3 &
         - frem * det * bioom5 &
         - fremDOM * dom * bioOM5
   _SET_ODE_(self%id_no3, rhs)

   ! ammonium
   rhs = -(up_nh4+0.5d-10)/(up_n+1.0d-10)*Prod &
         + BioC(18) * microzoo &
         + BioC(17) * mesozoo &
         + frem * det &
         + fremDOM * dom - bioom1 * nh4
   _SET_ODE_(self%id_nh4, rhs)

   ! phosphate
   _SET_ODE_(self%id_pho, -Prod -BioC(28)*bg*Bg_fix + BioC(18) * microzoo + BioC(17) * mesozoo + frem*det + fremDOM*dom)

   ! silicate
   _SET_ODE_(self%id_sil, -BioC(1)*Pl_prod*dia + BioC(27)*opa)

   ! opal
   _SET_ODE_(self%id_opa, BioC(9)*dia + ZsonPl*microzoo + ZlonPl*mesozoo - BioC(27)*opa)

   ! oxygen
   rhs = ((6.625*up_nh4 + 8.125*up_no3+1.d-10)/(up_n+1.d-10)*Prod &
         -bioom6*6.625*(BioC(18)*microzoo + BioC(17)*mesozoo) &
         -frem*det*(bioom6+bioom7)*6.625 &
         -(bioom6+bioom7)*6.625*fremDOM*dom &
         -2.0_rk*bioom1*nh4)*redf(11)*redf(16)
   _SET_ODE_(self%id_oxy, rhs)

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_denit,(frem*det*bioom5+fremDOM*dom*bioom5)*redf(11)*redf(16))
   _SET_DIAGNOSTIC_(self%id_primprod, Prod + BioC(28)*bg*Bg_fix )
   _SET_DIAGNOSTIC_(self%id_secprod, Zl_prod*mesozoo + Zs_prod*microzoo)
   _SET_DIAGNOSTIC_(self%id_parmean_diag, mean_par)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC






!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Surface fluxes for the ecosmo model
!
! !INTERFACE:

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
   class (type_hzg_ecosmo),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_
!
! !LOCAL VARIABLES:
   real(rk) :: o2flux, T, tr, S, o2sat, oxy
   real(rk) :: no3flux, phoflux
   real(rk) :: pho,par,bg,blight,tbg,up_pho,prod
!EOP
!-----------------------------------------------------------------------
!BOC
   _HORIZONTAL_LOOP_BEGIN_

   _GET_(self%id_temp,T)
   _GET_(self%id_salt,S)
   _GET_(self%id_oxy,oxy)
   _GET_(self%id_par,par)
   _GET_(self%id_pho,pho)
   _GET_(self%id_bg,bg)

   ! Oxygen saturation micromol/liter__(Benson and Krause, 1984)
   tr = 1.0_rk/(T + 273.15_rk)
   o2sat= exp(- 135.90205_rk              &
       + (1.575701d05 ) * tr               &
       - (6.642308d07 ) * tr**2            &
       + (1.243800d10) * tr**3            &
       - (8.621949d11) * tr**4            &
       - S*(0.017674_rk-10.754_rk*tr+2140.7_rk*tr**2)  )

   o2flux = 5._rk/secs_pr_day * (o2sat - oxy)

   _SET_SURFACE_EXCHANGE_(self%id_oxy,o2flux)

   _SET_SURFACE_EXCHANGE_(self%id_no3,self%surface_deposition_no3)
   _SET_SURFACE_EXCHANGE_(self%id_nh4,self%surface_deposition_nh4)
   _SET_SURFACE_EXCHANGE_(self%id_pho,self%surface_deposition_pho)
   _SET_SURFACE_EXCHANGE_(self%id_sil,self%surface_deposition_sil)

#if 0
   ! calculate cyanobacteria surface production
   if (S <= 10.0) then
     tbg = 1.0_rk/(1.0_rk + exp(BioC(29)*(BioC(30)-T)))
   else
     tbg = 0.0_rk
   end if
   
   blight=max(tanh(BioC(3)*par),0.)
   up_pho = pho/(BioC(25)+pho)
   prod = BioC(28) * bg * Tbg * min(blight, up_pho) ! cyanobacteria production

   !_SET_ODE_(self%id_bg,  prod)
   !_SET_ODE_(self%id_pho, -prod)
   !_SET_SURFACE_ODE_(id_oxy, ) ! not included in the modular ECOSMO version
   !_SET_SURFACE_ODE_(id_dic, -Prod) 
#endif

   ! Leave spatial loops over the horizontal domain (if any)
   _HORIZONTAL_LOOP_END_

   end subroutine do_surface
!EOC





!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Bottom fluxes for the ecosmo model
!
! !INTERFACE:

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
   class (type_hzg_ecosmo),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !LOCAL VARIABLES:
   real(rk) :: temp, tbs, oxy, no3, det, opa, sed1, sed2, sed3
   real(rk) :: pho, Rds, Rsd, Rsa, Rsdenit, Rsa_p, yt1, yt2
   real(rk) :: rhs, flux
   real(rk) :: bioom1, bioom2, bioom3, bioom4, bioom5, bioom6, bioom7, bioom8
!EOP
!-----------------------------------------------------------------------
!BOC
   _HORIZONTAL_LOOP_BEGIN_

   _GET_(self%id_temp,temp)
   _GET_(self%id_oxy,oxy)
   _GET_(self%id_det,det)
   _GET_(self%id_opa,opa)
   _GET_(self%id_no3,no3)
   _GET_HORIZONTAL_(self%id_sed1,sed1)
   _GET_HORIZONTAL_(self%id_sed2,sed2)
   _GET_HORIZONTAL_(self%id_sed3,sed3)
   _GET_HORIZONTAL_(self%id_tbs,tbs)

   bioom1 = 0.0_rk
   bioom2 = 0.0_rk
   bioom3 = 0.0_rk
   bioom4 = 0.0_rk
   bioom5 = 0.0_rk
   bioom6 = 0.0_rk
   bioom7 = 0.0_rk
   bioom8 = 0.0_rk
   if (oxy > 0) then
     bioom1 = 0.1/secs_pr_day * exp(temp*0.11_rk) * oxy/((0.1_rk*redf(7))+oxy)
     bioom2 = bioom1
     bioom6 = 1.0_rk
   else
     if (no3>0) then
       bioom5 = 5.0_rk
       bioom8 = 1.0_rk
     else
       bioom7 = 1.0_rk
     end if
   end if

!----citical bottom shear stress
        if (tbs.ge.BioC(34)) then
          Rsd=BioC(35)
          Rds=0.0_rk
        else if (tbs.lt.BioC(34)) then
          Rsd=0.0_rk
          Rds=BioC(36)
        end if

!---------------------------------------------------------------
!----denitrification parameter in dependence of available oxygen
        if (oxy .gt. 0.0) then
          Rsa=BioC(38)*exp(BioC(39)*temp)*1.0_rk
          Rsdenit=0.0_rk
        else if (oxy .le. 0.0) then
          Rsdenit=BioC(38)*exp(BioC(39)*temp)*2.0_rk
          Rsa=0.0_rk
        end if
                
        !--- sediment 1 total sediment biomass and nitrogen pool
        rhs = Rds*det - Rsd*sed1 - 2.0_rk*Rsa*sed1 - Rsdenit*sed1 - BioC(37)*sed1
        _SET_BOTTOM_ODE_(self%id_sed1, rhs)

        ! oxygen
        flux = -(BioOM6*6.625_rk*2.0_rk*Rsa*sed1 &
                 +BioOM7*6.625_rk*Rsdenit*sed1 &
                 +2.0_rk*BioOM1*Rsa*sed1) &
                *REDF(11)*REDF(16)
        _SET_BOTTOM_EXCHANGE_(self%id_oxy, flux)        

        ! nitrate
        _SET_BOTTOM_EXCHANGE_(self%id_no3, -BioOM5*Rsdenit*sed1)

        ! detritus
        _SET_BOTTOM_EXCHANGE_(self%id_det, Rsd*sed1 - Rds*det)

        ! ammonium
        _SET_BOTTOM_EXCHANGE_(self%id_nh4, (Rsdenit+Rsa)*sed1)


        !--try out for phosphate ute 2.6.2010        
        Rsa_p=BioC(38)*exp(BioC(39)*temp)*2.0_rk

        if (oxy.gt.0.0) then
          yt2=oxy/375.0_rk                  !normieren des wertes wie in Neumann et al 2002
          yt1=yt2**2.0_rk/(BioC(41)**2.0_rk+yt2**2.0_rk)
        
          _SET_BOTTOM_EXCHANGE_(self%id_pho, Rsa_p*(1.0_rk-BioC(40)*yt1)*sed3)

          !--sed 3 phosphate pool sediment+remineralization-P release 
          _SET_BOTTOM_ODE_(self%id_sed3, 2.0_rk*Rsa*sed1 - Rsa_p*(1.0_rk-BioC(40)*yt1)*sed3)

        else if (oxy.le.0.0) then
          _SET_BOTTOM_EXCHANGE_(self%id_pho, Rsa_p*sed3)
          _SET_BOTTOM_ODE_(self%id_sed3, Rsdenit*sed1 - Rsa_p*sed3)
        end if

        ! sediment opal(Si)
        _SET_BOTTOM_ODE_(self%id_sed2, Rds*opa - Rsd*sed2 - BioC(42)*sed2 - BioC(37)*sed2)
        _SET_BOTTOM_EXCHANGE_(self%id_opa, Rsd*sed2 - Rds*opa)
        _SET_BOTTOM_EXCHANGE_(self%id_sil, BioC(42)*sed2) 


   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom
!EOC


   subroutine get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
   class (type_hzg_ecosmo), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_

   real(rk)                     :: dom,det,dia,fla,bg

   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.

   _GET_(self%id_dia, dia)
   _GET_(self%id_fla, fla)
   _GET_(self%id_det, det)
   _GET_(self%id_dom, dom)
   _GET_(self%id_bg, bg)

   _SET_EXTINCTION_( self%BioC(5)*(dia+fla+bg) + self%BioC(4) )

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine get_light_extinction



   end module fabm_hzg_ecosmo


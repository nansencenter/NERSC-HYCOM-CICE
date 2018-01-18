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
   !real(rk)            :: PrmGI(2,6)=0.0_rk,GI(2,6)=0.0_rk
   !real(rk)            :: CHLtoC(3,2)=0.0_rk,ALFA(3)=0.0_rk
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model) :: type_hzg_ecosmo
!     Variable identifiers
      type (type_state_variable_id)         :: id_no3, id_nh4, id_pho, id_sil
      type (type_state_variable_id)         :: id_opa, id_det, id_dia, id_fla
      type (type_state_variable_id)         :: id_diachl, id_flachl, id_bgchl
      type (type_state_variable_id)         :: id_mesozoo, id_microzoo, id_bg, id_dom, id_oxy
      type (type_bottom_state_variable_id)  :: id_sed1, id_sed2, id_sed3
      type (type_dependency_id)             :: id_temp, id_salt, id_par
      type (type_dependency_id)             :: id_parmean
      type (type_horizontal_dependency_id)  :: id_tbs
      type (type_horizontal_dependency_id)  :: id_sfpar, id_meansfpar
      type (type_diagnostic_variable_id)    :: id_denit, id_primprod, id_secprod
      type (type_diagnostic_variable_id)    :: id_parmean_diag

      type (type_diagnostic_variable_id)    :: id_c2chl_fla, id_c2chl_dia, id_c2chl_bg

!     Model parameters
      real(rk) :: BioC(45)
      real(rk) :: zpr, frr
      real(rk) :: prefZsPs
      real(rk) :: prefZsPl
      real(rk) :: prefZsBG
      real(rk) :: prefZsD
      real(rk) :: prefZlPs
      real(rk) :: prefZlPl
      real(rk) :: prefZlBG
      real(rk) :: prefZlD
      real(rk) :: prefZlZs
      real(rk) :: surface_deposition_no3
      real(rk) :: surface_deposition_nh4
      real(rk) :: surface_deposition_pho
      real(rk) :: surface_deposition_sil
      real(rk) :: nfixation_minimum_daily_par
      real(rk) :: bg_growth_minimum_daily_rad
      real(rk) :: MAXchl2nPs, MINchl2nPs
      real(rk) :: MAXchl2nPl, MINchl2nPl
      real(rk) :: MAXchl2nBG, MINchl2nBG
      real(rk) :: alfaPl, alfaPs, alfaBG
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
!  Here, the ecosmo yaml is read and the variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_hzg_ecosmo),intent(inout),target  :: self
   integer,                intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister
!
!  Caglar Yumruktepe:
!  Added fabm.yaml support: parameters from yaml file are copied to
!                           BioC array. Eventually, BioC array maybe dropped
!                           from the model where parameter names from the yaml
!                           file will be used.
!  Added dynamic chlorophyll-a from Geider etal., 1997
!



! !LOCAL VARIABLES:
  ! real(rk) :: no3_init, nh4_init,pho_init,sil_init,oxy_init, &
  !             dia_init,fla_init,bg_init, &
  !             mesozoo_init,microzoo_init, &
  !             det_init, dom_init, opa_init, &
  !             sed1_init, sed2_init, sed3_init,&
  !             diachl_init,flachl_init,bgchl_init
   real(rk) :: sed1_init, sed2_init, sed3_init
   !real(rk) :: zpr=0.001_rk
   !real(rk) :: frr=0.4_rk
   !real(rk) :: surface_deposition_no3!=0.0_rk
   !real(rk) :: surface_deposition_nh4!=0.0_rk
   !real(rk) :: surface_deposition_pho!=0.0_rk
   !real(rk) :: surface_deposition_sil!=0.0_rk
   !real(rk) :: nfixation_minimum_daily_par!=40.0_rk
   !real(rk) :: bg_growth_minimum_daily_rad!=120.0_rk


   integer  :: i

!   namelist /hzg_ecosmo/  PrmBioC,GI,CHLtoC,ALFA,zpr,frr, &
!                          no3_init,nh4_init,pho_init, &
!                          sil_init,oxy_init, &
!                          dia_init,fla_init,bg_init, &
!                          diachl_init,flachl_init,bgchl_init, &
!                          mesozoo_init,microzoo_init, &
!                          det_init, dom_init, opa_init, &
!                          sed1_init, sed2_init, sed3_init, &
!                          surface_deposition_no3, surface_deposition_nh4, &
!                          surface_deposition_pho, surface_deposition_sil, &
!                          nfixation_minimum_daily_par, bg_growth_minimum_daily_rad

!EOP
!-----------------------------------------------------------------------
!BOC

   !include 'ECOSMparamNSBS.f'

   ! Read the namelist
   !if (configunit>0) read(configunit,nml=hzg_ecosmo,err=99,end=100)

   call self%get_parameter(self%zpr, 'zpr','zpr_long_name_needed', default=0.001_rk)
   call self%get_parameter(self%frr, 'frr','frr_long_name_needed', default=0.4_rk)
   call self%get_parameter(self%nfixation_minimum_daily_par, 'nfixation_minimum_daily_par', 'minimum daily PAR for N-fixation', default=40.0_rk) ! units needed
   call self%get_parameter(self%bg_growth_minimum_daily_rad, 'bg_growth_minimum_daily_rad', 'minimum daily rad for BG growth', default=120.0_rk) ! units needed

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
   call self%get_parameter( self%surface_deposition_no3, 'surface_deposition_no3', 'mmolN/m**2 d', 'surface deposition no3', default=0.0_rk, scale_factor=redf(1)*redf(6)/sedy0 )
   call self%get_parameter( self%surface_deposition_nh4, 'surface_deposition_nh4', 'mmolN/m**2 d', 'surface deposition nh4', default=0.0_rk, scale_factor=redf(1)*redf(6)/sedy0 )
   call self%get_parameter( self%surface_deposition_pho, 'surface_deposition_pho', 'mmolN/m**2 d', 'surface deposition pho', default=0.0_rk, scale_factor=redf(2)*redf(6)/sedy0 )
   call self%get_parameter( self%surface_deposition_sil, 'surface_deposition_sil', 'mmolN/m**2 d', 'surface deposition sil', default=0.0_rk, scale_factor=redf(3)*redf(6)/sedy0 )

   !  change units 1/day to 1/sec and mmolN,P,Si to mmolC
   call self%get_parameter( self%BioC(1) , 'muPl',        '1/day',      'max growth rate for Pl',          default=1.30_rk,  scale_factor=1.0_rk/sedy0             )
   call self%get_parameter( self%BioC(2) , 'muPs',        '1/day',      'max growth rate for Ps',          default=1.10_rk,  scale_factor=1.0_rk/sedy0             )
   call self%get_parameter( self%BioC(3) , 'aa',          'm**2/W',     'photosynthesis ef-cy',            default=0.04_rk                                         )
   call self%get_parameter( self%BioC(4) , 'EXw',         '1/m',        'light extinction',                default=0.041_rk                                        )
   call self%get_parameter( self%BioC(5) , 'Exphy',       'm**2/mmolN', 'phyto self-shading',              default=0.04_rk,  scale_factor=1.0_rk/(redf(1)*redf(6)) )
   call self%get_parameter( self%BioC(6) , 'rNH4',        'mmolN/m**3', 'NH4 half saturation',             default=0.20_rk,  scale_factor=redf(1)*redf(6)          )
   call self%get_parameter( self%BioC(7) , 'rNO3',        'mmolN/m**3', 'NO3 half saturation',             default=0.50_rk,  scale_factor=redf(1)*redf(6)          )
   call self%get_parameter( self%BioC(8) , 'psi',         'm**3/mmolN', 'NH4 inhibition',                  default=3.0_rk,   scale_factor=1.0_rk/(redf(1)*redf(6)) )
   call self%get_parameter( self%BioC(9) , 'mPl',         '1/day',      'Pl mortality rate',               default=0.04_rk,  scale_factor=1.0_rk/sedy0             )
   call self%get_parameter( self%BioC(10), 'mPs',         '1/day',      'Ps mortality rate',               default=0.08_rk,  scale_factor=1.0_rk/sedy0             )
   call self%get_parameter( self%BioC(11), 'GrZlP',       '1/day',      'Grazing rate Zl on Phyto',        default=0.80_rk,  scale_factor=1.0_rk/sedy0             )
   call self%get_parameter( self%BioC(12), 'GrZsP',       '1/day',      'Grazing rate Zs on Phyto',        default=1.00_rk,  scale_factor=1.0_rk/sedy0             )
   call self%get_parameter( self%BioC(13), 'GrZlZ',       '1/day',      'Grazing rate Zl on Zs',           default=0.50_rk,  scale_factor=1.0_rk/sedy0             )
   call self%get_parameter( self%BioC(14), 'Rg',          'mmolN/m**3', 'Zs, Zl half saturation',          default=0.50_rk,  scale_factor=redf(1)*redf(6)          )
   call self%get_parameter( self%BioC(15), 'mZl',         '1/day',      'Zl mortality rate',               default=0.10_rk,  scale_factor=1.0_rk/sedy0             )
   call self%get_parameter( self%BioC(16), 'mZs',         '1/day',      'Zs mortality rate',               default=0.20_rk,  scale_factor=1.0_rk/sedy0             )
   call self%get_parameter( self%BioC(17), 'excZl',       '1/day',      'Zl excretion rate',               default=0.06_rk,  scale_factor=1.0_rk/sedy0             )
   call self%get_parameter( self%BioC(18), 'excZs',       '1/day',      'Zs excretion rate',               default=0.08_rk,  scale_factor=1.0_rk/sedy0             )
   call self%get_parameter( self%BioC(19), 'gammaZlp',    '1',          'Zl assim. eff. on plankton',      default=0.75_rk                                         )
   call self%get_parameter( self%BioC(20), 'gammaZsp',    '1',          'Zs assim. eff. on plankton',      default=0.75_rk                                         )
   call self%get_parameter( self%BioC(21), 'gammaZd',     '1',          'Zl & Zs assim. eff. on det',      default=0.75_rk                                         )
   call self%get_parameter( self%BioC(22), 'reminD',      '1/day',      'Detritus remin. rate',            default=0.003_rk, scale_factor=1.0_rk/sedy0             )
   call self%get_parameter( self%BioC(23), 'sinkDet',     'm/day',      'Detritus sinking rate',           default=5.00_rk,  scale_factor=1.0_rk/sedy0             )
   call self%get_parameter( self%BioC(24), 'Wa',          'm/day',      '???',                             default=1.00_rk,  scale_factor=1.0_rk/sedy0             )
! set as mmolN/m**3 =0.25  !!!! need a correction
   call self%get_parameter( self%BioC(25),  'rPO4',       'mmolP/m**3', 'PO4 half saturation',             default=0.05_rk,  scale_factor=redf(2)*redf(6)          )
   call self%get_parameter( self%BioC(26),  'rSi',        'mmolSi/m**3','SiO2 half saturation',            default=0.50_rk,  scale_factor=redf(3)*redf(6)          )
   call self%get_parameter( self%BioC(27),  'regenSi',    '1/day',      'Si regeneration rate',            default=0.015_rk, scale_factor=1.0_rk/sedy0             )
   call self%get_parameter( self%BioC(28),  'muBG',       '1/day',      'max growth rate for BG',          default=1.00_rk,  scale_factor=1.0_rk/sedy0             )
   call self%get_parameter( self%BioC(29),  'TctrlBG',    '1/degC',     'BG T control beta',               default=1.00_rk                                         )
   call self%get_parameter( self%BioC(30),  'TrefBG',     'degC',       'BG reference temperature',        default=0.00_rk                                         )
   call self%get_parameter( self%BioC(31),  'GrBG',       '1/day',      'BG max grazing rate',             default=0.30_rk,  scale_factor=1.0_rk/sedy0             )
   call self%get_parameter( self%BioC(32),  'mBG',        '1/day',      'BG mortality rate',               default=0.08_rk,  scale_factor=1.0_rk/sedy0             )
   call self%get_parameter( self%BioC(33),  'upliftBG',   'm/day',      'BG uplifting rate',               default=0.10_rk,  scale_factor=1.0_rk/sedy0             )
   call self%get_parameter( self%BioC(34),  'crBotStr',   'N/m**2',     'critic. bot. stress for resusp.', default=0.007_rk                                        )
   call self%get_parameter( self%BioC(35),  'resuspRt',   '1/day',      'resuspension rate',               default=25.0_rk,  scale_factor=1.0_rk/sedy0             )
   call self%get_parameter( self%BioC(36),  'sedimRt',    'm/day',      'sedimentation rate',              default=3.5_rk,   scale_factor=1.0_rk/sedy0             )
   call self%get_parameter( self%BioC(37),  'burialRt',   '1/day',      'burial rate',                     default=1E-5_rk,  scale_factor=1.0_rk/sedy0             )
   call self%get_parameter( self%BioC(38),  'reminSED',   '1/day',      'sediment remineralization rate',  default=0.001_rk, scale_factor=1.0_rk/sedy0             )
   call self%get_parameter( self%BioC(39),  'TctrlDenit', '1/degC',     'temp. control denitrification',   default=0.15_rk                                         )
   call self%get_parameter( self%BioC(40),  'RelSEDp1',   'units??',    'P sedim. rel. p1',                default=0.15_rk                                         )
   call self%get_parameter( self%BioC(41),  'RelSEDp2',   'units??',    'P sedim. rel. p2',                default=0.10_rk                                         )
   call self%get_parameter( self%BioC(42),  'reminSEDsi', '1/day',      'sed. remineralization rate Si',   default=0.0002_rk,scale_factor=1.0_rk/sedy0             )
   call self%get_parameter( self%BioC(43),  'sinkOPAL',   'm/day',      'OPAL sinking rate',               default=5.0_rk,   scale_factor=1.0_rk/sedy0             )
   call self%get_parameter( self%BioC(44),  'sinkBG',     'm/day',      'BG sinking rate',                 default=-1.0_rk,  scale_factor=1.0_rk/sedy0             )
   call self%get_parameter( self%BioC(45),  'sinkDia',    'm/day',      'Diatom sinking rate',             default=0.0_rk,   scale_factor=1.0_rk/sedy0             )

   !! Store parameter values in our own derived type
   !do i=1,size(BioC)
    ! self%BioC(i) = BioC(i)
   !end do
   !  growth fractions
   call self%get_parameter( self%prefZsPs,  'prefZsPs',   '1',          'Grazing preference Zs on Ps',     default=0.70_rk                                         )
   call self%get_parameter( self%prefZsPl,  'prefZsPl',   '1',          'Grazing preference Zs on Pl',     default=0.25_rk                                         )
   call self%get_parameter( self%prefZsD,   'prefZsD',    '1',          'Grazing preference Zs on Det.',   default=0.00_rk                                         )
   call self%get_parameter( self%prefZsBG,  'prefZsBG',   '1',          'Grazing preference Zs on BG',     default=0.30_rk                                         )
   call self%get_parameter( self%prefZlPs,  'prefZlPs',   '1',          'Grazing preference Zl on Ps',     default=0.10_rk                                         )
   call self%get_parameter( self%prefZlPl,  'prefZlPl',   '1',          'Grazing preference Zl on Pl',     default=0.85_rk                                         )
   call self%get_parameter( self%prefZlZs,  'prefZlZs',   '1',          'Grazing preference Zl on Zs',     default=0.15_rk                                         )
   call self%get_parameter( self%prefZlD,   'prefZlD',    '1',          'Grazing preference Zl on Det.',   default=0.00_rk                                         )
   call self%get_parameter( self%prefZlBG,  'prefZlBG',   '1',          'Grazing preference Zl on BG',     default=0.30_rk                                         )
   ! chlorophyll-a constants
   call self%get_parameter( self%MINchl2nPs, 'MINchl2nPs', 'mgChl/mmolN', 'minimum Chl to N ratio Ps', default=0.50_rk, scale_factor=redf(11)*redf(16)             )
   call self%get_parameter( self%MAXchl2nPs, 'MAXchl2nPs', 'mgChl/mmolN', 'maximum Chl to N ratio Ps', default=3.83_rk, scale_factor=redf(11)*redf(16)             )
   call self%get_parameter( self%MINchl2nPl, 'MINchl2nPl', 'mgChl/mmolN', 'minimum Chl to N ratio Pl', default=0.50_rk, scale_factor=redf(11)*redf(16)             )
   call self%get_parameter( self%MAXchl2nPl, 'MAXchl2nPl', 'mgChl/mmolN', 'maximum Chl to N ratio Pl', default=2.94_rk, scale_factor=redf(11)*redf(16)             )
   call self%get_parameter( self%MINchl2nBG, 'MINchl2nBG', 'mgChl/mmolN', 'minimum Chl to N ratio BG', default=0.50_rk, scale_factor=redf(11)*redf(16)             )
   call self%get_parameter( self%MAXchl2nBG, 'MAXchl2nBG', 'mgChl/mmolN', 'maximum Chl to N ratio BG', default=3.83_rk, scale_factor=redf(11)*redf(16)             )
   call self%get_parameter( self%alfaPs,     'alfaPs', 'mmolN m2/(mgChl day W)**-1', 'initial slope P-I curve Ps', default=0.0393_rk, scale_factor=redf(1)*redf(6) )
   call self%get_parameter( self%alfaPl,     'alfaPl', 'mmolN m2/(mgChl day W)**-1', 'initial slope P-I curve Pl', default=0.0531_rk, scale_factor=redf(1)*redf(6) )
   call self%get_parameter( self%alfaBG,     'alfaBG', 'mmolN m2/(mgChl day W)**-1', 'initial slope P-I curve BG', default=0.0393_rk, scale_factor=redf(1)*redf(6) )


!write(*,*) 'chltoc',CHLtoC
!write(*,*) 'alfa',ALFA
!write(*,*) 'MINchl2nPs',self%MINchl2nPs
!write(*,*) 'MINchl2nPl',self%MINchl2nPl
!write(*,*) 'MAXchl2nPs',self%MAXchl2nPs
!write(*,*) 'MAXchl2nPl',self%MAXchl2nPl
!write(*,*) 'alfa_dia',self%alfaPl
!write(*,*) 'alfa_fla',self%alfaPs

   ! Register state variables
   call self%register_state_variable( self%id_no3,      'no3',    'mgC/m3',    'nitrate',                   minimum=0.0_rk,        vertical_movement=0.0_rk )
   call self%register_state_variable( self%id_nh4,      'nh4',    'mgC/m3',    'ammonium',                  minimum=0.0_rk,        vertical_movement=0.0_rk )
   call self%register_state_variable( self%id_pho,      'pho',    'mgC/m3',    'phosphate',                 minimum=0.0_rk,        vertical_movement=0.0_rk )
   call self%register_state_variable( self%id_sil,      'sil',    'mgC/m3',    'silicate',                  minimum=0.0_rk,        vertical_movement=0.0_rk )
   call self%register_state_variable( self%id_oxy,      'oxy',    'mmolO2/m3', 'oxygen',                    minimum=-100000._rk,   vertical_movement=0.0_rk )
   call self%register_state_variable( self%id_fla,      'fla',    'mgC/m3',    'small phytoplankton',       minimum=1.0e-7_rk,     vertical_movement=0.0_rk )
   call self%register_state_variable( self%id_dia,      'dia',    'mgC/m3',    'large phytoplankton',       minimum=1.0e-7_rk,     vertical_movement=-self%BioC(45) )
   call self%register_state_variable( self%id_bg,       'bg',     'mgC/m3',    'cyanobacteria',             minimum=1.0e-7_rk,     vertical_movement=-self%BioC(44) )
   call self%register_state_variable( self%id_diachl,   'diachl', 'mgChl/m3',  'large phytoplankton chl-a', minimum=1.0e-7_rk/27., vertical_movement=-self%BioC(45) )
   call self%register_state_variable( self%id_flachl,   'flachl', 'mgChl/m3',  'small phytoplankton chl-a', minimum=1.0e-7_rk/20., vertical_movement=0.0_rk )
   call self%register_state_variable( self%id_bgchl,    'bgchl',  'mgChl/m3',  'cyanobacteria chl-a',       minimum=1.0e-7_rk/20., vertical_movement=-self%BioC(44) )

   call self%register_state_variable( self%id_microzoo, 'microzoo','mgC/m3','microzooplankton',     &
                                     vertical_movement=0.0_rk, &
                                     minimum=1.0e-7_rk)

   call self%register_state_variable(self%id_mesozoo,'mesozoo','mgC/m3','mesozooplankton',     &
                                     vertical_movement=0.0_rk, &
                                     minimum=1.0e-7_rk)

   call self%register_state_variable(self%id_det,'det','mgC/m3','detritus',     &
                                     vertical_movement=-self%BioC(23), &
                                     minimum=0.0_rk)

   call self%register_state_variable(self%id_opa,'opa','mgC/m3','opal',     &
                                     vertical_movement=-self%BioC(43), &
                                     minimum=0.0_rk)

   call self%register_state_variable(self%id_dom,'dom','mgC/m3','labile dissolved om',     &
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

   ! calls outputs - simulated Carbon to chlorophyll-a ratio
   call self%register_diagnostic_variable(self%id_c2chl_fla,'c2chl_fla','mgC/mgCHL', &
         'daily-mean C to CHL ratio for flagellates', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_c2chl_dia,'c2chl_dia','mgC/mgCHL', &
         'daily-mean C to CHL ratio for diatoms', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_c2chl_bg,'c2chl_bg','mgC/mgCHL', &
         'daily-mean C to CHL ratio for cyanobacteria', output=output_time_step_averaged)

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

!99 call self%fatal_error('hzg_ecosmo_initialize','Error reading namelist hzg_ecosmo.')

!100 call self%fatal_error('hzg_ecosmo_initialize','Namelist hzg_ecosmo was not found.')

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
   real(rk) :: flachl,diachl,bgchl,chl2c_fla,chl2c_dia,chl2c_bg
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
   real(rk) :: fla_loss=1.0_rk
   real(rk) :: dia_loss=1.0_rk
   real(rk) :: bg_loss=1.0_rk
   real(rk) :: mic_loss=1.0_rk
   real(rk) :: mes_loss=1.0_rk
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
   _GET_(self%id_diachl,diachl)
   _GET_(self%id_flachl,flachl)
   _GET_(self%id_bgchl,bgchl)
   _GET_(self%id_microzoo,microzoo)
   _GET_(self%id_mesozoo,mesozoo)
   _GET_(self%id_det,det)
   _GET_(self%id_dom,dom)
   _GET_(self%id_opa,opa)
   _GET_(self%id_oxy,oxy)
   _GET_(self%id_parmean,mean_par)
   _GET_HORIZONTAL_(self%id_meansfpar,mean_surface_par)

   ! CAGLAR
   ! checks - whether the biomass of plankton is below a predefined threshold,
   !          where below the threshold, loss terms are removed from the RHS of
   !          the equations. The idea is to keep plankton safe from extinction.
   ! loss terms are multiplied by the constants below, which can only be set
   ! by the model to 0 or 1.

   fla_loss = max(sign(-1.0_rk,fla-0.01_rk),0.0_rk)      ! flagellates
   dia_loss = max(sign(-1.0_rk,dia-0.01_rk),0.0_rk)      ! diatoms
   bg_loss  = max(sign(-1.0_rk,bg-0.01_rk),0.0_rk)       ! cyanobacteria
   mic_loss = max(sign(-1.0_rk,microzoo-0.001_rk),0.0_rk) ! microzooplankton
   mes_loss = max(sign(-1.0_rk,microzoo-0.001_rk),0.0_rk) ! mesozooplankton


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

   ! chlorophyll-a to C change
   chl2c_fla = self%MAXchl2nPs * max(0.1,Ps_prod) * BioC(2) * sedy0 * fla / &
               (self%alfaPs * par * flachl)
   chl2c_dia = self%MAXchl2nPl * max(0.1,Pl_prod) * BioC(1) * sedy0 * dia / &
               (self%alfaPl * par * diachl)
   chl2c_bg = self%MAXchl2nBG * max(0.1,Bg_prod) * BioC(28) * sedy0 * bg / &
               (self%alfaBG * par * bgchl)

            chl2c_fla = max(self%MINchl2nPs,chl2c_fla)
            chl2c_fla = min(self%MAXchl2nPs,chl2c_fla)
            chl2c_dia = max(self%MINchl2nPl,chl2c_dia)
            chl2c_dia = min(self%MAXchl2nPl,chl2c_dia)
            chl2c_bg = max(self%MINchl2nBG,chl2c_bg)
            chl2c_bg = min(self%MAXchl2nBG,chl2c_bg)

   !         chl2c_fla = 1./20.
   !         chl2c_dia = 1./27.

   ! grazing
   ! gf denotes grazing fraction

   Fs = self%prefZsPs*fla + self%prefZsPl*dia + self%prefZsD*det + self%prefZsBG*bg
   Fl = self%prefZlPs*fla + self%prefZlPl*dia + self%prefZlZs*microzoo + &
          self%prefZlD*det + self%prefZlBG*bg

   ZsonPs = fla_loss * BioC(12) * self%prefZsPs * fla/(BioC(14) + Fs)
   ZsonPl = dia_loss * BioC(12) * self%prefZsPl * dia/(BioC(14) + Fs)
   ZsonD  =            BioC(12) * self%prefZsD * det/(BioC(14) + Fs)
   ZsonBg = bg_loss  * BioC(31) * self%prefZsBG * bg/(BioC(14) + Fs)

   ZlonPs = fla_loss * BioC(11) * self%prefZlPs * fla/(BioC(14) + Fl)
   ZlonPl = dia_loss * BioC(11) * self%prefZlPl * dia/(BioC(14) + Fl)
   ZlonD =             BioC(11) * self%prefZlD * det/(BioC(14) + Fl)
   ZlonZs = mic_loss * BioC(13) * self%prefZlZs * microzoo/(BioC(14) + Fl)
   ZlonBg = bg_loss  * BioC(31) * self%prefZlBG * bg/(BioC(14) + Fl)

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

   _SET_ODE_(self%id_fla, (BioC(2)*Ps_prod - BioC(10)*fla_loss)*fla - ZsonPs*microzoo - ZlonPs*mesozoo)
   _SET_ODE_(self%id_dia, (BioC(1)*Pl_prod - BioC(9)*dia_loss)*dia - ZsonPl*microzoo - ZlonPl*mesozoo)
   _SET_ODE_(self%id_bg,  (BioC(28)*(Bg_prod + Bg_fix) - BioC(32)*bg_loss)*bg - ZsonBg*microzoo - ZlonBg*mesozoo)

  ! for chlorophyll-a

   rhs = BioC(2)*Ps_prod*chl2c_fla*fla - ( (BioC(10)*fla_loss*fla + ZsonPs*microzoo + ZlonPs*mesozoo)*flachl/fla )
   _SET_ODE_(self%id_flachl,rhs)
   rhs = BioC(1)*Pl_prod*chl2c_dia*dia - ( (BioC(9)*dia*dia_loss + ZsonPl*microzoo + ZlonPl*mesozoo)*diachl/dia )
   _SET_ODE_(self%id_diachl,rhs)
   rhs = BioC(28)*(Bg_prod + Bg_fix)*chl2c_bg*bg - ( (BioC(32)*bg*bg_loss + ZsonBg*microzoo + ZlonBg*mesozoo)*bgchl/bg )
   _SET_ODE_(self%id_bgchl,rhs)


   ! microzooplankton

   Zs_prod = BioC(20)*(ZsonPs + ZsonPl + ZsonBg) + BioC(21)*ZsonD
   rhs = (Zs_prod - (BioC(16) + BioC(18) + self%zpr)*mic_loss) * microzoo &
         - ZlonZs * mesozoo
   _SET_ODE_(self%id_microzoo, rhs)

   ! mesozooplankton
   Zl_prod = BioC(19)*(ZlonPs + ZlonPl + ZlonBg + ZlonZs) + BioC(21)*ZlonD
   rhs = (Zl_prod - (BioC(15) + BioC(17) + self%zpr)*mes_loss) * mesozoo
   _SET_ODE_(self%id_mesozoo, rhs)

   ! detritus
   dxxdet = (  ((1.0_rk-BioC(20))*(ZsonPs + ZsonPl + ZsonBg) &
              + (1.0_rk-BioC(21))*ZsonD) * microzoo &
              + ((1.0_rk-BioC(19))*(ZlonPs + ZlonPl + ZlonBg + ZlonZs) &
              + (1.0_rk-BioC(21))*ZlonD) * mesozoo &
              + BioC(16) * microzoo * mic_loss &
              + BioC(15) * mesozoo * mes_loss &
              + BioC(10) * fla * fla_loss &
              + BioC(9)  * dia * dia_loss &
              + BioC(32) * bg * bg_loss )

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
         + BioC(18) * microzoo * mic_loss &
         + BioC(17) * mesozoo * mes_loss &
         + frem * det &
         + fremDOM * dom - bioom1 * nh4
   _SET_ODE_(self%id_nh4, rhs)

   ! phosphate

   rhs = -Prod -BioC(28)*bg*Bg_fix &
         + BioC(18) * microzoo * mic_loss &
         + BioC(17) * mesozoo * mes_loss &
         + frem*det + fremDOM*dom
   _SET_ODE_(self%id_pho, rhs)


   ! silicate
   _SET_ODE_(self%id_sil, -BioC(1)*Pl_prod*dia + BioC(27)*opa)

   ! opal
   _SET_ODE_(self%id_opa, BioC(9)*dia*dia_loss + ZsonPl*microzoo + ZlonPl*mesozoo - BioC(27)*opa)

   ! oxygen
   rhs = ((6.625*up_nh4 + 8.125*up_no3+1.d-10)/(up_n+1.d-10)*Prod &
         -bioom6*6.625*(BioC(18)*microzoo*mic_loss + BioC(17)*mesozoo*mes_loss) &
         -frem*det*(bioom6+bioom7)*6.625 &
         -(bioom6+bioom7)*6.625*fremDOM*dom &
         -2.0_rk*bioom1*nh4)*redf(11)*redf(16)
   _SET_ODE_(self%id_oxy, rhs)

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_denit,(frem*det*bioom5+fremDOM*dom*bioom5)*redf(11)*redf(16))
   _SET_DIAGNOSTIC_(self%id_primprod, Prod + BioC(28)*bg*Bg_fix )
   _SET_DIAGNOSTIC_(self%id_secprod, Zl_prod*mesozoo + Zs_prod*microzoo)
   _SET_DIAGNOSTIC_(self%id_parmean_diag, mean_par)

   _SET_DIAGNOSTIC_(self%id_c2chl_fla, 1.0_rk/chl2c_fla)
   _SET_DIAGNOSTIC_(self%id_c2chl_dia, 1.0_rk/chl2c_dia)
   _SET_DIAGNOSTIC_(self%id_c2chl_bg, 1.0_rk/chl2c_bg)
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

   real(rk)                     :: dom,det,diachl,flachl,bgchl

   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.

   _GET_(self%id_diachl, diachl)
   _GET_(self%id_flachl, flachl)
   _GET_(self%id_det, det)
   _GET_(self%id_dom, dom)
   _GET_(self%id_bgchl, bgchl)

   _SET_EXTINCTION_( self%BioC(5)*(diachl+flachl+bgchl) + self%BioC(4) )

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine get_light_extinction



   end module fabm_hzg_ecosmo

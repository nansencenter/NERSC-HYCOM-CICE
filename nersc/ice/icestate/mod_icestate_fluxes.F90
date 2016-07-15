! ======================================================================
!  This module contains all fluxes passed between the icestate modules
!  and other ice/ocean model components
! ======================================================================
!Written by: David Salas y Melia  and  Ragnhild Hansen
!Last changed 03.01.2002 by Knut Arild Liseter
!Changed 24.09.2004 to set up for MPI

MODULE mod_icestate_fluxes
   use mod_icestate , only : nthick
   use mod_xc
   implicit none

   real, parameter :: windi=0.006 ! Minimum friction velocity under ice


   !  - Fluxes to MICOM
   real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::  &
      Isalflx,Isalflx2,Isurflx,Ibuoyfl



   !  - These fluxes are only used in HYCOM
   real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
      Ibuoylw, Ibuoyt, Ibuoysw, Isswflx



   !  - Stress going into the ocean. ustar is used in MICOM and HYCOM
   real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
      Iustar


   !  - Divergence of the ice field. Used in mechanical redistribution
   real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
      divu


   !  - The atmospheric forcing - Kill this and use mod_forcing ?
   real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),save :: &
      Iwndspd, Irelhum, Iairtmp, Iprecip,  &
      Iclouds, Islp, Itaux, Itauy, Itauxice, Itauyice, &
      Iradflx, ivapmix, Iswflx
    
   real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),save :: cawdir_day
   real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),save :: radfl_day

   ! Diagnostic fluxes
   real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,nthick),save ::          &
     icestate_swfl, icestate_swtr, icestate_trb, icestate_lw, icestate_brfl, &
     icestate_mlfl, icestate_ctop, icestate_cbot, icestate_ntop,             &
     icestate_nbot, icestate_grw, icestate_lgrw
   real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),save ::          &
     icestate_lead_tot, icestate_lead_sw, icestate_lead_lw, icestate_lead_trb

end module mod_icestate_fluxes

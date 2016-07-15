      module mod_icestate_diag
      contains

! ===================================================================================
! =========================== icestate_area =========================================
! ===================================================================================
! Routine calculates the area of ocean covered by ice

      real*8 function icestate_area()
      use mod_xc
      use mod_icestate  , only : icestate
      implicit none
      integer i,j,l
      real icearea(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      include 'common_blocks.h'

C$OMP PARALLEL DO PRIVATE(j,l,i)
C$OMP&SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do l=1,isp(j)
      do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
         icearea(i,j)=  sum(icestate(i,j)%ice%fice)*scpx(i,j)*scpy(i,j)
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO
      call xcsum(icestate_area,icearea,ip)
      end function icestate_area


! ===================================================================================
! =========================== icestate_extent =======================================
! ===================================================================================
! Routine calculates the area of ocean with an ice concentration > treshold

      real*8 function icestate_extent(treshold)
      use mod_xc
      use mod_icestate , only : icestate
      use mod_icestate_tools
      implicit none
      integer i,j,l
      real, intent(in) :: treshold
      real iceextent(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      include 'common_blocks.h'

C$OMP PARALLEL DO PRIVATE(j,l,i) 
C$OMP&SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do l=1,isp(j)
      do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
         iceextent(i,j)= scpx(i,j)*scpy(i,j)*
     &        step(sum(icestate(i,j)%ice%fice),treshold)
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO
      call xcsum(icestate_extent,iceextent,ip)
      end function icestate_extent


! ===================================================================================
! =========================== icestate_volume========================================
! ===================================================================================
! Routine calculates the volume of ice

      real*8 function icestate_volume()
      use mod_xc
      use mod_icestate , only : icestate
      implicit none
      integer i,j,l
      real volume(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      include 'common_blocks.h'

C$OMP PARALLEL DO PRIVATE(j,l,i) 
C$OMP&SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do l=1,isp(j)
      do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
         volume(i,j)=
     &              sum(icestate(i,j)%ice%fice*icestate(i,j)%ice%hice)*
     &              scpx(i,j)*scpy(i,j)
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO
      call xcsum(icestate_volume,volume,ip)
      end function icestate_volume

! ===================================================================================
! =========================== icestate_salfl2 =======================================
! ===================================================================================

      real*8 function icestate_salfl2(dts)
      use mod_xc
      use mod_icestate , only: icestate, rhoref
      use mod_icestate_fluxes
      implicit none
      integer i,j,l
      real,intent(in) ::dts
      real salfl(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      include 'common_blocks.h'

C$OMP PARALLEL DO PRIVATE(j,l,i) 
C$OMP&SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do l=1,isp(j)
      do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
         salfl(i,j)=
     &      ((Isalflx(i,j)-Isalflx2(i,j))*scpx(i,j)*scpy(i,j))/ 
     &      (rhoref*depths(i,j))
      enddo
      enddo
      enddo
C$OMP END PARALLEL DO
      call xcsum(icestate_salfl2,salfl,ip)
      end function icestate_salfl2






! ===================================================================================
! ============================== Energy_difference ==================================
! ===================================================================================
! Routine to get thermo energy changes in ice and ocean (atm heat flux)
!
      real function energy_difference(old,new,frac)
      use mod_icestate
      use mod_icestate_tools
      implicit none
      type (t_istate_cell), intent(in) :: old,new
      real,                 intent(in) :: frac

      real :: ice_delta_e, ml_delta_e,rhosw_old,rhosw_new,delta_q
      !real :: qthref 
      !include 'stmt_fns.h'
      !qthref=1./thref ! For stmt_fns

      ! Get new an old enthalpies
      ice_delta_e = sum(ice_enthalpy(new%ice)) - 
     &              sum(ice_enthalpy(old%ice))

      ! Heat stored in ice
      delta_q = sum( new%ice%fice*new%ice%qstore  - 
     &               old%ice%fice*old%ice%qstore )

      ! Get new and old energy content in mixed layer
      rhosw_old = 1000. + sig0(old%tml-t0deg,old%sml)
      rhosw_new = 1000. + sig0(new%tml-t0deg,new%sml)
      ml_delta_e  = cpsw*(rhosw_new*new%hml*new%tml - 
     &                    rhosw_old*old%hml*old%tml)

      !print *,ice_delta_e,ml_delta_e*frac
      energy_difference = ice_delta_e + ml_delta_e*frac + delta_q
      

      end function energy_difference






! ============================== Ice_enthalpy  ======================================
! Returns enthalpy of the ice. (Heat needed to heat ice to melting pt and melt it)
!
      elemental function ice_enthalpy(icem)
      use mod_icestate 
      implicit none
      type(t_ice), intent(in) :: icem
      real :: ice_enthalpy
      real :: lath,spch,hofusn

      ! Effective latent heat of fusion
      hofusn=hofusn0
      if (icem%qstore>epsil1.and.icem%hice>epsil1) then
         hofusn = hofusn0 - icem%qstore/icem%hice
      end if

      
      ! Sum up latent heat in ice/snow 
      lath = icem%hsnw*icem%fice*icem%rhosnw/rhoice*hofusn0
      lath = lath + hofusn*icem%hice*icem%fice

      ! Sum up specific heat in ice (snow not accounted for)
      spch =    cpice * rhoice * icem%hice * icem%fice 
     &          * sum(tice_m-icem%vtp(1:icem%nlay)) / icem%nlay

      ! Enthalpy is energy needed to heat ice to the melting point and melt it
      ice_enthalpy = lath + spch
      end function ice_enthalpy



! ===================================================================================
! ============================== Salinity_difference ================================
! ===================================================================================
! Routine to get difference in salinity between old and new state (should be zero)

      real function salinity_difference(old,new,frac)
      use mod_icestate
      use mod_icestate_tools
      implicit none
      type (t_istate_cell), intent(in) :: old,new
      real, intent(in) :: frac
      real :: ice_ds,ml_ds,rhosw_old,rhosw_new
      !real :: qthref
      !include 'stmt_fns.h'
      !qthref=1./thref ! For stmt_fns

      ! Get new an old enthalpies
      ice_ds = rhoice * sice * 
     &         (sum(new%ice%hice*new%ice%fice) - 
     &          sum(old%ice%hice*old%ice%fice))

      ! Get new and old salinity content in mixed layer
      rhosw_old = 1000. + sig0(old%tml-t0deg,old%sml)
      rhosw_new = 1000. + sig0(new%tml-t0deg,new%sml)
      ml_ds  =  ( rhosw_new*new%hml*new%sml - rhosw_old*old%hml*old%sml)

      salinity_difference = ice_ds + ml_ds*frac

      ! Relative sal difference :
      salinity_difference = salinity_difference 

      !print *,ml_ds,ice_ds

      end function salinity_difference




! ===================================================================================
! ============================== mass_difference ====================================
! ===================================================================================
! Routine to get difference in mass between old and new state (evap)

      real function mass_difference(old,new,frac)
      use mod_icestate
      use mod_icestate_tools
      implicit none
      type (t_istate_cell), intent(in) :: old,new
      real, intent(in) :: frac
      real :: ice_dm,ml_dm,rhosw_old,rhosw_new,snw_dm
      !real :: qthref
      !include 'stmt_fns.h'
      !qthref=1./thref ! For stmt_fns


      ! Get new an old enthalpies
      ice_dm = rhoice * sum(new%ice%hice*new%ice%fice - old%ice%hice*
     &   old%ice%fice)
      snw_dm = sum(new%ice%hsnw*new%ice%fice*new%ice%rhosnw  
     &           - old%ice%hsnw*old%ice%fice*old%ice%rhosnw)

      ! Get new and old mass of  "mixed" layer
      rhosw_old = 1000. + sig0(old%tml-t0deg,old%sml)
      rhosw_new = 1000. + sig0(new%tml-t0deg,new%sml)
      ml_dm  =  new%hml*rhosw_new - old%hml*rhosw_old

      mass_difference = ice_dm + frac*ml_dm + snw_dm

      !print *,ice_dm,ml_dm,snw_dm

      end function mass_difference


! ===================================================================================
! ============================== show_state =========================================
! ===================================================================================
! Subroutine 

      subroutine show_state(state)
      use mod_icestate
      implicit none
      
      type(t_istate_cell), intent(in) :: state

      ! Set output format
      character (len= 2)      :: cnn,cll
      character (len=20)      :: frmt,frmt2,frmt3

      write (cnn,'(i2.2)') nthick
      frmt= cnn//'f15.7'
      write (cll,'(i2.2)') nlaymax
      frmt2= trim(cll)//'f8.3'


      print '(a)',''
      print '(a)','Mixed Layer'
      print '(a)','==========='
      print '(a,f12.6)','Temperature:',state%tml
      print '(a,f12.6)','Salinity   :',state%sml
      print '(a,f12.6)','Depth      :',state%hml

      call show_ice(state%ice)

      end subroutine show_state


      subroutine show_ice(icem)
      use mod_icestate
      implicit none
      type(t_ice), intent(in) :: icem(nthick)

      ! Set output format
      character (len= 2)      :: cnn,cll
      character (len=20)      :: frmt,frmt2,frmt3

      write (cnn,'(i2.2)') nthick
      frmt= cnn//'f15.7'
      write (cll,'(i2.2)') nlaymax
      frmt2= trim(cll)//'f8.3'
      

      print '(a)',''
      print '(a)','ICE'
      print '(a)','==='
      print '(a,'//trim(frmt)//')','Fraction   :',icem%fice
      print '(a,'//trim(frmt)//')','Thickness  :',icem%hice
      print '(a,'//trim(frmt)//')','Snow thk   :',icem%hsnw
      print '(a,'//trim(frmt)//')','Snow dens  :',icem%rhosnw
      print '(a,'//trim(frmt)//')','Snow alb   :',icem%albs
      print '(a,'//trim(frmt)//')','Surf temp  :',icem%tsrf
      end subroutine show_ice


      end module mod_icestate_diag

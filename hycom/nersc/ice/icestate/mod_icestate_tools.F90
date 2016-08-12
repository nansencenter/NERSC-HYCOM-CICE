! =======================================================================
! ============================ MODULE TOOLS ============================= 
! =======================================================================
! This module contains various .. tools ...
MODULE mod_icestate_tools
  IMPLICIT NONE 

  INTERFACE step
     module procedure step_scal_scal
     module procedure step_scal_vec
     module procedure step_vec_vec
  end INTERFACE



 interface operator (+)
      module procedure add_ice
   end interface

   interface clear_ice
      module procedure clear_ice_scal
      module procedure clear_ice_vec
   end interface clear_ice
#if defined (SSNOWD)
   interface snow_frac
      module procedure snow_frac_scal
      module procedure snow_frac_vec
   end interface snow_frac

   interface average_depth
      module procedure average_depth_scal
      module procedure average_depth_vec
   end interface average_depth
#endif
   real, private, parameter :: &
      c1=-1.36471E-01, c2= 4.68181E-02, c3= 8.07004E-01, &
      c4=-7.45353E-03, c5=-2.94418E-03, c6= 3.43570E-05, &
      c7= 3.48658E-05


 CONTAINS 


! ======================================================================
! ================== Freezing point temperature=========================
! ======================================================================
   pure real function freeze_temp(saln)
      implicit none
      real, intent(in) :: saln
      freeze_temp = 273.216 - .057 * saln
   end function freeze_temp


! ======================================================================
! ================== Sigma 0 functions from hycom ======================
! ======================================================================
   pure real function sig0(t,s)
      implicit none
      real, intent(in) :: t,s
      sig0=(c1+c3*s+t*(c2+c5*s+t*(c4+c7*s+c6*t)))
   end function sig0

   pure real function dsig0ds(t,s)
      real, intent(in) :: t,s
      dsig0ds=(c3+t*(c5+t*c7))
   end function dsig0ds

   pure real function dsig0dt(t,s)
      real, intent(in) :: t,s
      dsig0dt=(c2+c5*s+2.*t*(c4+c7*s+1.5*c6*t))
   end function dsig0dt



! =======================================================================
! ========================== step functions =============================
! =======================================================================
!
! The step function is set to 1 if the first parameter (x) is bigger 
! than the second parameter (y), otherwise zero. Below are three 
! "versions". Note that none of these should be called directly, they are 
! called "step", and the interface defined in the start of the module 
! takes care of the rest.
! 
! =======================================================================
! Both inputs are scalar values in this version. The result is scalar.
   PURE FUNCTION step_scal_scal(x,y)
     implicit none
     REAL, INTENT(in) :: x,y
     REAL :: step_scal_scal
     step_scal_scal= (1. + SIGN(1.,x-y)) * .5
   END FUNCTION step_scal_scal
! =======================================================================
! First input is a vector, second input is a scalar in this version. The
! result is a vector.
   FUNCTION step_scal_vec(x,y)
     implicit none
     REAL, INTENT(in) :: x(:),y
     REAL :: step_scal_vec(size(x))
     integer i
     forall (i=1:size(x))
        step_scal_vec(i) = (1. + SIGN(1.,x(i)-y)) * .5
     end forall
   END FUNCTION step_scal_vec
! =======================================================================
! Both inputs are vectors in this version. The result is a vector.
   FUNCTION step_vec_vec(x,y)
     implicit none
     REAL, INTENT(in) :: x(:),y(:)
     REAL :: step_vec_vec(size(x))
     step_vec_vec(:) = (1. + SIGN(1.,x(:)-y(:))) * .5
   END FUNCTION step_vec_vec


















! -----------------------------------------------------------------------
! ------------------------- FUNCTION VTPLIN -----------------------------
! -----------------------------------------------------------------------
pure function VtpLin(icem,tsrf,tbot)
  ! Set up a linear temperature profile in  the ice layers.
  USE mod_icestate
  IMPLICIT none
  REAL,         INTENT(in)  :: tsrf, tbot
  type (t_ice), intent(in)  :: icem
  type (t_ice)              :: VtpLin

  real    :: rkeff, a, bs, bi, dz,z0,rksn
  integer :: hk,hl

  VtpLin=icem

  ! Get snow conductivity
  rksn  = rkice * (icem%rhosnw / rhow)**1.885

  ! Get effective conductivity
#if defined(SSNOWD) && defined (HEAT_CONDUC)
  VtpLin%hsnw = average_depth(VtpLin)
  rkeff = corr_factor(VtpLin)*rksn*rkice/(rkice*VtpLin%hsnw + rksn*VtpLin%hice)
#else 
  rkeff = rksn*rkice/(rkice*icem%hsnw + rksn*icem%hice)
#endif  

  ! Proportionality constant in ice  and constant value.
  a     = ( tsrf - tbot ) * rkeff
  bi    = tbot

  ! Spacing increment
  dz    = icem%hice/icem%nlay

  ! Position of first layer midpoint
  z0    =  dz/2

  ! Regardless of snow or not the following gives the temp. profile:
  do hl=1,icem%nlay

     VtpLin%vtp(hl) = a*((hl-1)*dz + z0)/rkice + bi 
  end do

  ! Set remaining layers  to uppermost layer temp.
  VtpLin%vtp((icem%nlay+1):nlaymax) = VtpLin%vtp(icem%nlay)

end function VtpLin

! -----------------------------------------------------------------------
! ------------------------- FUNCTION   VTP_LIN --------------------------
! -----------------------------------------------------------------------
PURE function vtp_lin(tice_fx,ficex,hicex,hsnwx,rhosnwx,rksnwx,tsrfx,vtpx,numlay)
  ! Set up a linear temperature in each of the ice/snow layers. This does
  ! not mean that the overal TP is linear ... (different cond in snow
  ! and ice). 
  use mod_icestate
  IMPLICIT none
  REAL,                       INTENT(in) :: tice_fx 
  REAL,                       INTENT(in) :: ficex,hicex,hsnwx,rhosnwx,rksnwx,tsrfx
  REAL, DIMENSION(1:nlaymax), INTENT(in) :: vtpx
  integer,                    intent(in) :: numlay
  real, dimension(1:nlaymax)             :: vtp_lin
  real    :: rkeff, a, bs, bi, dz,z0
  integer :: hk,hl

  ! Get effective conductivity 
  rkeff = rksnwx*rkice/(rkice*hsnwx + rksnwx*hicex)

  ! Proportionality constant in ice and snow ( modified by conductivity )
  a     = ( tsrfx - tice_fx ) * rkeff

  ! ... and constant values for ice and snow.
  bs    = tsrfx - a * (hicex + hsnwx) / rksnwx
  bi    = tice_fx

  ! Spacing increment
  dz    = hicex/numlay

  ! Position of first layer midpoint
  z0    = dz/2

  ! Regardless of snow or not the following gives the temp. profile:
  do hl=1,numlay

     vtp_lin(hl) = a*((hl-1)*dz + z0)/rkice + bi 
  end do


  ! Set remaining layers  to uppermost layer temp.
  vtp_lin(numlay+1:nlaymax) = vtp_lin(numlay)

end function vtp_lin






! -----------------------------------------------------------------------
! --------------------- FUNCTION CHANGECLASS ----------------------------
! -----------------------------------------------------------------------
! If ice initially in one class through some process (ice growth/melt, 
! advection, redistribution) moves into another class, this routine should
! be called. It ensures that the temperature profile with the old number of
! layers is represented with the new number of layers. This is done through
! conservation of energy.
pure function changeclass(icem,newnlay,tice_fXXX)
   use mod_icestate
   implicit none
   integer,     intent(in)  :: newnlay
   type(t_ice), intent(in)  :: icem
   type(t_ice)              :: changeclass
   real,        intent(in)  :: tice_fXXX
   integer :: i,oind_lo, oind_hi, j, oldnlay, curr_old_lay
   real    :: oind, dh_old, e, tmpdh, T(newnlay), dh, up, dn, zup, zdn, dheff
   real    :: tbot(nlaymax), Tgrad(nlaymax), rksn, rkeff

   changeclass = icem

   ! Old number of layers
   oldnlay = changeclass%nlay

   ! Old vertical increment
   dh_old = icem%hice/changeclass%nlay

   ! New vertical increment
   dh     = icem%hice/newnlay

   if ((oldnlay.ne.newnlay).and.icem%fice>epsil1 ) then 

      if (oldnlay>1) then

         do j = 1,oldnlay-1
            Tgrad(j)  = (icem%vtp(j+1) - icem%vtp(j))/dh_old
            tbot(j) = icem%vtp(j) - Tgrad(j)*dh_old/2
         end do
         Tgrad(oldnlay)  = Tgrad(oldnlay-1)
         tbot(oldnlay)   = icem%vtp(oldnlay) - Tgrad(oldnlay)*dh_old/2
      else if (oldnlay==1.and.icem%hice>epsil1) then

         rksn = rkice * (icem%rhosnw/rhow)**1.885
         rkeff = (rksn*rkice)/(icem%hice/2*rksn + icem%hsnw*rkice)
         Tgrad(oldnlay) = rkeff * (icem%tsrf -icem%vtp(oldnlay)) / &
                          (rkice*(icem%hice/2+icem%hsnw))
         tbot(oldnlay)  = icem%vtp(oldnlay) - Tgrad(oldnlay)*dh_old/2
      else
         tbot = tice_m
         Tgrad = 0.
      endif

      !Calculate temperature of new layers
      T=0.
      do i=1,newnlay

         ! Get indices of old layers entirely containing new layer
         oind    = real(i*oldnlay)/real(newnlay)
         oind_hi = min(ceiling(oind),oldnlay)
         oind    = real((i-1)*oldnlay)/real(newnlay)
         oind_lo = max(ceiling(oind),1)


         e=0.
         ! calculate energy of the new layer
         do j=oind_lo,oind_hi

            ! Upper and lower coordinate of "sublayer"
            up = min(j*dh_old,i*dh)
            dn = max((j-1)*dh_old,(i-1)*dh)

            ! zdn, zup are coordinates relative to bottom of old layer
            zup   = up - (j-1)*dh_old
            zdn   = dn - (j-1)*dh_old

            ! Vertically integrated energy of the new layer
            e = e + (tbot(j) + .5*Tgrad(j)*zup)*zup &
                  - (tbot(j) + .5*Tgrad(j)*zdn)*zdn

         end do

         ! Calculate temperature of new layer. Density and specific
         ! heat capacity falls out in the division and is omitted
         ! here and above.
         T(i) = e / dh

      end do
      changeclass%vtp(1:newnlay)         = T(1:newnlay)
      changeclass%vtp(newnlay+1:nlaymax) = icem%vtp(newnlay)
      changeclass%nlay                   = newnlay
   else if (icem%fice<=epsil1) then
      changeclass%nlay                   = newnlay
      changeclass%vtp                    = tice_m
   end if

end function changeclass





subroutine snow_albedo_update(rain,snow,rainfall,qtotice,albs)
   use mod_icestate
   implicit none
   logical, intent(in)    :: rain,snow
   real,    intent(in)    :: rainfall,qtotice(nthick)
   real,    intent(inout) :: albs(nthick)

   ! Always -- albedo decays linearly
   albs = albs - taua * dtt / 86400.

   ! New snowfall increases albedo, while rain decreases it
   if (rain) then 
      albs = (albs - albsmin) * exp(- tauf * dtt / 86400.) + albsmin
   else if (snow) then
      albs = albs + rainfall  * (albsmax -albs) / wnew
   end if

   ! Melt decreases albedo
   if (.not. rain) then
      where (qtotice>0.)
         albs = (albs - albsmin) * exp(- tauf * dtt / 86400.) + albsmin
      end where
   end if

   ! Avoid albedo values outside of limits
   albs = min(max(albs,albsmin),albsmax)

end subroutine




! ==================================================================================
! ==============================function clear_ice==================================
! ==================================================================================
#if defined (SSNOWD)

function clear_ice_scal(icex,numlay,ncv)
      use mod_icestate
      implicit none
      type(t_ice), intent(in) :: icex 
      integer,     intent(in) :: numlay
      real,        intent(in) :: ncv
      type(t_ice)             :: clear_ice_scal

      clear_ice_scal%qstore   = 0.
      clear_ice_scal%albs     = albsmin
      clear_ice_scal%fice     = 0.
      clear_ice_scal%hice     = 0.
      clear_ice_scal%hsnw     = 0.
      clear_ice_scal%rhosnw   = rhosnwmin
      clear_ice_scal%tsrf     = tice_m
      clear_ice_scal%vtp      = tice_m
      clear_ice_scal%nlay     = numlay
      clear_ice_scal%hprcp    = 0. 
      clear_ice_scal%hmelt    = 0.
      clear_ice_scal%cv       = ncv
#if defined(ICEAGE)
     clear_ice_scal%age      = 0.
#endif      
   end function clear_ice_scal
!-----------------------------------------------------------------------------------
  function clear_ice_vec(icex,numlay,ncv)
      use mod_icestate
      implicit none
      type(t_ice), intent(in), dimension(:)          :: icex
      integer,     intent(in), dimension(:)          :: numlay
      real,        intent(in), dimension(:)          :: ncv
      type(t_ice),             dimension(size(icex)) :: clear_ice_vec
      integer :: hk

      clear_ice_vec%qstore   = 0.
      clear_ice_vec%albs     = albsmin
      clear_ice_vec%fice     = 0.
      clear_ice_vec%hice     = 0.
      clear_ice_vec%hsnw     = 0.
      clear_ice_vec%rhosnw   = rhosnwmin
      clear_ice_vec%tsrf     = tice_m
      forall (hk=1:size(icex))
         clear_ice_vec(hk)%vtp   = tice_m
      end forall
      clear_ice_vec%nlay     = numlay
      clear_ice_vec%hprcp    = 0. 
      clear_ice_vec%hmelt    = 0.
      clear_ice_vec%cv       = ncv
#if defined(ICEAGE)
     clear_ice_vec%age       = 0.
#endif      
   end function clear_ice_vec

#else

   function clear_ice_scal(icex,numlay)
      use mod_icestate
      implicit none
      type(t_ice), intent(in) :: icex 
      integer,     intent(in) :: numlay
      type(t_ice)             :: clear_ice_scal

      clear_ice_scal%qstore   = 0.
      clear_ice_scal%albs     = albsmin
      clear_ice_scal%fice     = 0.
      clear_ice_scal%hice     = 0.
      clear_ice_scal%hsnw     = 0.
      clear_ice_scal%rhosnw   = rhosnwmin
      clear_ice_scal%tsrf     = tice_m
      clear_ice_scal%vtp      = tice_m
      clear_ice_scal%nlay     = numlay
#if defined(ICEAGE)
     clear_ice_scal%age      = 0.
#endif      
   end function clear_ice_scal
!-----------------------------------------------------------------------------------
   function clear_ice_vec(icex,numlay)
      use mod_icestate
      implicit none
      type(t_ice), intent(in), dimension(:)          :: icex
      integer,     intent(in), dimension(:)          :: numlay
      type(t_ice),             dimension(size(icex)) :: clear_ice_vec
      integer :: hk

      clear_ice_vec%qstore   = 0.
      clear_ice_vec%albs     = albsmin
      clear_ice_vec%fice     = 0.
      clear_ice_vec%hice     = 0.
      clear_ice_vec%hsnw     = 0.
      clear_ice_vec%rhosnw   = rhosnwmin
      clear_ice_vec%tsrf     = tice_m
      forall (hk=1:size(icex))
         clear_ice_vec(hk)%vtp   = tice_m
      end forall
      clear_ice_vec%nlay     = numlay
#if defined(ICEAGE)
     clear_ice_vec%age       = 0.
#endif      
   end function clear_ice_vec

#endif



! ==================================================================================
! ==============================ice operator (+) ===================================
! ==================================================================================
   function add_ice(ice1,ice2)
      use mod_icestate
      implicit none
      type(t_ice), intent(in) :: ice1,ice2
      type(t_ice)             :: add_ice

      integer :: i
      real    :: dh, dh1, dh2

      ! addition of ice with different # of layers is not supported.
      if (ice1%nlay.ne.ice2%nlay) then 
         print *,ice1%nlay,ice2%nlay
         stop 'ice operator (+), number of layers different. mod_icestate_tools.F90'
      end if 

      add_ice%fice   = ice1%fice + ice2%fice
      IF (add_ice%fice > epsil1) THEN
        add_ice%hice     = (ice1%fice*ice1%hice    + ice2%fice*ice2%hice )  /add_ice%fice
        add_ice%hsnw     = (ice1%fice*ice1%hsnw    + ice2%fice*ice2%hsnw )  /add_ice%fice
        add_ice%albs     = (ice1%fice*ice1%albs    + ice2%fice*ice2%albs)   /add_ice%fice
        add_ice%nlay     = ice1%nlay

        ! NB:qstore is not distributed in the vertical, it is an integrand for
        ! the ice layer
        !add_ice%qstore   = (ice1%fice*ice1%hice*ice1%qstore  + ice2%fice*ice2%hice*ice2%qstore) &
        !                   /(add_ice%hice*add_ice%fice)
        add_ice%qstore   = (ice1%fice*ice1%qstore  + ice2%fice*ice2%qstore) &
                           /(add_ice%fice)
        add_ice%tsrf     = (ice1%fice*ice1%hice*ice1%tsrf + ice2%fice*ice2%hice*ice2%tsrf) &
                           /(add_ice%hice*add_ice%fice)
#if defined(SSNOWD)
       !CV depend on ice type and is fixed for a given ice category.
       add_ice%cv        = ice1%cv
       add_ice%hprcp     = (ice1%fice*ice1%hprcp    + ice2%fice*ice2%hprcp )  /add_ice%fice
       add_ice%hmelt     = (ice1%fice*ice1%hmelt    + ice2%fice*ice2%hmelt )  /add_ice%fice
#endif   
#if defined (ICEAGE)
       add_ice%age      = (ice1%fice*ice1%hice*ice1%age + ice2%fice*ice2%hice*ice2%age) &
                           /(add_ice%hice*add_ice%fice)
#endif
      ELSE
        add_ice%hice     = 0.
        add_ice%hsnw     = 0.
        add_ice%albs     = albsmin
        add_ice%qstore   = 0.
        add_ice%tsrf     = tice_m
        add_ice%nlay     = ice1%nlay
#if defined(SSNOWD)
       !CV depend on ice type and is fixed for a given ice category 
       add_ice%cv        = ice1%cv
       add_ice%hprcp     = 0.
       add_ice%hmelt     = 0.       
#endif       
#if defined (ICEAGE)
       add_ice%age      = 0.
#endif
      ENDIF

      IF (add_ice%fice * add_ice%hsnw > epsil1) THEN
        add_ice%rhosnw = (ice1%fice*ice1%hsnw*ice1%rhosnw +&
                          ice2%fice*ice2%hsnw*ice2%rhosnw) /(add_ice%hsnw*add_ice%fice)
      ELSE
        add_ice%rhosnw = rhosnwmin
      ENDIF
     
      ! Get heat capacity factors for all layers
      dh1 = ice1%hice/ice1%nlay
      dh2 = ice2%hice/ice2%nlay
      dh  = add_ice%hice/add_ice%nlay 

      ! This way energy is conserved (rhoice and cpice falls out)
      if (add_ice%fice > epsil1) then
         add_ice%vtp(:) = ( ice1%fice*ice1%vtp(:)*dh1 + ice2%fice*ice2%vtp(:)*dh2 ) / &
                           ( add_ice%fice*dh )
      else
         add_ice%vtp(:) = tice_m
      end if

   end function add_ice




   subroutine combine_ice(icem,psim,tml)
      use mod_icestate
      implicit none

      type(t_ice), dimension(nthick), intent(inout) :: icem
      type(t_ice), dimension(nthick), intent(inout) :: psim
      real, intent(in) :: tml

      integer :: hh,hk

      DO hk = 1, nthick
        hh = nthick
        do while (hh >= 1)

           ! Add together iceclasses fallen into same thickness category.
           if ( thickl(hh)<psim(hk)%hice.and.psim(hk)%fice>epsil1 ) then
              ! This ice initially came from class hk, the vertical temperature profile
              ! should be changed to that of class hh if necessary (number of layers)
#if defined(SSNOWD)
              psim(hk)   = changeclass_SSD(psim(hk),nlay(hh),tml,cvsnw(hh))
#else
              psim(hk)   = changeclass(psim(hk),nlay(hh),tml)
#endif             
              icem(hh)   = icem(hh) + psim(hk)
              hh = 0  ! Exit inner loop
           endif

           hh = hh - 1
        enddo    
      enddo

   end subroutine combine_ice


#if defined(SSNOWD)
   

   ! -----------------------------------------------------------------------
! --------------------- FUNCTION CHANGECLASS ----------------------------
! -----------------------------------------------------------------------
! If ice initially in one class through some process (ice growth/melt, 
! advection, redistribution) moves into another class, this routine should
! be called. It ensures that the temperature profile with the old number of
! layers is represented with the new number of layers. This is done through
! conservation of energy. CV of snow distribution is adjusted. 
   function changeclass_SSD(icem,newnlay,tice_fXXX,newcv)
   use mod_icestate
   use m_get_erfc 
   implicit none
   integer,     intent(in)  :: newnlay
   type(t_ice), intent(in)  :: icem
   type(t_ice)              :: changeclass_SSD
   real,        intent(in)  :: tice_fXXX
   real,        intent(in)  :: newcv
   integer :: i,oind_lo, oind_hi, j, oldnlay, curr_old_lay,passage
   real    :: oind, dh_old, e, tmpdh, T(newnlay), dh, up, dn, zup, zdn, dheff
   real    :: tbot(nlaymax), Tgrad(nlaymax), rksn, rkeff
   real    :: fsnwcov,xi1,xi2,c1,c2,lam1,zz1,bb

   changeclass_SSD = icem

   !if(.not.(changeclass_SSD%hmelt>=0)) then
   !        write(*,*) 'hmelt NaN av cc',changeclass_SSD%hmelt, changeclass_SSD%hprcp
   !        write(*,*) 'hmelt NaN av cc',newnlay,newcv,icem%cv
    !endif
   ! passage=0

   ! Old number of layers
   oldnlay = changeclass_SSD%nlay

   ! Old vertical increment
   dh_old = icem%hice/changeclass_SSD%nlay

   ! New vertical increment
   dh     = icem%hice/newnlay

   if ((oldnlay.ne.newnlay).and.icem%fice>epsil1 ) then 

      if (oldnlay>1) then

         do j = 1,oldnlay-1
            Tgrad(j)  = (icem%vtp(j+1) - icem%vtp(j))/dh_old
            tbot(j) = icem%vtp(j) - Tgrad(j)*dh_old/2
         end do
         Tgrad(oldnlay)  = Tgrad(oldnlay-1)
         tbot(oldnlay)   = icem%vtp(oldnlay) - Tgrad(oldnlay)*dh_old/2
      else if (oldnlay==1.and.icem%hice>epsil1) then

         rksn = rkice * (icem%rhosnw/rhow)**1.885
#if defined(HEAT_CONDUC)
         changeclass_SSD%hsnw = average_depth(changeclass_SSD)
         rkeff = corr_factor(changeclass_SSD) *   &
         &           (rksn*rkice)/(changeclass_SSD%hice/2*rksn + changeclass_SSD%hsnw*rkice)
#else
         rkeff = (rksn*rkice)/(icem%hice/2*rksn + icem%hsnw*rkice)
#endif         
         Tgrad(oldnlay) = rkeff * (icem%tsrf -icem%vtp(oldnlay)) / &
                          (rkice*(icem%hice/2+icem%hsnw))
         tbot(oldnlay)  = icem%vtp(oldnlay) - Tgrad(oldnlay)*dh_old/2
      else
         tbot = tice_m
         Tgrad = 0.
      endif

      !Calculate temperature of new layers
      T=0.
      do i=1,newnlay

         ! Get indices of old layers entirely containing new layer
         oind    = real(i*oldnlay)/real(newnlay)
         oind_hi = min(ceiling(oind),oldnlay)
         oind    = real((i-1)*oldnlay)/real(newnlay)
         oind_lo = max(ceiling(oind),1)


         e=0.
         ! calculate energy of the new layer
         do j=oind_lo,oind_hi

            ! Upper and lower coordinate of "sublayer"
            up = min(j*dh_old,i*dh)
            dn = max((j-1)*dh_old,(i-1)*dh)

            ! zdn, zup are coordinates relative to bottom of old layer
            zup   = up - (j-1)*dh_old
            zdn   = dn - (j-1)*dh_old

            ! Vertically integrated energy of the new layer
            e = e + (tbot(j) + .5*Tgrad(j)*zup)*zup &
                  - (tbot(j) + .5*Tgrad(j)*zdn)*zdn

         end do

         ! Calculate temperature of new layer. Density and specific
         ! heat capacity falls out in the division and is omitted
         ! here and above.
         T(i) = e / dh

      end do
      changeclass_SSD%vtp(1:newnlay)         = T(1:newnlay)
      changeclass_SSD%vtp(newnlay+1:nlaymax) = icem%vtp(newnlay)
      changeclass_SSD%nlay                   = newnlay
      ! CV is changed since snow distribution depends on ice class. 
      changeclass_SSD%cv                     = newcv
      ! If snow cover fraction is not equal to 1, hprcp and hmelt are adjusted
      ! so that snow cover fraction and snow volume are kept constant despite CV's adjustement
      !Snow distribution is adjusted if Cv is modified
      if(.not.(changeclass_SSD%cv == icem%cv )) then
      ! If snow cover fraction is not equal to 1, hprcp and hmelt are adjusted
      ! so that snow cover fraction and snow volume are kept constant despite CV's adjustement
      if(changeclass_SSD%hmelt>epsil1 .and. changeclass_SSD%hprcp>epsil1) then
         passage=1
         !fsnwcov=snow_frac(icem)
         fsnwcov=snow_frac(changeclass_SSD)
         xi1  = (log(1+icem%cv**2))**0.5
         lam1 = log(icem%hprcp)-0.5*xi1**2
         zz1  = (log(icem%hmelt)-lam1)/xi1

         xi2=(log(1+changeclass_SSD%cv**2))**0.5
         c1 = 0.5*get_erfc((zz1-xi1)/2**0.5)
         c2 = 0.5*get_erfc((zz1-xi2)/2**0.5)

         bb=(icem%hmelt/icem%hprcp*sqrt(1+icem%cv**2))**(xi2/xi1)/sqrt(1+changeclass_SSD%cv**2)
         
         changeclass_SSD%hprcp = max(0.,(icem%hprcp*c1-icem%hmelt*fsnwcov)/(c2-fsnwcov*bb))
         changeclass_SSD%hmelt = bb*changeclass_SSD%hprcp
      end if
      end if
   else if (icem%fice<=epsil1) then
      changeclass_SSD%nlay                   = newnlay
      changeclass_SSD%vtp                    = tice_m
      changeclass_SSD%cv                     = newcv
   end if

end function changeclass_SSD

   FUNCTION snow_frac_scal(icem)
   ! Compute sonw cover fraction based on ssnowd parameters
     use m_get_erfc
     use mod_icestate
     implicit none
     type(t_ice), intent(in)  :: icem
     real lam,zz,xi
     real snow_frac_scal
     
! No snow on the ground at the previous time step     
      if(icem%hprcp<epsil1) then
         snow_frac_scal = 0.
! Snow on the ground at the previous time step
      else
! Melt has occured
         if(icem%hmelt > epsil1) then
             xi=(log(1+icem%cv**2))**(0.5)
             lam=log(icem%hprcp)-0.5*xi**2
             zz=(log(icem%hmelt)-lam)/xi
! Compute the snow cover fraction under the assumption of the erf distribution curve
             snow_frac_scal = 0.5*get_erfc(zz/2**0.5)
! The ground is covered by snow
          else 
             snow_frac_scal = 1.
          endif
      endif 
     
   END FUNCTION snow_frac_scal

   FUNCTION snow_frac_vec(icem)
   ! Compute sonw cover fraction based on ssnowd parameters
     use m_get_erfc
     use mod_icestate
     implicit none
     type(t_ice),dimension(nthick), intent(in)  :: icem
     real lam,zz,xi
     real snow_frac_vec(nthick)
     integer hk

     do hk=1,nthick
     ! No snow on the ground at the previous time step
       if(icem(hk)%hprcp<epsil1) then
         snow_frac_vec(hk) = 0.
       ! Snow on the ground at the previous time step 
       else
          ! Melt has occured
          if(icem(hk)%hmelt>epsil1)then
            xi=(log(1+icem(hk)%cv**2))**(0.5)
            lam=log(icem(hk)%hprcp)-0.5*xi**2
            zz=(log(icem(hk)%hmelt)-lam)/xi
            !  Compute the snow cover fraction under the assumption of the erf distribution curve 
            snow_frac_vec(hk) = 0.5*get_erfc(zz/2**0.5)
          else
            snow_frac_vec(hk) = 1.
          endif
       end if
     end do
         
   END FUNCTION snow_frac_vec

   FUNCTION average_depth_scal(icem)
   ! Compute average snow depth based on ssnowd parameters
     use m_get_erfc
     use mod_icestate
     implicit none
     type(t_ice),  intent(in)  :: icem
     real lam,zz,xi
     real average_depth_scal

     if(icem%hprcp<=epsil1)then
       average_depth_scal = 0.
     else     
       if(icem%hmelt<=epsil1)then
          average_depth_scal = icem%hprcp
       else   
         xi=(log(1+icem%cv**2))**(0.5)
         lam=log(icem%hprcp)-0.5*xi**2
         zz=(log(icem%hmelt)-lam)/xi
         average_depth_scal = 0.5*icem%hprcp*get_erfc((zz-xi)/2**0.5)-  &
     &  icem%hmelt*0.5*get_erfc(zz/2**0.5)
       end if
    endif   
   END FUNCTION average_depth_scal

   FUNCTION average_depth_vec(icem)
   ! Compute average snow depth based on ssnowd parameters
     use m_get_erfc
     use mod_icestate
     implicit none
     type(t_ice),dimension(nthick),  intent(in)  :: icem
     real lam,zz,xi
     real average_depth_vec(nthick)
     integer hk

     do hk=1,nthick
     if(icem(hk)%hprcp<=epsil1)then
          average_depth_vec(hk)  = 0.
     else     
       if(icem(hk)%hmelt<=epsil1)then
          average_depth_vec(hk) = icem(hk)%hprcp
       else   
         xi=(log(1+icem(hk)%cv**2))**(0.5)
         lam=log(icem(hk)%hprcp)-0.5*xi**2
         zz=(log(icem(hk)%hmelt)-lam)/xi
         average_depth_vec(hk) = 0.5*icem(hk)%hprcp*get_erfc((zz-xi)/2**0.5)-  &
     &  icem(hk)%hmelt*0.5*get_erfc(zz/2**0.5)
       end if
    endif  
     enddo
     
   END FUNCTION average_depth_vec


   function fg(xx,Da,l,ex)
   use m_get_erfc
   use mod_icestate
   real :: xx,Da,l,ex
   real :: zz,gam

   real :: fg
   !define fg function, used with fixe point method
     
     zz=(log(xx)-l)/ex
     gam=0.5*get_erfc(zz/2**0.5)
     fg=(0.5*exp(l+0.5*ex**2)*get_erfc((zz-ex)/(2**0.5))-Da)/gam

   end function fg

   function fixed_point(N,icem)
   use mod_icestate

   integer, intent(in)                           :: N
   type(t_ice),dimension(nthick), intent(inout)  :: icem

   real     :: p0(nthick)
   real     :: hmoy(nthick)
   real     :: lam(nthick)
   real     :: xi(nthick)
   real     :: fixed_point(nthick)


   !local variable
   integer i
   real, dimension(N,nthick) :: p
   real :: tol
   integer hk

   p0   = icem%hmelt
   hmoy = icem%hsnw
   xi   = (log(1+icem%cv**2))**(0.5)
   lam  = log(icem%hprcp)-0.5*xi**2

   
   tol=1e-5
   p(1,:)=p0(:)

  do hk=1,nthick
    if(hmoy(hk)==0. .and. icem(hk)%hprcp ==0.) then
       fixed_point(hk)=0.
    else 
       lam(hk)=log(icem(hk)%hprcp)-0.5*xi(hk)**2
       i=1
       do while(i<N)
        p(i+1,hk)=fg(p(i,hk),hmoy(hk),lam(hk),xi(hk))
        if(abs(p(i+1,hk)-p(i,hk))<tol) then
          fixed_point(hk)=p(i+1,hk)
          exit
        endif
       i=i+1
       enddo
    endif
  enddo

  do hk=1,nthick
    fixed_point(hk)=max(0.,fixed_point(hk))
  end do     
  end function fixed_point

#if defined (HEAT_CONDUC)

  FUNCTION corr_factor(icem)
   ! Compute conductive correction factor G depending on snow heterogeneities
   !       F = G ks.ki/(hi.ks+hs.ki)(Tf-Tsrf) 
   ! where :  hs : average snow depth (m)
   !          hi : ice thickness (m)
   !          Tf : freezing temperature of sea water (K)
   !          Tsrf : surface temperature (K)
   !          ks and ki : snow and ice thermal conductivity (W/(K.m))
   ! G depends on ice thickness, snow depth, CV
     use m_get_erfc
     use mod_icestate
     implicit none
     type(t_ice),  intent(in)  :: icem
     real lam,zz,xi,hlmin,zdm,gam,hmoy,rksnw,rkeq,epss,m_temp,I
     real int_max
     real corr_factor
     integer error

     m_temp=icem%hmelt
     xi=(log(1+icem%cv**2))**0.5
     lam=log(icem%hprcp)-0.5*xi**2
     zdm=(log(icem%hmelt)-lam)/xi;
     gam=snow_frac(icem)
     hmoy=average_depth(icem)

     int_max = 4  ! Maximum depth for integration

     ! Snow conductivity
      rksnw   = rkice*(icem%rhosnw / rhow)**1.885
      rkeq    = rksnw*rkice / (rksnw+rkice)

      hlmin= rksnw/(rksnw+rkice)*icem%hice
    
      !Effect of snow cover heterogeneities are included if snow depth is higher
      !than 1 cm. No correction is applied if hs<1 cm.
      if(hmoy>0.01) then
       if(hlmin<0.1) then
         epss=(0.1/rkeq-icem%hice/rkice)*rksnw;
         zz=(log(icem%hmelt+epss)-lam)/xi;

         ! Integral is solved numericallly using the 8 panel Newton-Cotes's rule
         call QANC8(epss+icem%hmelt,int_max,1e-4,1e-4,I,icem)

         corr_factor=max(1.,(rkice*hmoy+rksnw*icem%hice)*((1-0.5*get_erfc(zz/2**0.5))/((rksnw+rkice)*0.1)+I))
         if(.not.(corr_factor>=0)) then
           write(*,*) 'NaN corfac <0',corr_factor           
         endif
      else
 
         ! Integral is solved numericallly using the 8 panel Newton-Cotes's rule
         call QANC8(icem%hmelt+1e-15,4.,1e-4,1e-4,I,icem)

         corr_factor=max(1.,(rkice*hmoy/rksnw+icem%hice)*((1.-gam)/icem%hice+rksnw*I))
         if(.not.(corr_factor>=0)) then
           write(*,*) 'NaN corfac <0',corr_factor           
          endif
      end if
     else
      corr_factor=1.
     endif

     
   END FUNCTION corr_factor

real Function func(x,icem)
!  Function used for integration
       use mod_icestate
       real x
       type(t_ice),  intent(in)  :: icem
       real lam,xi,b,rksnw

       ! Snow conductivity
       rksnw   = rkice*(icem%rhosnw / rhow)**1.885
       
       xi=(log(1+icem%cv**2))**0.5
       lam=log(icem%hprcp)-0.5*xi**2
       
       func = 1./(sqrt(2*pi)*xi*(rkice*(x-icem%hmelt)+rksnw*icem%hice)*x)*exp(-(log(x)-lam)**2./(2*xi**2.))
       return
End   



     SUBROUTINE QANC8 (A,B,AERR,RERR,RES,icem)
!
!     INTEGRATE A REAL FUNCTION FCT(X) FROM X=A TO X=B, WITH
!     GIVEN ABSOLUTE AND RELATIVE PRECISIONS, AERR, RERR. 

!     INPUTS:
!     FCT     EXTERNAL USER-DEFINED FUNCTION FOR ANY X VALUE
!             IN INTERVAL (A,B)
!     A,B     LIMITS OF INTERVAL
!     AERR,RERR   RESPECTIVELY ABSOLUTE ERROR AND RELATIVE ERROR
!                 REQUIRED BY USER
!
!     OUTPUTS:
!     RES     VALUE OF INTEGRAL
!     ERR     ESTIMATED ERROR
!     NBF     NUMBER OF NECESSARY FCT(X) EVALUATIONS
!     FLG     INDICATOR
!             = 0.0       CORRECT RESULT
!             = NNN.RRR   NO CONVERGENCE DU TO A SINGULARITY.
!             THE SINGULAR POINT ABCISSA IS GIVEN BY FORMULA:
!             XS = B-.RRR*(B-A)

!     REFERENCE :
!     FORSYTHE,G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
!     COMPUTATIONS. PRENTICE-HALL, INC.
! -----------------------------------------------------------------------
 use mod_icestate
      
      implicit none
      type(t_ice),  intent(in)  :: icem
      real AERR,RERR,RES


      !REAL *8 (A-H,O-Z)
      REAL *8 A,B
      REAL QR(31),F(16),X(16),FS(8,30),XS(8,30)
      REAL W0,W1,W2,W3,W4,COR,TOL1,TEMP
      REAL SUM,X0,QP,F0,PAS1,PAS
      REAL ERR,QL,QN,QD,ERR1
      INTEGER LMIN,LMAX,LOUT,NMAX
      INTEGER NFIN,NIM,I,J,L 
      REAL NBF,FLG
      LMIN = 1
      LMAX = 30
      LOUT = 6
      NMAX = 5000
      NFIN = NMAX-8*(LMAX-LOUT+2**(LOUT+1))
      W0   =   3956.D0/14175.D0
      W1   =  23552.D0/14175.D0
      W2   =  -3712.D0/14175.D0
      W3   =  41984.D0/14175.D0
      W4   = -18160.D0/14175.D0
      FLG  = 0.D0
      RES  = 0.D0
      COR  = 0.D0
      ERR  = 0.D0
      SUM  = 0.D0
      NBF  = 0
      IF (A.EQ.B) RETURN
      L = 0
      NIM = 1
      X0  = A
      X(16) = B
      QP  = 0.D0
      F0   = func(X0,icem)
      PAS1  = (B-A)/16.D0
      X(8)  = (X0+X(16))*.5D0
      X(4)  = (X0+X(8))*.5D0
      X(12) = (X(8)+X(16))*.5D0
      X(2)  = (X0+X(4))*.5D0
      X(6)  = (X(4)+X(8))*.5D0
      X(10) = (X(8)+X(12))*.5D0
      X(14) = (X(12)+X(16))*.5D0
      DO 25 J = 2,16,2
      F(J) = func(X(J),icem)
   25 CONTINUE
      NBF = 9
   30 X(1)  = (X0+X(2))*.5D0
      F(1) = func(X(1),icem)
      DO 35 J = 3,15,2
      X(J)  = (X(J-1)+X(J+1))*.5D0
   35 F(J) = func(X(J),icem)
      NBF = NBF+8
      PAS = (X(16)-X0)/16.D0
      QL  = (W0*(F0+F(8))+W1*(F(1)+F(7))+W2*(F(2)+F(6))        &
                                 +W3*(F(3)+F(5))+W4*F(4))*PAS
      QR(L+1) = (W0*(F(8)+F(16))+W1*(F(9)+F(15))               &
             +W2*(F(10)+F(14))+W3*(F(11)+F(13))+W4*F(12))*PAS
      QN = QL + QR(L+1)
      QD = QN - QP
      SUM = SUM + QD
      ERR1 = DABS(QD)/1023.D0
      TOL1 = DMAX1(AERR,RERR*DABS(SUM))*(PAS/PAS1)
      IF (L.LT.LMIN) GO TO 50
      IF (L.GE.LMAX) GO TO 62
      IF (NBF.GT.NFIN) GO TO 60
      IF (ERR1.LE.TOL1) GO TO 70
   50 NIM = 2*NIM
      L = L+1
      DO 52 I = 1,8
      FS(I,L) = F(I+8)
      XS(I,L) = X(I+8)
   52 CONTINUE
      QP = QL
      DO 55 I = 1,8
      F(18-2*I) = F(9-I)
      X(18-2*I) = X(9-I)
   55 CONTINUE
      GO TO 30
   60 NFIN = 2*NFIN
      LMAX = LOUT
      FLG = FLG + (B-X0)/(B-A)
      GO TO 70
   62 FLG = FLG + 1.D0
   70 RES = RES + QN
      ERR = ERR + ERR1
      COR = COR + QD/1023.D0
   72 IF (NIM.EQ.2*(NIM/2)) GO TO 75
      NIM = NIM/2
      L = L-1
      GO TO 72
   75 NIM = NIM+1
      IF (L.LE.0) GO TO 80
      QP = QR(L)
      X0 = X(16)
      F0 = F(16)
      DO 78 I = 1,8
      F(2*I) = FS(I,L)
      X(2*I) = XS(I,L)
   78 CONTINUE
      GO TO 30
   80 RES = RES + COR
      IF (ERR.EQ.0.D0) RETURN
   82 TEMP = DABS(RES) + ERR
      IF (TEMP.NE.DABS(RES)) RETURN
      ERR = 2.D0*ERR
      GO TO 82
      END


#endif
  
#endif

#if defined (ICEAGE)
subroutine increment_age(icem,rt)
!     Increase ice age by time step length
      use mod_year_info, only : year_info
      use mod_icestate
      implicit none

      type(t_ice), dimension(nthick), intent(inout) :: icem
      type(year_info),intent(in) :: rt

      icem%age = icem%age +dtt/(86400.*real(rt%daysinyear))
      
end subroutine increment_age       
#endif

END MODULE mod_icestate_tools

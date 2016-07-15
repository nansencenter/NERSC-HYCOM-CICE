! ======================================================================
!  This module contains all definitions and parameters required to run 
!  the sea ice redistribution process (through freezing, opening,
!  ridging and rafting).
! ======================================================================
!Written by: David Salas y Melia  and  Ragnhild Hansen
!Last changed 03.01.2002 by Knut Arild Liseter
! Set up for MPI 24.09.2004

MODULE mod_icestate_redist
   implicit none

   ! Parameters for the redistribution
   real,    save              :: kri            ! random coeff of ridging
   !REAL,    PARAMETER         :: avkri = 3.    ! Mean value of kri
   !REAL,    PARAMETER         :: devkri = 2.    ! Std dev of kri
   !REAL,    PARAMETER         :: avkri = 10.    ! Mean value of kri
   !REAL,    PARAMETER         :: devkri = 2.    ! Std dev of kri
   REAL,    PARAMETER         :: avkri = 6.    ! Mean value of kri Following Babko et al 2002
   REAL,    PARAMETER         :: devkri = 1.5    ! Std dev of kri
   REAL,    PARAMETER         :: gstar = .1     ! The gstar fraction of thin ice are included in ridging.
   REAL,    PARAMETER         :: h0 = .25       ! Maximum ice thickness for rafting, m
   INTEGER, PARAMETER         :: stepfunct = 50 ! Classes for the ridging coefficient repartition function

   real,    save, dimension(stepfunct) :: repfunct   ! Ridging coefficient repartition function 
   integer, save, dimension(stepfunct) :: stat       ! # times the ridging coeff. is found in k'th class

   private :: int_b, b

contains

! -----------------------------------------------------------------------
! ----------------------- SUBROUTINE CALC_REP ---------------------------
! -----------------------------------------------------------------------
! Subroutine that computes the ridging coefficient repartition function
!
   SUBROUTINE calc_rep(avkri,devkri,stat,repfunct,stepfunct)
      implicit none
      integer,                       intent(in)  :: stepfunct
      REAL,                          INTENT(in)  :: avkri,devkri
      INTEGER, DIMENSION(stepfunct), INTENT(out) :: stat 
      REAL,    DIMENSION(stepfunct), INTENT(out) :: repfunct

      INTEGER  :: i
      real     :: pi2

      pi2=4*acos(0.)

      ! repfunct: Ridging coefficient repartition function
      stat     = 0
      repfunct = 0.
      DO i = 2,stepfunct
         repfunct(i) = repfunct(i-1)+4./real(stepfunct)*1./sqrt(pi2*devkri)*&
                  exp(-(2.+(real(i)-1.)/real(stepfunct-1)*4.-4.)**2/  &
                  (2*devkri**2))   
      END DO 
      repfunct = repfunct/repfunct(stepfunct)   
   END SUBROUTINE calc_rep  


   function icestate_calcah(icem)
      use mod_icestate
      implicit none

      type(t_ice),intent(in) :: icem(nthick)
      real :: icestate_calcah(0:nthick)

      real :: sigf,ficetot,gh 
      integer :: cl_lim,hk


      ficetot = SUM(icem%fice)

      ! note that after advection ficetot might be > fice_max  here
      ! sigf denotes this fraction (which must be redistributed in
      ! some way)
      sigf = max(0.,fice_max - ficetot)


      ! Calculate classes included in redistribution  (cl_lim)
      DO hk = 1,nthick
         sigf = sigf + icem(hk)%fice
         IF (sigf > gstar) EXIT 
      END DO
      cl_lim = hk - 1  


      ! Calculation of a(h), probability that ice of thickness h should
      ! participate to rafting or ridging. Notice that the following was
      ! adapted from Thorndike, 1975 (JGR), to a discrete case. Density 
      ! probabilities were actually integrated on ice classes.
      ! Function int_b is the integral of function b(h) as defined by 
      ! Thorndike (relation (7)). (See further down in code)
      icestate_calcah = 0.
      gh             = max(0.,fice_max - ficetot)


      ! No redistribution of open water area ....
      icestate_calcah(0)          = int_b(min(gh,gstar))


      ! icestate_calcah(hk) for every class. All of these classes are entirely involved in ridging.
      DO hk = 1, cl_lim     
         icestate_calcah(hk)    = - int_b(gh)
         gh                     = gh + icem(hk)%fice
         icestate_calcah(hk)    = icestate_calcah(hk) + int_b(gh)
      END DO


      ! icestate_calcah(hk) for last class (which could be partially involved in ridging)
      icestate_calcah(cl_lim+1) = int_b(gstar) - int_b(min(gh,gstar)) 

   end function icestate_calcah



   function icestate_calcahcum(icem)
      use mod_icestate
      implicit none

      type(t_ice),intent(in) :: icem(nthick)
      real :: icestate_calcahcum(0:nthick)

      real :: sigf,ficetot,gh 
      integer :: cl_lim,hk


      ficetot = SUM(icem%fice)

      ! note that after advection ficetot might be > fice_max  here
      ! sigf denotes this fraction (which must be redistributed in
      ! some way)
      sigf = max(0.,fice_max - ficetot)


      ! Calculate classes included in redistribution  (cl_lim)
      DO hk = 1,nthick
         sigf = sigf + icem(hk)%fice
         IF (sigf > gstar) EXIT 
      END DO
      cl_lim = hk - 1  


      ! Calculation of a(h), probability that ice of thickness h should
      ! participate to rafting or ridging. Notice that the following was
      ! adapted from Thorndike, 1975 (JGR), to a discrete case. Density 
      ! probabilities were actually integrated on ice classes.
      ! Function int_b is the integral of function b(h) as defined by 
      ! Thorndike (relation (7)). (See further down in code)
      icestate_calcahcum = 0.
      gh             = max(0.,fice_max - ficetot)


      ! No redistribution of open water area ....
      icestate_calcahcum(0)          = b(min(gh,gstar))


      ! icestate_calcahcum(hk) for every class. All of these classes are entirely involved in ridging.
      DO hk = 1, cl_lim     
         icestate_calcahcum(hk)    = b(gh) * icem(hk)%fice
         gh                     = gh + icem(hk)%fice
      END DO


      ! icestate_calcahcum(hk) for last class (which could be partially involved in ridging)
      icestate_calcahcum(cl_lim+1) = b(gstar)*min(icem(cl_lim+1)%fice,sigf-gstar)

   end function icestate_calcahcum




   subroutine icestate_mechred
      use mod_icestate
      use mod_icestate_tools
      use mod_icestate_fluxes
      implicit none

      type(t_ice), dimension(nthick) :: icem, psim
      REAL :: gh,harvest,sigf,ficetot,arealoss,toomuch, ah(0:nthick)
      real :: psi_f1m(nthick), divg, tml
      logical :: ridgm
      integer :: i,j,l,hk
      integer, dimension(nthick) :: kri_eff

#if defined (MPI)
      integer :: ierr
      real*4  :: r4harvest
      include 'mpif.h'
#endif

                                                  



      ! Calculation of a random ridging coefficient (between 2 and 6). It is
      ! assumed that ridging coefficients are distributed according to a  
      ! Gaussian law, centered in avkri and of standard deviation devkri.
      ! Statistics are performed in order to check the repartition of kri
      ! among different possible classes : stat.
      !harvest=.5 ! For testing if OMP is ok
      CALL random_number(harvest)
#if defined (MPI)
      r4harvest=harvest
#if defined (NOMPIR8)
      call mpi_bcast(r4harvest,1,MPI_REAL , 0, MPI_COMM_WORLD,ierr)
#else
      call mpi_bcast(r4harvest,1,MPI_REAL4, 0, MPI_COMM_WORLD,ierr)
#endif
      harvest=r4harvest
#endif

#ifdef MP_TEST_DUMP
      harvest=0.5
#endif
      !print *,harvest
      DO hk = 1, stepfunct-1
         IF (repfunct(hk) .LE. harvest .AND. harvest .LE. repfunct(hk+1)) EXIT 
      END DO
      stat(hk) = stat(hk) + 1
      kri     = 2. + 4. * (repfunct(hk) + repfunct(hk+1)) / 2.
      !print *,'kri is ',kri,h0




      !$OMP  PARALLEL DO PRIVATE(j,l,i,icem,divg,ficetot,sigf,hk,   & 
      !$OMP   ah,ridgm,gh,toomuch,psi_f1m,psim, arealoss,kri_eff,tml) &
      !$OMP SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do l=1,isp(j)
      do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))

         ! Set things from state etc....
         icem = icestate(i,j)%ice
         tml  = icestate(i,j)%tml
         divg    = divu(i,j) * dtt
         ficetot = SUM(icem%fice)


         ! This routine is called directly after advection. Make sure that 
         ! the thickness distribution is consistent
#if defined(SSNOWD)
         psim = clear_ice(psim,nlay,cvsnw)
#else
         psim = clear_ice(psim,nlay)
#endif         
         call combine_ice(psim,icem,tml)
         icem=psim



         ! Calculation of probability that ice of thickness h should
         ! participate to rafting or ridging.  (ah)
         if ( (divg<epsil1.and.fice_max-ficetot<gstar) .or. (ficetot > fice_max) ) then
            ridgm=.true.
            ah=icestate_calcah(icem)
         else
            ah=0.
            ridgm=.false.
         end if


        
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !---------------START of redistribution calculations---------------!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! First attempt at redistribution. Use the amount specified in the 
         ! ah - distribution
         psi_f1m = 0.
#if defined(SSNOWD)
         psim = clear_ice(psim,nlay,cvsnw)
#else
         psim = clear_ice(psim,nlay)
#endif         
         if (ridgm) then
            where (icem%hice <= h0) 
               kri_eff = 2
            elsewhere
               kri_eff = kri
            endwhere

            ! Area removed from this iceclass ...
            psi_f1m     = - kri_eff / (kri_eff - 1.) * ah(1:nthick) * abs(divg)

            ! ... results in this iceclass
            psim        = icem
            psim%hice   = kri_eff * icem%hice
!            psim%hsnw   = kri_eff * icem%hsnw
#if defined (SSNOWD)
            ! Overlying snow thickness remains the same as that over original sea
            ! ice category. To conserve mass, we suppose the rest of the snow is
            ! compacted to ice and incorporated to the ice layer
            ! We assume that ridging doesn't affectthe snow cover fraction. 
            ! hmelt/hprcp is conseved when ice is ridging
            
            psim%hsnw   = icem%hsnw
            psim%hice   = psim%hice + icem%rhosnw/rhoice*(kri_eff-1)*icem%hsnw
            psim%hprcp  = icem%hprcp
            psim%hmelt  = icem%hmelt
#else
            psim%hsnw   = icem%hsnw
            psim%hice   = psim%hice + icem%rhosnw/rhoice*(kri_eff-1)*icem%hsnw
#endif             
            psim%fice   = 1. / (kri_eff - 1.) * ah(1:nthick) * abs(divg)

            ! Check that we do not remove too much ice ...
            where (psi_f1m  < - icem%fice )  
               psi_f1m   = - icem%fice 
               psim%fice = - psi_f1m / kri_eff
            endwhere

         end if




         ! Now check wether the net area loss is sufficient to 
         ! keep ficetot below fice_max:
         arealoss      =  - sum( psi_f1m + psim%fice )  
         toomuch       = max(0.,ficetot+epsil1-fice_max)
         if (toomuch>epsil1.and.arealoss<toomuch) then

            ! Revert to traditional scheme, icestate properties are rebuilt later
            ridgm = .true.
            psi_f1m   = - icem%fice   ! All ice will first be removed in redice
            psim      = icem          ! Rebuilt ice based on icem
            psim%fice = psim%fice * fice_max / ficetot
            psim%hice = psim%hice * ficetot / fice_max
            psim%hsnw = psim%hsnw * ficetot / fice_max
#if defined (SSNOWD)
            psim%hprcp  = psim%hprcp * ficetot / fice_max
            psim%hmelt  = psim%hmelt * ficetot / fice_max
#endif            
         endif





        ! Finally, for totfice < 0.01 we adjust the thickness so that totfice = 0.01
        if ( ficetot > epsil1  .and. ficetot < fice_min .and. sum(icem%fice*icem%hice) > epsil1 ) then
           where (icem%fice > epsil1)
              psim      = icem          ! New "unridged" ice
              psi_f1m   = -icem%fice    ! All ice will be rebuilt 

              ! "Unridging" process
              psim%hice = psim%hice * psim%fice * ficetot / fice_min
              psim%hsnw = psim%hsnw * psim%fice * ficetot / fice_min
#if defined (SSNOWD)
            psim%hprcp  = psim%hprcp * psim%fice * ficetot / fice_min
            psim%hmelt  = psim%hmelt * psim%fice * ficetot / fice_min
#endif              
              psim%fice =             psim%fice * fice_min/ ficetot
           end where
           ridgm     = .true.        ! yes for ridging 

           ! Clear ice classes if necessary
           where (psim%fice<epsil1.or.psim%hice<epsil1) 
#if defined(SSNOWD)
             psim=clear_ice(psim,nlay,cvsnw)
#else
             psim=clear_ice(psim,nlay)
#endif             
           endwhere
        else if ( ficetot < fice_min .and. sum(icem%fice*icem%hice) < epsil1 ) then
#if defined(SSNOWD)
         icem = clear_ice(icem,nlay,cvsnw)
#else
         icem = clear_ice(icem,nlay)
#endif           
        end if





        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !-------------------START of actual redistribution-----------------!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (ridgm) then
           
            ! Apply psi_f1, which indicates the quantity of each ice type that 
            ! should disappear :
            icem(:)%fice = icem(:)%fice + psi_f1m(:)

            ! Combine ice in icem and psim distribution to a final 
            ! distribution in icem
            call combine_ice(icem,psim,tml)


            ! Remove empty classes
            where (icem%fice<epsil1.or.icem%hice<epsil1)
#if defined(SSNOWD)
               icem = clear_ice(icem,nlay,cvsnw)
#else
               icem = clear_ice(icem,nlay)
#endif               
            endwhere


            icestate(i,j)%ice = icem
        endif


      enddo
      enddo
      enddo
      !OMP END PARALLEL DO


   end subroutine icestate_mechred



   ! =======================================================================
   ! ==============================   b   ==================================
   ! =======================================================================
   !
   ! Function b is the function b(x) as defined by 
   ! Thorndike (relation (7)). b(x) is a weighting distribution which
   ! is used to determine how much of g(h) is involved in redistrinution
   ! a(h) = b(h) * g(h)
   !
   real function  b(x)
      implicit none
      real, intent(in) :: x
      b= (1. - x/ gstar)*2./gstar
   end function b





   ! =======================================================================
   ! ============================== int_b ==================================
   ! =======================================================================
   !
   ! Function int_b is the integral of function b(h) as defined by 
   ! Thorndike (relation (7)). b(h) is a weighting distribution which
   ! is used to determine how much of g(h) is involved in redistrinution
   ! a(h) = d/dh{int_b(G(h))} = b(h) * g(h)
   !
   real function  int_b(x)
      implicit none
      real, intent(in) :: x
      int_b= 2. * x / gstar - x**2 / gstar**2
   end function int_b



end module mod_icestate_redist

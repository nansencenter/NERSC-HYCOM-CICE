
!=======================================================================
!BOP
!
! !IROUTINE: init_evp - initialize parameters needed for evp dynamics
!
! !INTERFACE:
!
      subroutine evp_init()
!
! !DESCRIPTION:
!
! Initialize parameters and variables needed for the evp dynamics
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke
!         Fluid Dynamics Group, Los Alamos National Laboratory
!
! !USES:
!
      use mod_xc
      use mod_evp
#if defined(ICE_DYN_DIAG)
       use mod_common_ice, only : strainI, strainII
#endif
      implicit none
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer  :: i, j, k


      real ::  
     &   dte             ! subcycling timestep for EVP dynamics, s
     &,  ecc             ! (ratio of major to minor ellipse axes)^2
     &,  tdamp2          ! 2(wave damping time scale T)
#if defined PERTURB_ICE
     &,  rdm             ! random number
      integer  :: io,state(16),seed(1),clock,sze
      integer, allocatable, dimension(:):: pt
#endif
#if defined METNO
      external gsl_drandinitialize,gsl_drandgamma
#else
      external DRANDINITIALIZE, DRANDGAMMA
#endif

      !ndte=120 Now read from infile.evp

      ! dynamics time step
      ! KAL - clean up the dt stuff....
      ! KAL -Note that ndyn_dt is set to one
      dyn_dt = dt/ndyn_dt        ! s



      ! elastic time step
      dte = dyn_dt/real(ndte)    ! s
      dtei = 1.0/dte              ! 1/s
      tdamp2 = 2.0*eyc*dyn_dt     ! s
      if (mnproc==1) 
     & write(*,*) 'dyn_dt = ',dyn_dt,'  dte = ',dte,
     &                  '  tdamp = ', 0.5*tdamp2
      ! major/minor axis length ratio, squared
      ecc  = 4.0
#if defined PERTURB_ICE
      call RANDOM_SEED(size=sze)
         allocate(pt(sze))
      call RANDOM_SEED(GET=pt)
      seed(1)=pt(1)
#if defined METNO
      CALL gsl_drandinitialize(seed(1),state)
      CALL gsl_drandgamma(1,5.,1.,state,rdm)
#else
      CALL DRANDINITIALIZE(1,1,seed(1),1,state,16,io)
      !CALL DRANDGAUSSIAN(1,0.,1.,state,rdm,io)
      CALL DRANDGAMMA(1,5.,1.,state,rdm,io)
#endif
      ecc=rdm
      !ecc = 4+2*exp(rdm-1.4764)
      if (mnproc==1) 
     & print *, 'the perturb value is',ecc,rdm
#endif
      ecci = 1/ecc
      ! constants for stress equation
      dte2T = dte/tdamp2                    ! unitless
      denom = 1.0/(1.0+dte2T)
c      denom2 = 1.0/(1.0+dte2T*ecc) ! follow Bouillon et al. OM 2013
      rcon_evp = 1230.*eyc*dyn_dt*dtei**2  ! kg/s
      rcon_miz = 4.*1230*eyc*dyn_dt*dtei**2/(3.141592*sqrt(2.0))

      imargin=nbdy
      do j=1-imargin,jj+imargin
       do i=1-imargin,ii+imargin
        ! Coriolis parameter
c        fcor(i,j) = 1.46e-4 ! Hibler 1979, Northern Hemisphere; 1/s
        fcor(i,j) = 2.0*omega*sin(ULAT(i,j))  ! 1/s
      if (e_itst==i .and. e_jtst==j) then
         print *,'evp_init',ULAT(i,j)
         print *,'evp_init',fcor(i,j)
         print *,'evp_init',dxt(i,j)
      end if

        ! velocity
        uvel(i,j) = 0.0       ! m/s
        vvel(i,j) = 0.0       ! m/s

       enddo
      enddo

      ! stress tensor,  kg/s^2
      do j=1-imargin,jj+imargin
       do i=1-imargin,ii+imargin
         stressp_1 (i,j) = 0.0
         stressp_2 (i,j) = 0.0
         stressp_3 (i,j) = 0.0
         stressp_4 (i,j) = 0.0
         stressm_1 (i,j) = 0.0
         stressm_2 (i,j) = 0.0
         stressm_3 (i,j) = 0.0
         stressm_4 (i,j) = 0.0
         stress12_1(i,j) = 0.0
         stress12_2(i,j) = 0.0
         stress12_3(i,j) = 0.0
         stress12_4(i,j) = 0.0

         ! initialize other variables prior to writing i.c. history file
         strtltx(i,j) = 0.0
         strtlty(i,j) = 0.0
         strintx(i,j) = 0.0
         strinty(i,j) = 0.0
         fm     (i,j) = 0.0
         prs_sig(i,j) = 0.0
         divu   (i,j) = 0.0
         shear  (i,j) = 0.0
#if defined(ICE_DYN_DIAG)
         strainI(i,j)  = 0.0
         strainII(i,j) = 0.0
#endif
         ! Init mask 
         iceumask(i,j)=.false.

       enddo
      enddo

      ! KAL - TODO - xctilr, but margin already filled here
      !KAL call bound(uvel)
      !KAL call bound(vvel)


      ! KAL - read infile.evp here

      end subroutine evp_init

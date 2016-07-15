C!****f* HYCOM_2.1.03/icemodels_step
C!
C! NAME
C!  icemodels_step - Step ice model one time step forward
C!
C! SYNOPSIS
C!  subroutine icemodels_step(m,n,restart,rt,rforce,dtime)
C!
C!
C! DESCRIPTION
C!  This is the main ice-related routine in NERSC-HYCOM. It is mostly 
C!  a wrapper for other routines which do ice thermodynamics and dynamics.
C!
C!  Integrate ice model one step forward in time. Thermodynamics are
C!  computed every time step, while dynamics can be computed with a
C!  different time step. Sea ice advection is also handled in this
C!  routine.
C!
C! INPUTS   
C!  m,n             Old and new time step indices in HYCOM
C!
C!  restart         True if starting from restart file
C!
C!  rt              Variable holding time info
C!
C!  rforce          Flag for forcing option
C!
C!  rt              Variable holding time info
C!
C!  dtime           Time step   - MORE HERE.
C!
C! SIDE EFFECTS
C!  Changes are made to ice variables, but also to surface fluxes (stress, salinity 
C!  flux and heat flux) going to hycom. Note that most of these changes are made
C!  in the subroutines called below. The only thing being changed locally
C!  are the advection when using EVP_MPI. This will change (See TODO).
C!   
C!  The variables AVEu1 and AVEv1 are average variables, they are
C!  updated between time steps, and set to zero after EVP/EVP_MPI has been
C!  called.
C!
C! WARNINGS
C!  Different CPP options will change this routine. CPP options currently
C!  implemented are ICE, EVP, EVP_MPI and ICESTATE.
C!
C!
C! PARAMETERS
C!
C!  fice_max        Maximum allowed ice concentration. From module
C!                  mod_common_ice, set in subroutine icedat.
C!
C!  fice_min        Set locally.
C!
C!
C! AUTHOR
C!  Knut Arild Liseter
C!
C! CREATION DATE
C!  Sometime in 2003
C!
C! HISTORY 
C!  Jan 16th 2006, Knut Lisaeter: Added robodoc documentation header
C! 
C!
C! SOURCE
C!

      module m_icemodels_step
      use mod_xc 
      real,save,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),private::
     &   AVEu1,AVEv1,AVEsurfU,AVEsurfV,AVEssh
      integer,save,private :: AVE1cnt
      contains


      subroutine icemodels_step(m,n,restart,rt,rforce,dtime)
      use mod_common_ice
      use mod_year_info, only: year_info
c --- --------------------------------------------------------
c --- ------------------- EVP modules-------------------------
c --- --------------------------------------------------------
#if defined (EVP)
      use m_depthave
      use mod_advem
      use m_evp_next_step
      use mod_evp, only : dragw
#endif
c
#if defined (ICESTATE) && defined (ICE)
#error CPP FLAGS You have defined both ICE and ICESTATE
#fi
c
#elif defined (ICESTATE)
c --- ------------------------------------------------------
c --- ---------------ICESTATE modules ----------------------
c --- ------------------------------------------------------
      use mod_common_ice , only :icestate2ice
      use mod_icestate_diag     ! ICESTATE diagnostic routines
      use mod_icestate_redist   ! ICESTATE redistribution routines
      use mod_icestate_transfer ! module contains funcs to send to/from ICESTATE
      use m_icestate_thermf     ! thermodynamic module
#else
c --- ------------------------------------------------------
c --- ----------------- ICE modules-------------------------
c --- ------------------------------------------------------
      use m_thermf_nersc
#endif
      use m_icemodels_advect
      implicit none 
c --- Input parameters
      type(year_info),       intent(in)    :: rt
      integer,               intent(in)    :: m,n
      logical,               intent(in)    :: restart
      character(len=5),      intent(in)    :: rforce
      real*8,                intent(in)    :: dtime
c --- Local parameters
      logical :: dynflg
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: 
     &   ouvel,ovvel
      integer :: i,j,l
      real :: offset,flxdiv, wdiff, tmpssh
      include 'common_blocks.h'
c
c --- This is true whenever EVP should run - set later
      dynflg=.false.
                             
                             
c --- ------------------------------------------------------
c --- ----------------- Ice dynamics part ------------------
c --- ------------------------------------------------------

#if defined (EVP)
c --- Compute average velocities in upper model layer(s) - (20 meters)
      if (restart) then
         AVESSH=0.
         AVEu1=0.
         AVEv1=0.
         AVEsurfU=0.
         AVEsurfV=0.
         AVE1cnt=0
      end if
c
c --- horizontal velocities integrated over top 20 meters
      ouvel=depthave(u(:,:,1:kdm,m),dpu(:,:,1:kdm,m),20.,'u') 
      ovvel=depthave(v(:,:,1:kdm,m),dpv(:,:,1:kdm,m),20.,'v') 
C$OMP PARALLEL DO PRIVATE(j,l,i,tmpssh) 
C$OMP&SCHEDULE(STATIC,jblk)
      do j=1-nbdy,jj+nbdy
C
      do l=1,isu(j)
      do i=max(1-nbdy,ifu(j,l)),min(ii+nbdy,ilu(j,l))
         AVEsurfU(i,j) = AVEsurfU(i,j) + u(i,j,1,m) + ubavg(i,j,m)
         AVEu1(i,j)    = AVEu1(i,j) + ouvel(i,j)
         AVEu1(i,j)    = AVEu1(i,j) + ubavg(i,j,m)
      end do
      end do
C
      do l=1,isv(j)
      do i=max(1-nbdy,ifv(j,l)),min(ii+nbdy,ilv(j,l))
         AVEsurfV(i,j) = AVEsurfV(i,j) + v(i,j,1,m) + vbavg(i,j,m)
         AVEv1(i,j)    = AVEv1(i,j) + ovvel(i,j)
         AVEv1(i,j)    = AVEv1(i,j) + vbavg(i,j,m)
      end do
      end do
C
      do l=1,isp(j)
      do i=max(1-nbdy,ifp(j,l)),min(ii+nbdy,ilp(j,l))
         tmpssh = (montg1(i,j)+thref*pbavg(i,j,m))/onem
         AVEssh(i,j) = AVEssh(i,j) + tmpssh
      end do
      end do
C
      end do
C$OMP END PARALLEL DO
      AVE1cnt = AVE1cnt + 1
c
c --- Hunke & Dukowicz EVP model - New MPI implementation
      dynflg = EVP_next_step(rt,baclin)
      if (dynflg) then
!         if (mnproc==1) write(lp,*) 'Calling EVP'
         AVEssh  =AVEssh  /AVE1cnt
         AVEu1   =AVEu1   /AVE1cnt
         AVEv1   =AVEv1   /AVE1cnt
         AVEsurfu=AVEsurfu/AVE1cnt
         AVEsurfv=AVEsurfv/AVE1cnt
         call xctilr(AVEu1( 1-nbdy,1-nbdy    ),1,   1, 6,6, halo_uv)
         call xctilr(AVEv1( 1-nbdy,1-nbdy    ),1,   1, 6,6, halo_vv)
         call xctilr(AVEsurfu( 1-nbdy,1-nbdy    ),1,   1, 6,6, halo_uv)
         call xctilr(AVEsurfv( 1-nbdy,1-nbdy    ),1,   1, 6,6, halo_vv)
         !call hycomtoevp(AVEsurfu,AVEsurfv,m,n)

c ---    NB - hycom2evp also transfers ICE and ICESTATE variables to EVP
         call hycomtoevp(AVEu1,AVEv1,m,n)
c
c ---    Called with Hibler type Strength parameterization (1)
         call evp(1) 
c
c ---    NB - evp2hycom transfers EVP to HYCOM
         call evptohycom
c
#if defined (ICESTATE)
c ---    NB - evp2istate transfers EVP variables to icestate
         call evp2istate
#endif
         AVEssh  =0.
         AVEu1   =0.
         AVEv1   =0.
         AVEsurfU=0.
         AVEsurfV=0.
         AVE1cnt =0
      end if
      margin=nbdy
C     Boundary were not filled in we retrieve it
      call xctilr(u(1-nbdy,1-nbdy,1,m),1,kk, nbdy,nbdy, halo_uv)
      call xctilr(v(1-nbdy,1-nbdy,1,m),1,kk, nbdy,nbdy, halo_uv)
C
c --- Update ice-ocean stresses
c --- TODO: tauxice/tauyice should perhaps be on p-grid
c --- TODO: use top 20 meters velocity here (AVEu1, AVEv1)
C$OMP PARALLEL DO PRIVATE(j,l,i,wdiff) 
C$OMP&SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         wdiff=sqrt( (iceU(i,j)-(u(i,j,1,m) + ubavg(i,j,m)))**2 +
     &               (iceV(i,j)-(v(i,j,1,m) + vbavg(i,j,m)))**2 )
         tauxice(i,j) = dragw*wdiff*
     &      (iceU(i,j) -(u(i,j,1,m) + ubavg(i,j,m)))
         tauyice(i,j) = dragw*wdiff*
     &      (iceV(i,j) -(v(i,j,1,m) + vbavg(i,j,m)))
      end do
      end do
C$OMP END PARALLEL DO
cdiag if ( itest /= -99 .and. jtest  /= -99) then
cdiag    print '(a,2f10.2)','icemodels_step - finel ice vel:',
cdiag&      iceu(itest,jtest), icev(itest,jtest)
cdiag end if

c --- ------------------------------------------------------
c --- --------------- Ice Advection part -------------------
c --- ------------------------------------------------------
c --- This became large and messy - moved it to a separate
c --- subroutine. Processes both ICE and ICESTATE advection
      call icemodels_advect()
#endif /* EVP*/


c --- ------------------------------------------------------
c --- --------------- Thermodynamic part -------------------
c --- ------------------------------------------------------
#if defined (ICESTATE)
c --- Transfer variables from hycom to icestate
      call hycom2istate(n,m)
c 
c --- Mechanical redistribution
      call icestate_mechred()  
c
c --- Thermodynamics routine
      if (mnproc==1) write(lp,'(a)') 'Running ICESTATE!!!'
      call icestate_thermf(rt,thermo,restart,dtime,rforce)
c
c --- Transfer fluxes from icestate to hycom
      call istate2hycom
c
c --- Transfer ICESTATE to common_ice for diagnostic
      call icestate2ice
#else
c --- Standard one-category thermodynamics routine
      call thermf_nersc(m,n,rt)
#endif

      end subroutine icemodels_step
      end module m_icemodels_step
!!***

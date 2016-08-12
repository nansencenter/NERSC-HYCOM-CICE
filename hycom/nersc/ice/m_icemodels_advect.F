!!****f* HYCOM_2.1.03/icemodels_advect
!!
!! NAME
!!  icemodels_advect - Step ice model one time step forward
!!
!! SYNOPSIS
!!  subroutine icemodels_advect()
!!
!!
!! DESCRIPTION
!!  This is the ice advection routine in NERSC-HYCOM. It is the top-
!!  level logic for advectong scalar fields (hicem, ficem etc).
!!
!!
!! INPUTS   
!!
!! SIDE EFFECTS
!!  called.
!!
!! WARNINGS
!!  Different CPP options will change this routine. CPP options currently
!!  implemented are ICE, EVP_MPI and ICESTATE.
!!
!!
!! PARAMETERS
!!  fice_max        Maximum allowed ice concentration. From module
!!                  mod_common_ice, set in subroutine icedat.
!!
!!  fice_min        Set locally.
!!
!!
!! AUTHOR
!!  Knut Arild Liseter
!!
!! CREATION DATE
!!  25th Feb 2008
!!
!! MODIFICATIONS
!! 9th Jun  2008 - Bug fix in diffusion
!! 9th Jun  2008 - multiplied advection variables with tenm to get
!!                 advected variables on order of 1 
!! SOURCE
!!
      module m_icemodels_advect
      use mod_xc 
      contains


      subroutine icemodels_advect()
#if defined (ICE)
      use mod_common_ice
#elif defined (ICESTATE)
      use mod_icestate   , only: icestate, nthick, nlaymax, epsil1,
     &                           epsil0,itst,jtst, 
     &                           rhosnwmin, tice_m, albi_mlt, nlay
#if defined (SSNOWD)
     &                            ,cvsnw
#endif     
#if defined(ICEAGE)
     &                            ,age_max
#endif      
      use mod_icestate_tools , only: clear_ice
#endif
      use mod_common_ice, only: iceu,icev,iceadv 
      use mod_advem
      implicit none 
      integer :: i,j,l,k,oldmargin,hk,hl,its,jts
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: 
     &   dummy1, dummy2, snwvol, icevol, snwvol0,
     &   iuflx, ivflx,newicevol, farea,
     &   ifld, ficem_old, s_iuflx, s_ivflx
#if defined(ALBSNW_EVOL)
     &   ,albsnw
#endif
#if defined(SSNOWD_ICE)
     &   ,hmlt,hpcp
#endif       
#if defined (TEST_ICE_AGE)
     &  ,crdg, c_l, ridgearea, fyage, fyageold, fyarea
      real, parameter :: kridge=4.
#endif
#if defined (ICESTATE)
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: 
     &   rhosnw,  qbrine, albs, tsrf, ficem, farea0, icevol0, 
     &   c_iuflx, c_ivflx, v_ivflx, v_iuflx
#if defined (SSNOWD)
     &   ,prcp,melt
#endif   
#if defined (ICEAGE)
     &   ,age
#endif      
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,nlaymax) :: 
     &   vtp
#endif
      real :: offset,flxdiv,factor, icedf
      real, parameter :: fice_min=0.05 ! Should be gotten from mod_common_ice
      include 'common_blocks.h'
C --- define step function;  step=1  if  xx1.ge.xx2 ; =0  if xx1.lt.xx2
      real :: xx1, xx2, step
      step(xx1,xx2)=(1.+sign(1.,xx1-xx2))*.5


c --- Keep old margin - reset at end of routine
      oldmargin=margin

c --- -------------------------------------------------------------------
c --- Setup of advection fluxes 
c --- -------------------------------------------------------------------
c
c --- Ice flux variables are now initialized on the full margin - we tile later
      margin=nbdy
C$OMP PARALLEL DO PRIVATE(j,l,i) 
C$OMP&SCHEDULE(STATIC,jblk)
      do j=1-nbdy,jdm+nbdy
         do i=1-nbdy,idm+nbdy
            dummy1(i,j)=tenm
            dummy2(i,j)=dummy1(i,j)
            iuflx(i,j)=0.
            ivflx(i,j)=0.
         end do
      end do
C$OMP END PARALLEL DO

  
c --- Set up ice variables as "layered" - so that integrating variable 
c --- over a layer retrieves origial fields. iuflx and ivflux Set up
c --- here
      margin=nbdy-1 ! due to i+1,i-1
C$OMP PARALLEL DO PRIVATE(j,l,i) 
C$OMP&SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
C
c ---    Set up u-fluxes
         do l=1,isu(j)
         do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            iuflx(i,j)=iceu(i,j)*scuy(i,j)*dummy1(i-1,j)
         end do
         end do
c ---    Set up v-fluxes
         do l=1,isv(j)
         do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            ivflx(i,j)=icev(i,j)*scvx(i,j)*dummy1(i,j-1)
         end do
         end do
      end do

c --- Set up before and after fake layer thickness 
      margin=nbdy-1 ! due to i+1,i-1
C$OMP PARALLEL DO PRIVATE(j,l,i,flxdiv,offset) 
C$OMP&SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do l=1,isp(j)
      do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
         flxdiv=((iuflx(i+1,j)  -iuflx(i,j)  )
     &          +(ivflx(i,j+1)  -ivflx(i,j)  ))*baclin*scp2i(i,j)
         dummy2(i,j)=dummy1(i,j)-flxdiv
         offset=min(0.,dummy1(i,j),dummy2(i,j))
         dummy2(i,j)=dummy2(i,j)-offset
         dummy1(i,j)=dummy1(i,j)-offset
      end do
      end do
      end do
C$OMP END PARALLEL DO


      call xctilr(dummy1(1-nbdy,1-nbdy),1,1, 6,6, halo_ps)
      call xctilr(dummy2(1-nbdy,1-nbdy),1,1, 6,6, halo_ps)
      call xctilr(iuflx (1-nbdy,1-nbdy),1,1, 6,6, halo_uv)
      call xctilr(ivflx (1-nbdy,1-nbdy),1,1, 6,6, halo_vv)


! --- -------------------------------------------------------------------
! --- Standard ice model advection + diffusion follows.
! --- -------------------------------------------------------------------
#if defined (ICE)
! APPROACH 1 : ice advection using MPDATA
! APPROACH 2 : ice advection using WENO
c
! WENO is the default advection sceheme      

#if ! defined (APPROACH_1)
#define APPROACH_2
#endif

#if defined(APPROACH_1)      
c
c --- Ice variables are now initialized on the full margin - we tile later
      margin=nbdy
C$OMP PARALLEL DO PRIVATE(j,l,i) 
C$OMP&SCHEDULE(STATIC,jblk)
      do j=1-nbdy,jdm+nbdy
         do i=1-nbdy,idm+nbdy
            icevol(i,j)=0.
            snwvol(i,j)=0.
            farea (i,j)=0.
#if defined (TEST_ICE_AGE)
            fyarea(i,j)=0.
            fyage (i,j)=0.
            fyageold(i,j)=0.
#endif
         end do
      end do
C$OMP END PARALLEL DO


c --- Set up ice variables as "layered" - so that integrating variable 
c --- over a layer retrieves origial fields. Also correct for flux
c --- divergence. Note that this may be better stated in terms of ice
c --- volume rather than a fake layer thickness.
      margin=nbdy-1 ! due to i+1,i-1
C$OMP PARALLEL DO PRIVATE(j,l,i) 
C$OMP&SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
C
c ---    Set up advection variables
         do l=1,isp(j)
         do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
C
c ---       Now these are "layered", eg farea*dummy1 gives ice conc.
            icevol(i,j)=ficem(i,j)*hicem(i,j)*tenm/dummy1(i,j)
            snwvol(i,j)=ficem(i,j)*hsnwm(i,j)*tenm/dummy1(i,j)
            farea (i,j)=ficem(i,j)*tenm/dummy1(i,j)
#if defined (TEST_ICE_AGE)
            fyarea   (i,j)=fy_frac (i,j)*tenm/dummy1(i,j)
            fyage    (i,j)=fy_age  (i,j)*fy_frac(i,j)*tenm/dummy1(i,j)
            fyageold (i,j)=fy_age  (i,j)
            ridgearea(i,j)=rdg_frac(i,j)*tenm/dummy1(i,j)
#endif /*TEST_ICE_AGE*/
         end do
         end do
      end do
C$OMP END PARALLEL DO
C

c --- Tile  and advect variables. For now use the mpdata routine - 
c --- which does not use middle time step. After advection "integrate"
c --- field to recover 
      margin=nbdy 
      call xctilr(icevol(1-nbdy,1-nbdy),1,1, 6,6, halo_ps)
      call xctilr(snwvol(1-nbdy,1-nbdy),1,1, 6,6, halo_ps)
      call xctilr(farea (1-nbdy,1-nbdy),1,1, 6,6, halo_ps)
      call advem(1,farea,farea,iuflx,ivflx,
     &           dummy1,dummy2,0.,scp2,scp2i,baclin )
      call advem(1,icevol,icevol,iuflx,ivflx,
     &           dummy1,dummy2,0.,scp2,scp2i,baclin )
      call advem(1,snwvol,snwvol,iuflx,ivflx,
     &           dummy1,dummy2,0.,scp2,scp2i,baclin )
      ficem=farea*dummy2/tenm
      icevol=icevol*dummy2/tenm
      snwvol=snwvol*dummy2/tenm

#if defined (TEST_ICE_AGE)
      call xctilr(ridgearea(1-nbdy,1-nbdy),1,1, 6,6, halo_ps)
      call xctilr(fyarea   (1-nbdy,1-nbdy),1,1, 6,6, halo_ps)
      call xctilr(fyage    (1-nbdy,1-nbdy),1,1, 6,6, halo_ps)
      call advem(1,fyarea,fyarea,iuflx,ivflx,
     &           dummy1,dummy2,0.,scp2,scp2i,baclin)
      call advem(1,fyage,fyage,iuflx,ivflx,
     &           dummy1,dummy2,0.,scp2,scp2i,baclin)
      call advem(1,ridgearea,ridgearea,iuflx,ivflx,
     &           dummy1,dummy2,0.,scp2,scp2i,baclin)
      ridgearea=ridgearea*dummy2/tenm
      fyarea   =fyarea   *dummy2/tenm
      fyage    =fyage    *dummy2/tenm
#endif /*TEST_ICE_AGE*/


      ! Artificial ridging
Cdiag tmp_new_icearea=ficem*scp2
      where (ficem > fice_max )
#if defined (TEST_ICE_AGE)

         ! This is the fraction to be ridged
         crdg=ficem-fice_max

         ! This is the level ice to be removed. 
         ! k*h_level is the thickness of ice ridged in this way.
         ! h_level is unknown. We only operate with mean ice thickness
         ! here.
         c_l=crdg*kridge/(kridge-1)

         ! Remove fy area. c_l is removed. max is needed in case
         ! my level ice (or ridged ice) is removed as well. 
         ! Note that both my and fy level
         ! ice is assumed to have the same thickness in this
         ! approximation.
         fyarea    =max(0.,fyarea -   c_l)

         ! increase in ridged area. min is needed in case ridged ice is
         ! re-ridged. (In this case ridged ice is not conserved)
         ridgearea =min(1.,ridgearea + c_l/kridge)
#endif
         ficem =fice_max
      else where (ficem<fice_min .and. ficem>0.)
#if defined (TEST_ICE_AGE)
         fyarea   =fyarea  *   (1. + fice_min-ficem)
         ridgearea=ridgearea * (1. + fice_min-ficem)
#endif /*TEST_ICE_AGE*/
         ficem=fice_min
      end where

      ! NB - Should check the logic here
      where (ficem >= fice_min )
         hicem=icevol/ficem
         hsnwm=snwvol/ficem
      end where

#if defined (TEST_ICE_AGE)
      fy_frac=fyarea
      rdg_frac=ridgearea
      where (fy_frac>1e-3)
         fy_age = fyage/fy_frac
      elsewhere
         fy_age = fyageold
      end where
#endif /*TEST_ICE_AGE*/


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Diffusion of variables starts here ..
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !  Use 1/5 "temdf2"  as diffusion velocity
      icedf=temdf2*0.2
      icevol=hicem*ficem
      snwvol=hsnwm*ficem
      ficem_old=ficem
      margin=nbdy-1 ! due to i+1, i-1
      do k=1,3
         if (k==1) then
            ifld=ficem
         elseif (k==2) then
            ifld=icevol
         elseif (k==3) then
            ifld=snwvol
         end if
!$OMP PARALLEL DO PRIVATE(j,l,i,factor)
!$OMP&         SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
           do l=1,isu(j)
             do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
               factor=icedf*aspux(i,j)*scuy(i,j)*
     &                min(step(ficem_old(i  ,j),fice_min),
     &                    step(ficem_old(i-1,j),fice_min))
               iuflx(i,j)=factor*(ifld(i-1,j)-ifld(i,  j))
             enddo
           enddo
           do l=1,isv(j)
             do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
               factor=icedf*aspvy(i,j)*scvx(i,j)*
     &          min(step(ficem_old(i,j  ),fice_min),
     &              step(ficem_old(i,j-1),fice_min))
               ivflx(i,j)=factor*(ifld(i,j-1)-ifld(i,j))
             enddo
           enddo
         enddo
C$OMP END PARALLEL DO


!$OMP PARALLEL DO PRIVATE(j,l,i,factor)
!$OMP&         SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
           do l=1,isp(j)
             do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
               factor=-delt1/scp2(i,j)
               util1(i,j)=(iuflx(i+1,j)-iuflx(i,j)
     &                 +ivflx(i,j+1)-ivflx(i,j))*factor
               ifld(i,j)=ifld(i,j)+util1(i,j)
             enddo
           enddo
         enddo
C$OMP END PARALLEL DO

         if (k==1) then
            ficem=ifld
         elseif (k==2) then
            icevol=ifld
         elseif (k==3) then
            snwvol=ifld
         end if
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Final adjustment of fields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Adjustment for small numbers after diffusion
      !where (ficem<fice_min .and. ficem>fice_tol)
      !   ficem=fice_min
      !elsewhere (ficem<=fice_tol)
      !   ficem=0.
      !   hsnwm=0.
      !   hicem=0.
      !end where
      where (ficem<fice_min .and. ficem>0.)
         ficem=fice_min
      end where
      where (ficem >= fice_min )
         hicem=icevol/ficem
         hsnwm=snwvol/ficem
      end where





      ! Tile final variables
      call xctilr(ficem( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(hicem( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(hsnwm( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
#if defined (TEST_ICE_AGE)
      call xctilr(fy_age  ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(fy_frac ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(rdg_frac( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
#endif /*TEST_ICE_AGE*/


! WENO : no diffusion is added to the transport fields      
#elif defined (APPROACH_2) 

      do j=1-nbdy,jdm+nbdy
         do i=1-nbdy,idm+nbdy
            icevol(i,j)=0.
            snwvol(i,j)=0.
            farea (i,j)=0.
#if defined (TEST_ICE_AGE)
            fyarea(i,j)=0.
            fyage (i,j)=0.
            fyageold(i,j)=0.
#endif
#if defined(ALBSNW_EVOL)
            albsnw(i,j)=0.
#endif  
#if defined(SSNOWD_ICE)
            hmlt(i,j)  =0.
            hpcp(i,j)  =0.
#endif
         end do
      end do

      margin=nbdy-1 ! due to i+1,i-1
C$OMP PARALLEL DO PRIVATE(j,l,i) 
C$OMP&SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
C
c ---    Set up advection variables
         do l=1,isp(j)
         do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
C
            icevol(i,j)=ficem(i,j)*hicem(i,j)
            snwvol(i,j)=ficem(i,j)*hsnwm(i,j)
            farea (i,j)=ficem(i,j)
#if defined (TEST_ICE_AGE)
            fyarea   (i,j)=fy_frac (i,j)
            fyage    (i,j)=fy_age  (i,j) * fy_frac(i,j)
            fyageold (i,j)=fy_age  (i,j)
            ridgearea(i,j)=rdg_frac(i,j)
#endif   
#if defined(ALBSNW_EVOL)
            albsnw(i,j)=ficem(i,j)*albsnwm(i,j)
#endif  
#if defined(SSNOWD_ICE)
           !if(.not.(hmelt(i,j).ge.0)) then
           ! write(*,*) 'hmelt NaN av ad',i0+i,j0+j
           ! write(*,*) 'hmelt NaN av ad',hmelt(i,j)
           !endif
            hmlt(i,j)  =ficem(i,j)*hmelt(i,j)
            hpcp(i,j)  =ficem(i,j)*hprcp(i,j)
#endif             
         end do
         end do
      end do
C$OMP END PARALLEL DO

c --- Tile  and advect variables. Use the WENO routine - 
c --- which does not use middle time step. 
      margin=nbdy 
      call xctilr(icevol(1-nbdy,1-nbdy),1,1, 6,6, halo_ps)
      call xctilr(snwvol(1-nbdy,1-nbdy),1,1, 6,6, halo_ps)
      call xctilr(farea (1-nbdy,1-nbdy),1,1, 6,6, halo_ps)
     
      call iceadv(farea ,iceu,icev,scuy,scvx,scp2i,scp2,baclin)
      call iceadv(icevol,iceu,icev,scuy,scvx,scp2i,scp2,baclin)
      call iceadv(snwvol,iceu,icev,scuy,scvx,scp2i,scp2,baclin)
#if defined (TEST_ICE_AGE)
      call xctilr(ridgearea(1-nbdy,1-nbdy),1,1, 6,6, halo_ps)
      call xctilr(fyarea   (1-nbdy,1-nbdy),1,1, 6,6, halo_ps)
      call xctilr(fyage    (1-nbdy,1-nbdy),1,1, 6,6, halo_ps)

      call iceadv(ridgearea,iceu,icev,scuy,scvx,scp2i,scp2,baclin)
      call iceadv(fyarea   ,iceu,icev,scuy,scvx,scp2i,scp2,baclin)
      call iceadv(fyage    ,iceu,icev,scuy,scvx,scp2i,scp2,baclin)
#endif
#if defined(ALBSNW_EVOL)
      call xctilr(albsnw(1-nbdy,1-nbdy),1,1, 6,6, halo_ps)
      call iceadv(albsnw,iceu,icev,scuy,scvx,scp2i,scp2,baclin)
#endif   
#if defined(SSNOWD_ICE)
      call xctilr(hmlt  (1-nbdy,1-nbdy),1,1, 6,6, halo_ps)
      call iceadv(hmlt  ,iceu,icev,scuy,scvx,scp2i,scp2,baclin)
      call xctilr(hpcp  (1-nbdy,1-nbdy),1,1, 6,6, halo_ps)
      call iceadv(hpcp  ,iceu,icev,scuy,scvx,scp2i,scp2,baclin)
#endif

      ficem=farea

      where (ficem > fice_max )
#if defined (TEST_ICE_AGE)

         ! This is the fraction to be ridged
         crdg=ficem-fice_max

         ! This is the level ice to be removed. 
         ! k*h_level is the thickness of ice ridged in this way.
         ! h_level is unknown. We only operate with mean ice thickness
         ! here.
         c_l=crdg*kridge/(kridge-1)

         ! Remove fy area. c_l is removed. max is needed in case
         ! my level ice (or ridged ice) is removed as well. 
         ! Note that both my and fy level
         ! ice is assumed to have the same thickness in this
         ! approximation.
         fyarea    =max(0.,fyarea -   c_l)

         ! increase in ridged area. min is needed in case ridged ice is
         ! re-ridged. (In this case ridged ice is not conserved)
         ridgearea =min(1.,ridgearea + c_l/kridge)
#endif         
         ficem =fice_max
      else where (ficem<fice_min .and. ficem>0.)
#if defined (TEST_ICE_AGE)
         fyarea   =fyarea  *   (1. + fice_min-ficem)
         ridgearea=ridgearea * (1. + fice_min-ficem)
#endif          
         ficem=fice_min
      end where

      ! NB - Should check the logic here
      where (ficem >= fice_min )
         hicem=icevol/ficem
         hsnwm=snwvol/ficem
#if defined(ALBSNW_EVOL)
         albsnwm=albsnw/ficem
#endif  
#if defined(SSNOWD_ICE)
         hmelt =max(hmlt/ficem,0.)
         hprcp =max(hpcp/ficem,0.)
#endif         
      end where
      
#if defined (TEST_ICE_AGE)
      fy_frac=fyarea
      rdg_frac=ridgearea
      where (fy_frac>1e-3)
         fy_age = fyage/fy_frac
      elsewhere
         fy_age = fyageold
      end where
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Final adjustment of fields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Adjustment for small numbers after diffusion
      !where (ficem<fice_min .and. ficem>fice_tol)
      !   ficem=fice_min
      !elsewhere (ficem<=fice_tol)
      !   ficem=0.
      !   hsnwm=0.
      !   hicem=0.
      !end where
      where (ficem<fice_min .and. ficem>0.)
         ficem=fice_min
      end where
      where (ficem >= fice_min )
         hicem=icevol/ficem
         hsnwm=snwvol/ficem
#if defined(ALBSNW_EVOL)
         albsnwm=albsnw/ficem
#endif
#if defined(SSNOWD_ICE)
         hmelt =max(hmlt/ficem,0.)
         hprcp =max(hpcp/ficem,0.)
#endif         
      end where

      ! Tile final variables
      call xctilr(ficem( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(hicem( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(hsnwm( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
#if defined (TEST_ICE_AGE)
      call xctilr(fy_age  ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(fy_frac ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(rdg_frac( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
#endif  
#if defined(ALBSNW_EVOL)
      call xctilr(albsnwm( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
#endif 
#if defined(SSNOWD_ICE)
      call xctilr(hmelt  ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(hprcp  ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
#endif     

#endif      
#elif defined (ICESTATE)
c --- -------------------------------------------------------------------
c --- ICESTATE ice model advection + diffusion follows.
c --- -------------------------------------------------------------------
c
! APPROACH 1 and 2 : uses MPADATA. Do not work with ICESTATE.
! APPROACH 3 : WENO. Must be used with ICESTATE 
#define ADV_APPROACH_3
#if defined (ADV_APPROACH_1)
      vtp=0.0
      do hk=1,nthick


c ---    NB:  In this approach  we use ice concentration as ``layer'' variable in 
c ---    the HYCOM advection routines.
c
         margin=nbdy
C$OMP PARALLEL DO PRIVATE(j,l,i) 
C$OMP&SCHEDULE(STATIC,jblk)
         do j=1-nbdy,jj+nbdy
         do i=1-nbdy,ii+nbdy
c ---       Used in advection
            farea0(i,j) =icestate(i,j)%ice(hk)%fice
            icevol0(i,j)=icestate(i,j)%ice(hk)%hice*
     &                   icestate(i,j)%ice(hk)%fice
            snwvol0(i,j)=icestate(i,j)%ice(hk)%hsnw*
     &                   icestate(i,j)%ice(hk)%fice
c ---       Conservation of snow mass
            rhosnw(i,j)=icestate(i,j)%ice(hk)%rhosnw
c ---       Conservation of energy
            vtp   (i,j,1:nlay(hk))=
     &         icestate(i,j)%ice(hk)%vtp(1:nlay(hk))
            qbrine(i,j)=icestate(i,j)%ice(hk)%qstore
c ---       Properties of the surface state. 
            albs  (i,j)=icestate(i,j)%ice(hk)%albs
            tsrf  (i,j)=icestate(i,j)%ice(hk)%tsrf
         end do
         end do
C$OMP END PARALLEL DO

  
c ---    Calculate fluxes needed for ice concentration, ice volume and
c ---    snow volume - related advecton.
         margin=nbdy-1 ! due to i+1,i-1
C$OMP PARALLEL DO PRIVATE(j,l,i) 
C$OMP&SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
C
c ---       Set up u-fluxes
            do l=1,isu(j)
            do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
               iuflx(i,j)=iceu(i,j)*scuy(i,j)
               c_iuflx(i,j)=iceu(i,j)*scuy(i,j)*farea0 (i-1,j)
               v_iuflx(i,j)=iceu(i,j)*scuy(i,j)*icevol0(i-1,j)
               s_iuflx(i,j)=iceu(i,j)*scuy(i,j)*snwvol0(i-1,j)
            end do
            end do
c ---       Set up v-fluxes
            do l=1,isv(j)
            do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
               ivflx(i,j)=icev(i,j)*scvx(i,j)
               c_ivflx(i,j)=icev(i,j)*scvx(i,j)*farea0(i,j-1)
               v_ivflx(i,j)=icev(i,j)*scvx(i,j)*icevol0(i,j-1)
               s_ivflx(i,j)=icev(i,j)*scvx(i,j)*snwvol0(i,j-1)
            end do
            end do
         end do

c ---    No need to limit fluxes here, done in adv schemes. Probably need
c ---    explicit diffusion at some stage.

!c --- Set up before and after fake layer thickness 
!         margin=nbdy-1 ! due to i+1,i-1
!C$OMP PARALLEL DO PRIVATE(j,l,i,flxdiv,offset) 
!C$OMP&SCHEDULE(STATIC,jblk)
!         do j=1-margin,jj+margin
!         do l=1,isp(j)
!         do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
!            flxdiv=((c_iuflx(i+1,j)  -c_iuflx(i,j)  )
!     &             +(c_ivflx(i,j+1)  -c_ivflx(i,j)  ))*baclin*scp2i(i,j)
!            dummy2(i,j)=dummy1(i,j)-flxdiv
!            offset=min(0.,dummy1(i,j),dummy2(i,j))
!            dummy2(i,j)=dummy2(i,j)-offset
!            dummy1(i,j)=dummy1(i,j)-offset
!         end do
!         end do
!         end do
!C$OMP END PARALLEL DO
         dummy1=1.
         dummy2=1.

         call xctilr(farea0  (1-nbdy,1-nbdy  ),1,       1, 6,6, halo_ps)
         call xctilr(icevol0 (1-nbdy,1-nbdy  ),1,       1, 6,6, halo_ps)
         call xctilr(snwvol0 (1-nbdy,1-nbdy  ),1,       1, 6,6, halo_ps)

c ---    Set up before and after ice thickness and ice volume
         farea =farea0
         icevol=icevol0
         snwvol=snwvol0

c ---    1) Advect the main conserved variables.
         call advem(1,snwvol,snwvol,c_iuflx,c_ivflx,dummy1,dummy2,
     &              0.,scp2,scp2i,baclin )
         call advem(1,icevol,icevol,c_iuflx,c_ivflx,dummy1,dummy2,
     &              0.,scp2,scp2i,baclin )
         call advem(1,farea,farea  ,c_iuflx,c_ivflx,dummy1,dummy2,
     &              0.,scp2,scp2i,baclin )


c ---    2) Advect different fields using conserved variables as ``layer thickness''
         call xctilr(farea  (1-nbdy,1-nbdy  ),1,       1, 6,6, halo_ps)
         call xctilr(icevol (1-nbdy,1-nbdy  ),1,       1, 6,6, halo_ps)
         call xctilr(snwvol (1-nbdy,1-nbdy  ),1,       1, 6,6, halo_ps)
         call xctilr(rhosnw (1-nbdy,1-nbdy  ),1,       1, 6,6, halo_ps)
         call xctilr(qbrine (1-nbdy,1-nbdy  ),1,       1, 6,6, halo_ps)
         call xctilr(albs   (1-nbdy,1-nbdy  ),1,       1, 6,6, halo_ps)
         call xctilr(tsrf   (1-nbdy,1-nbdy  ),1,       1, 6,6, halo_ps)
         call xctilr(vtp    (1-nbdy,1-nbdy,1),1,nlay(hk), 6,6, halo_ps)
c
c ---    Conservation of snow mass. "dp" is snow volume
         call advem(1,rhosnw,rhosnw,s_iuflx,s_ivflx,
     &              snwvol0,snwvol,0.,scp2,scp2i,baclin )
c ---    Conservation of brine store energy. "dp" is ice volume (= ice mass)
         call advem(1,qbrine,qbrine,v_iuflx,v_ivflx,
     &              icevol0,icevol,0.,scp2,scp2i,baclin )
c ---    "Conservation" of surface albedo. "dp" is ice area fraction
         call advem(1,albs  ,albs  ,c_iuflx,c_ivflx,
     &              farea0,farea,0.,scp2,scp2i,baclin )
c ---    "Conservation" of surface temp. "dp" is ice area fraction
         call advem(1,tsrf  ,tsrf  ,c_iuflx,c_ivflx,
     &              farea0,farea,0.,scp2,scp2i,baclin )
         do hl=1,nlay(hk)
c ---       "Conservation" of ice temp. "dp" is ice volume
            call advem(1,vtp(1,1,hl)  ,vtp(1,1,hl) ,v_iuflx,v_ivflx,
     &                 icevol0,icevol,0.,scp2,scp2i,baclin )
         end do

c ---    Convert from physically conserved quantities to icestate fields
C$OMP PARALLEL DO PRIVATE(j,l,i) 
C$OMP&SCHEDULE(STATIC,jblk)
         do j=1-nbdy,jj+nbdy
         do i=1-nbdy,ii+nbdy
            if (farea(i,j)>epsil0 .and. icevol(i,j)>epsil0) then
               icestate(i,j)%ice(hk)%fice=farea (i,j)
c ---          Properties of the surface state -
               icestate(i,j)%ice(hk)%albs=albs(i,j)
               icestate(i,j)%ice(hk)%tsrf=tsrf(i,j)
c ---          Conservation of volume fields
               icestate(i,j)%ice(hk)%hice=icevol(i,j)/farea(i,j)
               icestate(i,j)%ice(hk)%hsnw=snwvol(i,j)/farea(i,j)
c ---          Conservation of mass fields
               icestate(i,j)%ice(hk)%rhosnw=rhosnw(i,j)
c ---          Conservation of energy fields
               icestate(i,j)%ice(hk)%qstore=qbrine(i,j)
               icestate(i,j)%ice(hk)%vtp(1:nlay(hk))= 
     &            vtp   (i,j,1:nlay(hk))
            else
               icestate(i,j)%ice(hk)=clear_ice(icestate(i,j)%ice(hk),
     &                                         nlay(hk))
            end if
         end do
         end do
C$OMP END PARALLEL DO
      end do !hk



#elif defined (ADV_APPROACH_2)
      vtp=0.0
      do hk=1,nthick


C$OMP PARALLEL DO PRIVATE(j,l,i) 
C$OMP&SCHEDULE(STATIC,jblk)
         do j=1-nbdy,jj+nbdy
         do i=1-nbdy,ii+nbdy
            farea (i,j)=icestate(i,j)%ice(hk)%fice*tenm/dummy1(i,j)
c ---       Conservation of volume
            icevol(i,j)=icestate(i,j)%ice(hk)%hice*
     &                  icestate(i,j)%ice(hk)%fice*tenm/dummy1(i,j)
            snwvol(i,j)=icestate(i,j)%ice(hk)%hsnw*
     &                  icestate(i,j)%ice(hk)%fice*tenm/dummy1(i,j)
c ---       Conservation of mass
            rhosnw(i,j)=icestate(i,j)%ice(hk)%rhosnw*
     &                  icestate(i,j)%ice(hk)%fice*
     &                  icestate(i,j)%ice(hk)%hsnw*
     &                  tenm/dummy1(i,j)
c ---       Conservation of energy
            vtp   (i,j,1:nlay(hk))=
     &                  icestate(i,j)%ice(hk)%hice*
     &                  icestate(i,j)%ice(hk)%fice*
     &                  icestate(i,j)%ice(hk)%vtp(1:nlay(hk))* 
     &                  tenm/dummy1(i,j)
            qbrine(i,j)=icestate(i,j)%ice(hk)%fice*
c     &                  icestate(i,j)%ice(hk)%hice*
     &                  icestate(i,j)%ice(hk)%qstore*
     &                  tenm/dummy1(i,j)
c ---       Properties of the surface state. 
            albs  (i,j)=icestate(i,j)%ice(hk)%albs*
     &                  icestate(i,j)%ice(hk)%fice*tenm/dummy1(i,j)
            tsrf  (i,j)=icestate(i,j)%ice(hk)%tsrf*
     &                  icestate(i,j)%ice(hk)%fice*tenm/dummy1(i,j)
         end do
         end do
C$OMP END PARALLEL DO

c ---    Tile  and advect variables. For now use the mpdata routine - 
c ---    which does not use middle time step. After advection "integrate"
c ---    field to recover 
         margin=nbdy 
         call xctilr(farea  (1-nbdy,1-nbdy  ),1,       1, 6,6, halo_ps)
         call xctilr(icevol (1-nbdy,1-nbdy  ),1,       1, 6,6, halo_ps)
         call xctilr(snwvol (1-nbdy,1-nbdy  ),1,       1, 6,6, halo_ps)
         call xctilr(rhosnw (1-nbdy,1-nbdy  ),1,       1, 6,6, halo_ps)
         call xctilr(qbrine (1-nbdy,1-nbdy  ),1,       1, 6,6, halo_ps)
         call xctilr(albs   (1-nbdy,1-nbdy  ),1,       1, 6,6, halo_ps)
         call xctilr(tsrf   (1-nbdy,1-nbdy  ),1,       1, 6,6, halo_ps)
         call xctilr(vtp    (1-nbdy,1-nbdy,1),1,nlay(hk), 6,6, halo_ps)

         call advem(1,farea,farea,iuflx,ivflx,
     &              dummy1,dummy2,0.,scp2,scp2i,baclin )
         call advem(1,icevol,icevol,iuflx,ivflx,
     &              dummy1,dummy2,0.,scp2,scp2i,baclin )
         call advem(1,snwvol,snwvol,iuflx,ivflx,
     &              dummy1,dummy2,0.,scp2,scp2i,baclin )
         call advem(1,rhosnw,rhosnw,iuflx,ivflx,
     &              dummy1,dummy2,0.,scp2,scp2i,baclin )
         call advem(1,qbrine,qbrine,iuflx,ivflx,
     &              dummy1,dummy2,0.,scp2,scp2i,baclin )
         call advem(1,albs  ,albs  ,iuflx,ivflx,
     &              dummy1,dummy2,0.,scp2,scp2i,baclin )
         call advem(1,tsrf  ,tsrf  ,iuflx,ivflx,
     &              dummy1,dummy2,0.,scp2,scp2i,baclin )
         do hl=1,nlay(hk)
            call advem(1,vtp(1,1,hl)  ,vtp(1,1,hl) ,iuflx,ivflx,
     &                 dummy1,dummy2,0.,scp2,scp2i,baclin )
         end do
c
c ---    Recover from "layer-weighted" to physical 
         farea =farea *dummy2/tenm
         icevol=icevol*dummy2/tenm
         snwvol=snwvol*dummy2/tenm
         rhosnw=rhosnw*dummy2/tenm
         qbrine=qbrine*dummy2/tenm
         albs  =albs  *dummy2/tenm
         tsrf  =tsrf  *dummy2/tenm
         do hl=1,nlay(hk)
            vtp(:,:,hl)   =vtp(:,:,hl)   *dummy2/tenm
         end do
c ---    Add diffusion at this point. Skipped for now
c
c ---    Convert from physically conserved quantities to icestate fields
C$OMP PARALLEL DO PRIVATE(j,l,i) 
C$OMP&SCHEDULE(STATIC,jblk)
         do j=1-nbdy,jj+nbdy
         do i=1-nbdy,ii+nbdy
            icestate(i,j)%ice(hk)%fice=farea (i,j)
            if (farea(i,j)>epsil0 .and. icevol(i,j)>epsil0) then
c ---          Properties of the surface state -
               icestate(i,j)%ice(hk)%albs=albs(i,j)/farea(i,j)
               icestate(i,j)%ice(hk)%tsrf=tsrf(i,j)/farea(i,j)
c ---          Conservation of volume fields
               icestate(i,j)%ice(hk)%hice=icevol(i,j)/farea(i,j)
               icestate(i,j)%ice(hk)%hsnw=snwvol(i,j)/farea(i,j)
               if (snwvol(i,j)>epsil0) then
c ---             Conservation of mass fields
                   icestate(i,j)%ice(hk)%rhosnw=rhosnw(i,j)/snwvol(i,j)
               else
                   icestate(i,j)%ice(hk)%rhosnw=rhosnwmin
               end if
c ---          Conservation of energy fields
               icestate(i,j)%ice(hk)%vtp(1:nlay(hk))=
     &           vtp   (i,j,1:nlay(hk))/icevol(i,j)
c               icestate(i,j)%ice(hk)%qstore=qbrine(i,j)/icevol(i,j)
                icestate(i,j)%ice(hk)%qstore=qbrine(i,j)/farea(i,j)
            else
               icestate(i,j)%ice(hk)=
     &            clear_ice(icestate(i,j)%ice(hk),nlay(hk))
            end if
         end do
         end do
C$OMP END PARALLEL DO
      end do !hk
#elif defined (ADV_APPROACH_3)
c --- Ice variables are now initialized on the full margin - we tile later
      margin=nbdy
C$OMP PARALLEL DO PRIVATE(j,l,i) 
C$OMP&SCHEDULE(STATIC,jblk)
      
      do j=1-nbdy,jdm+nbdy
         do i=1-nbdy,idm+nbdy
            icevol(i,j)=0.
            snwvol(i,j)=0.
            farea (i,j)=0.
            rhosnw(i,j)=0.
            qbrine(i,j)=0.
            albs  (i,j)=0.
            tsrf  (i,j)=0.
#if defined (SSNOWD)            
            melt  (i,j)=0.
            prcp  (i,j)=0.
#endif   
#if defined (ICEAGE)
            age   (i,j)=0.
#endif            
         end do
      end do

      vtp=0.0
      do hk=1,nthick


C$OMP PARALLEL DO PRIVATE(j,l,i) 
C$OMP&SCHEDULE(STATIC,jblk)
         do j=1-nbdy,jj+nbdy
         do i=1-nbdy,ii+nbdy
            farea (i,j)=icestate(i,j)%ice(hk)%fice
c ---       Conservation of volume
            icevol(i,j)=icestate(i,j)%ice(hk)%hice*
     &                  icestate(i,j)%ice(hk)%fice
            snwvol(i,j)=icestate(i,j)%ice(hk)%hsnw*
     &                  icestate(i,j)%ice(hk)%fice
c ---       Conservation of mass
            rhosnw(i,j)=icestate(i,j)%ice(hk)%rhosnw*
     &                  icestate(i,j)%ice(hk)%fice*
     &                  icestate(i,j)%ice(hk)%hsnw
c ---       Conservation of energy
            vtp   (i,j,1:nlay(hk))=
     &                  icestate(i,j)%ice(hk)%hice*
     &                  icestate(i,j)%ice(hk)%fice*
     &                  icestate(i,j)%ice(hk)%vtp(1:nlay(hk)) 
            qbrine(i,j)=icestate(i,j)%ice(hk)%fice*
c     &                  icestate(i,j)%ice(hk)%hice*
     &                  icestate(i,j)%ice(hk)%qstore
c ---       Properties of the surface state. 
            albs  (i,j)=icestate(i,j)%ice(hk)%albs*
     &                  icestate(i,j)%ice(hk)%fice
            tsrf  (i,j)=icestate(i,j)%ice(hk)%tsrf*
     &                  icestate(i,j)%ice(hk)%fice
#if defined (SSNOWD)
            prcp  (i,j)=icestate(i,j)%ice(hk)%hprcp*
     &                  icestate(i,j)%ice(hk)%fice
            melt  (i,j)=icestate(i,j)%ice(hk)%hmelt*
     &                  icestate(i,j)%ice(hk)%fice
#endif 
#if defined (ICEAGE)
            age(i,j)  = icestate(i,j)%ice(hk)%fice*
     &                  icestate(i,j)%ice(hk)%hice*
     &                  icestate(i,j)%ice(hk)%age
#endif
         end do
         end do
C$OMP END PARALLEL DO

c ---    Tile  and advect variables. Use the WENO routine - 
c ---    which does not use middle time step. After advection "integrate"
c ---    field to recover 
         margin=nbdy 
         call xctilr(farea  (1-nbdy,1-nbdy  ),1,       1, 6,6, halo_ps)
         call xctilr(icevol (1-nbdy,1-nbdy  ),1,       1, 6,6, halo_ps)
         call xctilr(snwvol (1-nbdy,1-nbdy  ),1,       1, 6,6, halo_ps)
         call xctilr(rhosnw (1-nbdy,1-nbdy  ),1,       1, 6,6, halo_ps)
         call xctilr(qbrine (1-nbdy,1-nbdy  ),1,       1, 6,6, halo_ps)
         call xctilr(albs   (1-nbdy,1-nbdy  ),1,       1, 6,6, halo_ps)
         call xctilr(tsrf   (1-nbdy,1-nbdy  ),1,       1, 6,6, halo_ps)
         call xctilr(vtp    (1-nbdy,1-nbdy,1),1,nlay(hk), 6,6, halo_ps)
#if defined (SSNOWD)
         call xctilr(prcp   (1-nbdy,1-nbdy  ),1,       1, 6,6, halo_ps)
         call xctilr(melt   (1-nbdy,1-nbdy  ),1,       1, 6,6, halo_ps)
#endif  
#if defined (ICEAGE)
         call xctilr(age    (1-nbdy,1-nbdy  ),1,       1, 6,6, halo_ps)
#endif         
c --- ------------------------------------------------------------------
c --- Advection is done with flux limited 3rd order WENO in space and
c --- 2nd order Runge-Kutta in time
c --- ------------------------------------------------------------------

         call iceadv(farea ,iceu,icev,scuy,scvx,scp2i,scp2,baclin)
         call iceadv(icevol,iceu,icev,scuy,scvx,scp2i,scp2,baclin)
         call iceadv(snwvol,iceu,icev,scuy,scvx,scp2i,scp2,baclin)
         call iceadv(rhosnw,iceu,icev,scuy,scvx,scp2i,scp2,baclin)
         call iceadv(qbrine,iceu,icev,scuy,scvx,scp2i,scp2,baclin)
         call iceadv(albs  ,iceu,icev,scuy,scvx,scp2i,scp2,baclin)
         call iceadv(tsrf  ,iceu,icev,scuy,scvx,scp2i,scp2,baclin)
         do hl=1,nlay(hk)
             call iceadv(vtp(:,:,hl),iceu,icev,scuy,scvx,scp2i,         &
     &         scp2,baclin)
         end do
#if defined (SSNOWD)
         call iceadv(prcp  ,iceu,icev,scuy,scvx,scp2i,scp2,baclin)
         call iceadv(melt  ,iceu,icev,scuy,scvx,scp2i,scp2,baclin)
#endif         
#if defined (ICEAGE)
         call iceadv(age   ,iceu,icev,scuy,scvx,scp2i,scp2,baclin)
#endif
c ---    Convert from physically conserved quantities to icestate fields
C$OMP PARALLEL DO PRIVATE(j,l,i) 
C$OMP&SCHEDULE(STATIC,jblk)
         do j=1-nbdy,jj+nbdy
         do i=1-nbdy,ii+nbdy
            icestate(i,j)%ice(hk)%fice=farea (i,j)
            if (farea(i,j)>epsil0 .and. icevol(i,j)>epsil0) then
c ---          Properties of the surface state -
               icestate(i,j)%ice(hk)%albs=albs(i,j)/farea(i,j)
               icestate(i,j)%ice(hk)%tsrf=tsrf(i,j)/farea(i,j)
c ---          Conservation of volume fields
               icestate(i,j)%ice(hk)%hice=icevol(i,j)/farea(i,j)
               icestate(i,j)%ice(hk)%hsnw=snwvol(i,j)/farea(i,j)
               if (snwvol(i,j)>epsil0) then
c ---             Conservation of mass fields
                   icestate(i,j)%ice(hk)%rhosnw=rhosnw(i,j)/snwvol(i,j)
               else
                  icestate(i,j)%ice(hk)%rhosnw=rhosnwmin
               end if
c ---          Conservation of energy fields
               icestate(i,j)%ice(hk)%vtp(1:nlay(hk))=
     &           vtp   (i,j,1:nlay(hk))/icevol(i,j)
c               icestate(i,j)%ice(hk)%qstore=qbrine(i,j)/icevol(i,j)
               icestate(i,j)%ice(hk)%qstore=qbrine(i,j)/farea(i,j)
#if defined (SSNOWD)
               icestate(i,j)%ice(hk)%hprcp=prcp(i,j)/farea(i,j)
               !max avoids apparition of very small value for hmelt
               ! (for eks. : -1e-38)
               icestate(i,j)%ice(hk)%hmelt=max(0.,melt(i,j)/farea(i,j))
#endif    
#if defined(ICEAGE)
               icestate(i,j)%ice(hk)%age = min(age_max,                  &
     &                max(0.,age(i,j)/icevol(i,j)))
#endif               
            else
#if defined(SSNOWD)               
               icestate(i,j)%ice(hk)=
     &            clear_ice(icestate(i,j)%ice(hk),nlay(hk),cvsnw(hk))
#else
               icestate(i,j)%ice(hk)=
     &            clear_ice(icestate(i,j)%ice(hk),nlay(hk))
#endif               
            end if
         end do
         end do
C$OMP END PARALLEL DO
      end do !hk

#else
#error - No advection method defined for ICESTATE in icemodels_advect
#endif
c
c --- Mechanical redistribution done in mod_icestate_mechred. 
c --- sum(fice) may be > 1 at this point
#endif /* ICE or ICESTATE */
      
      ! Restore old margin
      margin=oldmargin


      end subroutine icemodels_advect
      end module m_icemodels_advect
!!***

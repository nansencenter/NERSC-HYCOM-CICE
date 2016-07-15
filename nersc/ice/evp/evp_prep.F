!=======================================================================
!BOP
!
! !IROUTINE: evp_prep - compute quantities needed for stress tensor and mom eqns
!
! !INTERFACE:
!
      subroutine evp_prep(kstrngth)
!
! !DESCRIPTION:
!
! Computes quantities needed in the stress tensor (sigma)  \\
! and momentum (u) equations, but which do not change during  \\
! the thermodynamics/transport time step: \\
! --wind stress shift to U grid, \\
! --ice mass and ice extent masks, \\
! --pressure (strength), and part of the forcing stresses \\
! initializes ice velocity for new points to ocean sfc current  \\
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke
!         Fluid Dynamics Group, Los Alamos National Laboratory
!
! !USES:
!
      use mod_evp
#if defined(ICE_DYN_DIAG)
       use mod_common_ice, only : strainI, strainII
#endif
      implicit none
!cice      use ice_mechred_cice
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer, intent(in) ::
     &   kstrngth    
!
!EOP
!
      integer :: i, j, k, n


      real :: 
     &   umass(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ! ice mass on u-grid
     &,  tmp  (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
     &,  tmp2 (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)

      real, parameter ::
     &   a_min = .001  ! minimum ice area
     &,  m_min = .01   ! minimum ice mass

      logical ::
     &   tmphm       (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
     &,  iceumask_old(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)

      real, parameter ::
     &   rhoi=930.0  ! density of ice
     &  ,rhos=250.0  ! density of snow

      real :: wabs

      !-----------------------------------------------------------------
      ! total mass of ice and snow, centered in T-cell
      ! NOTE: vice and vsno must be up to date in all grid cells,
      !       including ghost cells
      !-----------------------------------------------------------------
      imargin = nbdy
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-imargin,jj+imargin
       do i=1-imargin,ii+imargin
        if (tmask(i,j)) then
         tmass(i,j) = (rhoi*vice(i,j) + rhos*vsno(i,j)) ! kg/m^2
        else
         tmass(i,j) = 0.0
        endif
       enddo
      enddo
!$OMP END PARALLEL DO 
      imargin = nbdy-1
      call to_ugrid(tmass,umass,imargin)
      call to_ugrid(aice,aiu,imargin)
      call xctilr(umass  ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_qs)
      call xctilr(aiu    ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_qs)

      !-----------------------------------------------------------------
      ! convert dynamics variables to U grid
      !-----------------------------------------------------------------
      imargin = nbdy
!$OMP PARALLEL DO PRIVATE(j,i,wabs)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-imargin,jj+imargin
       do i=1-imargin,ii+imargin
        ! Factor of aice needed for correct treatment of free drift
        !strairx(i,j) = strairxT(i,j)*aice(i,j) ! prep to convert to U grid
        !strairy(i,j) = strairyT(i,j)*aice(i,j)

        wabs         = draga*sqrt(uair(i,j)**2+vair(i,j)**2)
CKAL    strairx(i,j) = wabs*(cosa*uair(i,j)-sina*vair(i,j))
CKAL    strairy(i,j) = wabs*(sina*uair(i,j)+cosa*vair(i,j))
CKAL -- NB
        strairx(i,j) = wabs*(cosa*uair(i,j)-sina*vair(i,j))*aiu(i,j)
        strairy(i,j) = wabs*(sina*uair(i,j)+cosa*vair(i,j))*aiu(i,j)
       enddo
      enddo
!$OMP END PARALLEL DO 
      ! KAL TODO - fix these routines
      !call t2ugrid(strairx) -- KAL int. to v-grid in hycom2evp
      !call t2ugrid(strairy) -- KAL int. to v-grid in hycom2evp
      imargin=nbdy-1
      !call xctilr(umass  ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_qs)
      !call xctilr(aiu    ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_qs)
      call xctilr(strairx( 1-nbdy,1-nbdy),1,   1, 6,6, halo_uv)
      call xctilr(strairy( 1-nbdy,1-nbdy),1,   1, 6,6, halo_vv)


      !-----------------------------------------------------------------
      ! convenient variable for evp
      !-----------------------------------------------------------------
      imargin = nbdy
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-imargin,jj+imargin
       do i=1-imargin,ii+imargin
        umassdtei(i,j) = umass(i,j)*dtei               ! m/dte, kg/m^2 s
       enddo
      enddo
!$OMP END PARALLEL DO 

      !-----------------------------------------------------------------
      ! augmented masks (land + open ocean)
      !-----------------------------------------------------------------
      imargin = nbdy
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-imargin,jj+imargin
       do i=1-imargin,ii+imargin
        ! ice extent mask (T-cells)
        tmphm(i,j) = tmask(i,j) .and. (aice (i,j).gt.a_min)
     &                          .and. (tmass(i,j).gt.m_min)
       enddo
      enddo
!$OMP END PARALLEL DO 

      tmp=0.
      tmp2=0.
      icetmask=.false.
      iceumask_old=.false.
      imargin = nbdy-1 !NB
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-imargin,jj+imargin
       do i=1-imargin,ii+imargin
        ! extend ice extent mask (T-cells) to points around pack
        icetmask(i,j) =  
     &   tmphm(i-1,j+1) .or. tmphm(i,j+1) .or. tmphm(i+1,j+1) .or.
     &   tmphm(i-1,j)   .or. tmphm(i,j)   .or. tmphm(i+1,j)   .or.
     &   tmphm(i-1,j-1) .or. tmphm(i,j-1) .or. tmphm(i+1,j-1)
        icetmask(i,j) = icetmask(i,j) .and. tmask(i,j) ! remask land points
        ! ice extent mask (U-cells)
         iceumask_old(i,j) = iceumask(i,j)  ! save
         iceumask(i,j) = (umask(i,j)) .and. (aiu  (i,j).gt.a_min)   
     &                                .and. (umass(i,j).gt.m_min)
         if (iceumask(i,j)) tmp (i,j)=1.0
         if (icetmask(i,j)) tmp2(i,j)=1.0
       enddo
      enddo
!$OMP END PARALLEL DO 

      call xctilr(tmp   ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_qs)
      call xctilr(tmp2  ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      imargin=nbdy
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-imargin,jj+imargin
      do i=1-imargin,ii+imargin
         iceumask(i,j)=tmp (i,j)>.5
         icetmask(i,j)=tmp2(i,j)>.5
      end do
      end do
!$OMP END PARALLEL DO 

      !-----------------------------------------------------------------
      ! pressure and forcing terms; set sigma=0 for no ice;
      ! initialize ice velocity in cells previously empty to ocn current
      !-----------------------------------------------------------------
!cice      call ice_strength(kstrngth)
      call evp_ice_strength(kstrngth)

      imargin = nbdy-1 !NB
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-imargin,jj+imargin
       do i=1-imargin,ii+imargin
       if (icetmask(i,j)) then
        prss(i,j) = strength(i,j)
        fm(i,j) = fcor(i,j)*umass(i,j)               ! Coriolis * mass

        if (umask(i,j)) then
        ! for ocean stress
        waterx(i,j) = uocn(i,j)*cosw - vocn(i,j)*sinw
        watery(i,j) = vocn(i,j)*cosw + uocn(i,j)*sinw
        ! combine tilt with wind stress
!#ifndef coupled
!        ! calculate tilt from geostrophic currents if needed
!        strtltx(i,j) = -fm(i,j)*vocn(i,j)
!        strtlty(i,j) =  fm(i,j)*uocn(i,j)
!#else
        strtltx(i,j) = -gravit*umass(i,j)*ss_tltx(i,j)
        strtlty(i,j) = -gravit*umass(i,j)*ss_tlty(i,j)
!#endif
        forcex(i,j) = strairx(i,j) + strtltx(i,j)
        forcey(i,j) = strairy(i,j) + strtlty(i,j)
        endif  ! umask

       else    ! .not. icetmask
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

        divu(i,j) = 0.0
        Delta(i,j) = 0.0
        shear(i,j) = 0.0
#if defined(ICE_DYN_DIAG)
        strainI(i,j)  = 0.0
        strainII(i,j) = 0.0
#endif
       endif  ! icetmask

        ! initialize velocity for new ice points to ocean sfc current
        if( iceumask(i,j) .and. (.not. iceumask_old(i,j))) then
          uvel(i,j) = uocn(i,j)
          vvel(i,j) = vocn(i,j)
CKAL -- Reset velocities to zero outside of iceumask
        else if (.not.iceumask(i,j)) then
           uvel(i,j)=0.0
           vvel(i,j)=0.0
        endif 

       enddo
      enddo
!$OMP END PARALLEL DO 
      call xctilr(prss   ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(fm     ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_qs)
      call xctilr(strtltx( 1-nbdy,1-nbdy),1,   1, 6,6, halo_qs)
      call xctilr(strtlty( 1-nbdy,1-nbdy),1,   1, 6,6, halo_qs)
      call xctilr(forcex ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_qs)
      call xctilr(forcey ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_qs)
      call xctilr(uvel   ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_qs)
      call xctilr(vvel   ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_qs)

      if (e_itst>-99 .and. e_jtst >-99) then
         print *,imargin
         print *,'evp_prep',aice(e_itst,e_jtst)
         print *,'evp_prep',vice(e_itst,e_jtst)
         print *,'evp_prep',tmass(e_itst,e_jtst)
         print *,'evp_prep umass',umass(e_itst,e_jtst)
         print *,'evp_prep fcor',fcor(e_itst,e_jtst)
         print *,'evp_prep',fm(e_itst,e_jtst)
         print *,'evp_prep',strtltx(e_itst,e_jtst)
         print *,'evp_prep',strtlty(e_itst,e_jtst)
         print *,'evp_prep',forcex(e_itst,e_jtst)
         print *,'evp_prep',forcey(e_itst,e_jtst)
         print *,'evp_prep',uvel(e_itst,e_jtst)
         print *,'evp_prep',vvel(e_itst,e_jtst)
      end if

      !call xcstop ('(evp_pre)')

CKAL  ! KAL - For vectorized version - not used for now
CKAL  imargin = nbdy-1 !Reduce for initial derivation in "stress"
CKAL  icellt = 0
CKAL  do j=1-imargin,jj+imargin
CKAL   do i=1-imargin,ii+imargin
CKAL      if (icetmask(i,j)) then
CKAL        icellt = icellt + 1
CKAL        indxti(icellt) = i
CKAL        indxtj(icellt) = j
CKAL      endif
CKAL    enddo
CKAL  enddo
CKAL
CKAL  icellu = 0
CKAL  do j=1-imargin,jj+imargin
CKAL   do i=1-imargin,ii+imargin
CKAL      if (iceumask(i,j)) then
CKAL        icellu = icellu + 1
CKAL        indxui(icellu) = i
CKAL        indxuj(icellu) = j
Ctst        if (i==93.and.j==58) then
Ctst           print *,'i j test',iceumask(i,j),umass(i,j),mnproc
Ctst        end if
CKAL      endif
CKAL    enddo
CKAL  enddo
      end subroutine evp_prep

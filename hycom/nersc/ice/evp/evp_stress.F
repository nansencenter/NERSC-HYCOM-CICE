!!****f* HYCOM_2.1.03/EVP_MPI/evp_stress
!!
!! NAME 
!!    stress 
!!
!! SYNOPSIS
!!    
!!
!!
!! DESCRIPTION
!!    Computes strain rates and internal stress components.
!!    Computes the rates of strain and internal stress components for 
!!    each of the four corners on each T-grid cell
!!
!!
!! AUTHOR
!!    Elizabeth C. Hunke, Fluid Dynamics Group, Los Alamos National Laboratory
!!
!! HISTORY
!!    Oct 2006   - Modified to use HYCOM tiling logic - Knut Liseter
!!    12.03.2007 - Added OMP directives
!!    31.07.2014 - Modified according to Bouillon et al. OM 2013
!!
!! INPUT
!!    ksub -- Subcycling time step
!! SIDE EFFECTS
!!    Sets up the following variables for evp_stepu
!!       stressm_[1-4]
!!       stressp_[1-4]
!!       stress12_[1-4]
!!    Keeps the following variables for diagnostics and mechanical redist
!!       Divu
!!       Delta
!!       shear
!!    Because a grid difference is involved, the effective margin is
!!    reduced by one
!!
!! SOURCE
!!
      subroutine evp_stress(ksub)
      use mod_evp
#if defined(ICE_DYN_DIAG)
      use mod_common_ice, only : strainI, strainII
#endif
      implicit none
      integer, intent(in) :: ksub  ! subcycling step
      integer :: i, j

      real :: 
     &  divune, divunw, divuse, divusw             ! divergence
     &, tensionne, tensionnw, tensionse, tensionsw ! tension
     &, shearne, shearnw, shearse, shearsw         ! shearing
     &, Deltane, Deltanw, Deltase, Deltasw         ! Delta 
     &, c0ne, c0nw, c0se, c0sw                     ! useful combinations
     &, c1ne, c1nw, c1se, c1sw

      integer ::
     &  ij      ! loop index, combination of i and j loops

Coldcdir$ ivdep      !Cray
Cold!cdir nodep      !NEC
Cold!ocl novrec      !Fujitsu
Cold  do ij=1,icellt
Cold  !KAL - added margin check
Cold   i = indxti(ij)
Cold   j = indxtj(ij)
Cold  if (i.ge.1-imargin .and. i.le.ii+imargin .and.
Cold &    j.ge.1-imargin .and. j.le.jj+imargin ) then

!$OMP PARALLEL DO PRIVATE(j,i,
!$OMP&     divune,    divunw,   divusw,   divuse,
!$OMP&  tensionne, tensionnw,tensionsw,tensionse,
!$OMP&    shearne,   shearnw,  shearsw,  shearse,
!$OMP&    Deltane,   Deltanw,  Deltasw,  Deltase,
!$OMP&       c0ne,      c0nw,     c0sw,     c0se,
!$OMP&       c1ne,      c1nw,     c1sw,     c1se)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-imargin,jj+imargin
      do i=1-imargin,ii+imargin
      if (icetmask(i,j)) then

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------
      ! divergence  =  e_11 + e_22
      divune    = cyp(i,j)*uvel(i  ,j  ) - dyt(i,j)*uvel(i-1,j  )
     &          + cxp(i,j)*vvel(i  ,j  ) - dxt(i,j)*vvel(i  ,j-1)
      divunw    = cym(i,j)*uvel(i-1,j  ) + dyt(i,j)*uvel(i  ,j  )
     &          + cxp(i,j)*vvel(i-1,j  ) - dxt(i,j)*vvel(i-1,j-1)
      divusw    = cym(i,j)*uvel(i-1,j-1) + dyt(i,j)*uvel(i  ,j-1)
     &          + cxm(i,j)*vvel(i-1,j-1) + dxt(i,j)*vvel(i-1,j  )
      divuse    = cyp(i,j)*uvel(i  ,j-1) - dyt(i,j)*uvel(i-1,j-1)
     &          + cxm(i,j)*vvel(i  ,j-1) + dxt(i,j)*vvel(i  ,j  )

      ! tension strain rate  =  e_11 - e_22
      tensionne = -cym(i,j)*uvel(i  ,j  ) - dyt(i,j)*uvel(i-1,j  )
     &          +  cxm(i,j)*vvel(i  ,j  ) + dxt(i,j)*vvel(i  ,j-1)
      tensionnw = -cyp(i,j)*uvel(i-1,j  ) + dyt(i,j)*uvel(i  ,j  )
     &          +  cxm(i,j)*vvel(i-1,j  ) + dxt(i,j)*vvel(i-1,j-1)
      tensionsw = -cyp(i,j)*uvel(i-1,j-1) + dyt(i,j)*uvel(i  ,j-1)
     &          +  cxp(i,j)*vvel(i-1,j-1) - dxt(i,j)*vvel(i-1,j  )
      tensionse = -cym(i,j)*uvel(i  ,j-1) - dyt(i,j)*uvel(i-1,j-1)
     &          +  cxp(i,j)*vvel(i  ,j-1) - dxt(i,j)*vvel(i  ,j  )

      ! shearing strain rate  =  e_12
      shearne = -cym(i,j)*vvel(i  ,j  ) - dyt(i,j)*vvel(i-1,j  )
     &        -  cxm(i,j)*uvel(i  ,j  ) - dxt(i,j)*uvel(i  ,j-1)
      shearnw = -cyp(i,j)*vvel(i-1,j  ) + dyt(i,j)*vvel(i  ,j  )
     &        -  cxm(i,j)*uvel(i-1,j  ) - dxt(i,j)*uvel(i-1,j-1)
      shearsw = -cyp(i,j)*vvel(i-1,j-1) + dyt(i,j)*vvel(i  ,j-1)
     &        -  cxp(i,j)*uvel(i-1,j-1) + dxt(i,j)*uvel(i-1,j  )
      shearse = -cym(i,j)*vvel(i  ,j-1) - dyt(i,j)*vvel(i-1,j-1)
     &        -  cxp(i,j)*uvel(i  ,j-1) + dxt(i,j)*uvel(i  ,j  )

      ! Delta (in the denominator of zeta, eta)
      Deltane = sqrt(divune**2 + ecci*(tensionne**2 + shearne**2))
      Deltanw = sqrt(divunw**2 + ecci*(tensionnw**2 + shearnw**2))
      Deltase = sqrt(divuse**2 + ecci*(tensionse**2 + shearse**2))
      Deltasw = sqrt(divusw**2 + ecci*(tensionsw**2 + shearsw**2))

      !-----------------------------------------------------------------
      ! save quantities for mechanical redistribution
      !-----------------------------------------------------------------
      if (ksub.eq.ndte) then
      divu(i,j) = 0.25*(divune + divunw + divuse + divusw)*tarear(i,j)
      Delta(i,j) = 0.25*(Deltane+Deltanw+Deltase+Deltasw)*tarear(i,j)
      ! diagnostic only
      ! shear = sqrt(tension**2 + shearing**2) 
      shear(i,j) = 0.25*tarear(i,j)*sqrt(
     &  (tensionne + tensionnw + tensionse + tensionsw)**2
     & +(  shearne +   shearnw +   shearse +   shearsw)**2)
#if defined(ICE_DYN_DIAG)
      strainI(i,j)  = 0.125*(divune+divunw+divuse+divusw)*tarear(i,j)
      strainII(i,j) = 0.125*tarear(i,j)*sqrt(
     &  (tensionne + tensionnw + tensionse + tensionsw)**2
     & +(  shearne +   shearnw +   shearse +   shearsw)**2)
#endif
      endif

      !-----------------------------------------------------------------
      ! replacement pressure/Delta                   ! kg/s
      ! save replacement pressure for principal stress calculation
      !-----------------------------------------------------------------
      if (evp_damping) then
        ! enforce damping criterion
        c0ne=min(prss(i,j)/max(Deltane,4.*tinyarea(i,j)),rcon_evp(i,j))
        c0nw=min(prss(i,j)/max(Deltanw,4.*tinyarea(i,j)),rcon_evp(i,j))
        c0sw=min(prss(i,j)/max(Deltasw,4.*tinyarea(i,j)),rcon_evp(i,j))
        c0se=min(prss(i,j)/max(Deltase,4.*tinyarea(i,j)),rcon_evp(i,j))
        prs_sig(i,j) = prss(i,j)*Deltane/max(Deltane,4.*tinyarea(i,j)) ! ne
      else
        ! original version
        c0ne = prss(i,j)/max(Deltane,tinyarea(i,j))
        c0nw = prss(i,j)/max(Deltanw,tinyarea(i,j))
        c0sw = prss(i,j)/max(Deltasw,tinyarea(i,j))
        c0se = prss(i,j)/max(Deltase,tinyarea(i,j))
        prs_sig(i,j) = c0ne*Deltane ! northeast
      endif

      c1ne = c0ne*dte2T
      c1nw = c0nw*dte2T
      c1sw = c0sw*dte2T
      c1se = c0se*dte2T

      !-----------------------------------------------------------------
      ! the stresses                            ! kg/s^2
      ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
      !-----------------------------------------------------------------

      stressp_1(i,j) = (stressp_1(i,j) + c1ne*(divune - Deltane))*denom
      stressp_2(i,j) = (stressp_2(i,j) + c1nw*(divunw - Deltanw))*denom
      stressp_3(i,j) = (stressp_3(i,j) + c1sw*(divusw - Deltasw))*denom
      stressp_4(i,j) = (stressp_4(i,j) + c1se*(divuse - Deltase))*denom

      !Modified EVP according to Bouillon et al. OM 2013 
      stressm_1(i,j) = (stressm_1(i,j) + c1ne*ecci*tensionne)*denom
      stressm_2(i,j) = (stressm_2(i,j) + c1nw*ecci*tensionnw)*denom
      stressm_3(i,j) = (stressm_3(i,j) + c1sw*ecci*tensionsw)*denom
      stressm_4(i,j) = (stressm_4(i,j) + c1se*ecci*tensionse)*denom

      stress12_1(i,j) = (stress12_1(i,j) + c1ne*ecci*shearne*.5)*denom
      stress12_2(i,j) = (stress12_2(i,j) + c1nw*ecci*shearnw*.5)*denom
      stress12_3(i,j) = (stress12_3(i,j) + c1sw*ecci*shearsw*.5)*denom
      stress12_4(i,j) = (stress12_4(i,j) + c1se*ecci*shearse*.5)*denom
 
      endif
      enddo
      enddo
!$OMP END PARALLEL DO 

Cold  endif                    
Cold  enddo                     ! ij

      end subroutine evp_stress
!!******

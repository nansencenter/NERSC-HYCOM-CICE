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
!!    28.07.2009 - Corrected a few bugs and cleaned the code - Dany Dumont
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
      subroutine evp_stress_miz(ksub)
      use mod_evp
#if defined(ICE_DYN_DIAG)
      use mod_common_ice, only : MIZ_MASK,strainI,strainII
#endif
      use mod_common_ice, only : MIZ_MASK
      implicit none
      integer, intent(in) :: ksub  ! subcycling step
      integer :: i, j

      real :: 
     &  divune, divunw, divuse, divusw             ! divergence
     &, tensionne, tensionnw, tensionse, tensionsw ! tension
     &, shearne, shearnw, shearse, shearsw         ! shearing
     &, msq_tension,msq_shear,divu_mean            ! mean square tension, shear; mean of divu
     &, Deltane, Deltanw, Deltase, Deltasw         ! Delta 
     &, c0ne, c0nw, c0se, c0sw                     ! useful combinations
     &, c1ne, c1nw, c1se, c1sw
     &, detne, detnw, detse, detsw
     &, K2ne, K2nw, K2se, K2sw
     &, K3ne, K3nw, K3se, K3sw
     &, vdne, vdnw, vdse, vdsw
     &, etane, etanw, etase, etasw
     &, zetane, zetanw, zetase, zetasw
     &, prssne, prssnw, prssse, prsssw
     %, Eyoung, fac1, fac2, fac, pi, gamp, sep, hfloe
     &, maxstressm, maxstressp,maxstress12
     &, minstressm, minstressp,minstress12
     &, minvpD, maxvpD, cice, minvp

      integer ::
     &  ij      ! loop index, combination of i and j loops
      logical, parameter :: diag=.true.

c --- From the collisional rheology model of Shen et al. (1987)
      real, parameter :: e=0.1           ! restitution coefficient
      real, parameter :: diam=30.        ! diameter of floes
      real, parameter :: C0=1.1         ! max compactness of hexagonal packing
      real, parameter :: K1=(1-e**2)*0.5 ! factor fo ert. velocity assumption
      real, parameter :: rhoi=930.       ! density of sea ice
      real, parameter :: cs=0.89         ! cutoff concentration to turn on MIZ
      real, parameter :: hmin=1.0        ! minimum ice thickness for VP at c=1
      real, parameter :: hmax=2.0        ! maximum MIZ thickness at any concentration

c --- MIZ viscosity ++
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &   mizzeta
     &,  mizeta
     &,  mizprss
     &,  mizvd

      logical  :: USE_EVP!!criteria for EVP/MIZ stresses

      pi = 4.d0*atan(1.d0)
      maxstressm  = -1e14
      maxstressp  = -1e14
      maxstress12 = -1e14
      maxvpD      = -1e14
      minstressm  =  1e14
      minstressp  =  1e14
      minstress12 =  1e14
      minvpD      =  1e14

!$OMP PARALLEL DO PRIVATE(j,i,
!$OMP&     divune,    divunw,   divusw,   divuse,
!$OMP&  tensionne, tensionnw,tensionsw,tensionse,
!$OMP&    shearne,   shearnw,  shearsw,  shearse,
!$OMP&    Deltane,   Deltanw,  Deltasw,  Deltase,
!$OMP&       c0ne,      c0nw,     c0sw,     c0se,
!$OMP&       c1ne,      c1nw,     c1sw,     c1se,
!$OMP&      detne,     detnw,    detsw,    detse, 
!$OMP&       K2ne,      K2nw,     K2sw,     K2se, 
!$OMP&       K3ne,      K3nw,     K3sw,     K3se, 
!$OMP&  vdne, vdnw, vdse, vdsw,
!$OMP&  etane, etanw, etase, etasw,
!$OMP&  zetane, zetanw, zetase, zetasw,
!$OMP&  prssne, prssnw, prssse, prsssw,
!$OMP&       fac , dfloe, Eyoung, fac1, fac2,cice)
!$OMP&         SCHEDULE(STATIC,jblk)

      MIZ_MASK    = 0.0!!test variable stored in mod_common_ice.F
      msq_tension = 0.0
      msq_shear   = 0.0
      divu_mean   = 0.0

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

      ! shearing strain rate  =  2*e_12
      shearne = -cym(i,j)*vvel(i  ,j  ) - dyt(i,j)*vvel(i-1,j  )
     &        -  cxm(i,j)*uvel(i  ,j  ) - dxt(i,j)*uvel(i  ,j-1)
      shearnw = -cyp(i,j)*vvel(i-1,j  ) + dyt(i,j)*vvel(i  ,j  )
     &        -  cxm(i,j)*uvel(i-1,j  ) - dxt(i,j)*uvel(i-1,j-1)
      shearsw = -cyp(i,j)*vvel(i-1,j-1) + dyt(i,j)*vvel(i  ,j-1)
     &        -  cxp(i,j)*uvel(i-1,j-1) + dxt(i,j)*uvel(i-1,j  )
      shearse = -cym(i,j)*vvel(i  ,j-1) - dyt(i,j)*vvel(i-1,j-1)
     &        -  cxp(i,j)*uvel(i  ,j-1) + dxt(i,j)*uvel(i  ,j  )

      ! Define and store strain rate invariants for diagnostics
#if defined(ICE_DYN_DIAG)
      !!TW change: 20140307
      !!averaging of corners was done weirdly
      if (ksub.eq.ndte) then
!!    strainI(i,j)  = 0.5*0.25*tarear(i,j)*(divune + divunw 
!!   &              + divuse + divusw)
!!    strainII(i,j) = 0.5*0.25*tarear(i,j)*sqrt(
!!   &              (tensionne + tensionnw + tensionse + tensionsw)**2
!!   &            + (  shearne +   shearnw +   shearse +   shearsw)**2)
         divu_mean   = .25*(divune + divunw 
     &                    + divuse + divusw)
         msq_tension = .25*(tensionne**2+tensionnw**2
     &                      +tensionse**2+tensionsw**2);
         msq_shear   = .25*(shearne**2+shearnw**2
     &                      +shearse**2+shearsw**2);

         !!principal strains are strainI +/- strainII
         !!NB tension=e11-e22, shear=2*e12
         strainI(i,j)  = 0.5*tarear(i,j)*divu_mean
         strainII(i,j) = 0.5*tarear(i,j)*sqrt(msq_tension+msq_shear)
      endif
#endif

      ! Ice concentration
      cice  = max(aice(i,j),0.01)      ! always > 0.01 
      hfloe = max(vice(i,j)/cice,0.10) ! always > 0.10

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!Define when NOT to use MIZ rheology
      ! VERSION 0 : Unstable
      ! USE_EVP  = (aice(i,j).gt.cs)
      ! VERSION 1 : Stable but not realistic enough
      ! USE_EVP  = (cice.gt.cs .or. hfloe.gt.hmax)
      ! VERSION 2 : Based on a waves-in-ice model. To be tested.
      USE_EVP  = (hfloe .gt. (-1*hmin+((1-cice)*(hmax-hmin)/(1-cs))))
      USE_EVP  = (USE_EVP.or.(hfloe .gt. hmax))
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!Normal EVP stress
      if (USE_EVP) then

      ! Viscous-plastic regime computed with the EVP scheme

      ! Delta (in the denominator of zeta, eta)
      Deltane = sqrt(divune**2 + ecci*(tensionne**2 + shearne**2))
      Deltanw = sqrt(divunw**2 + ecci*(tensionnw**2 + shearnw**2))
      Deltase = sqrt(divuse**2 + ecci*(tensionse**2 + shearse**2))
      Deltasw = sqrt(divusw**2 + ecci*(tensionsw**2 + shearsw**2))

      if (evp_damping) then
        ! enforce damping criterion
        c0ne=min(prss(i,j)/max(Deltane,4.*tinyarea(i,j)),rcon_evp(i,j))
        c0nw=min(prss(i,j)/max(Deltanw,4.*tinyarea(i,j)),rcon_evp(i,j))
        c0se=min(prss(i,j)/max(Deltase,4.*tinyarea(i,j)),rcon_evp(i,j))
        c0sw=min(prss(i,j)/max(Deltasw,4.*tinyarea(i,j)),rcon_evp(i,j))
        prs_sig(i,j) = prss(i,j)*Deltane/max(Deltane,4.*tinyarea(i,j)) ! ne
      else
        ! original version
        c0ne = prss(i,j)/max(Deltane,tinyarea(i,j))
        c0nw = prss(i,j)/max(Deltanw,tinyarea(i,j))
        c0se = prss(i,j)/max(Deltase,tinyarea(i,j))
        c0sw = prss(i,j)/max(Deltasw,tinyarea(i,j))
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

      ! modified EVP as in Bouillon et al. OM 2013. 
      stressm_1(i,j) = (stressm_1(i,j) + c1ne*ecci*tensionne)*denom
      stressm_2(i,j) = (stressm_2(i,j) + c1nw*ecci*tensionnw)*denom
      stressm_3(i,j) = (stressm_3(i,j) + c1sw*ecci*tensionsw)*denom
      stressm_4(i,j) = (stressm_4(i,j) + c1se*ecci*tensionse)*denom

      stress12_1(i,j) = (stress12_1(i,j) + c1ne*ecci*shearne*.5)*denom
      stress12_2(i,j) = (stress12_2(i,j) + c1nw*ecci*shearnw*.5)*denom
      stress12_3(i,j) = (stress12_3(i,j) + c1sw*ecci*shearsw*.5)*denom
      stress12_4(i,j) = (stress12_4(i,j) + c1se*ecci*shearse*.5)*denom

      else
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Collisional regime computed following Shen et al. (1987)
      MIZ_MASK(i,j)  = 1.0!!test variable stored in mod_common_ice.F

      ! Convert strain rates units from m^2/s to 1/s
      ! (breaks curvilinear factors, ok for uniform grids)
      divune = divune*uarear(i,j)
      divunw = divunw*uarear(i,j)
      divusw = divusw*uarear(i,j)
      divuse = divuse*uarear(i,j)

      tensionne = tensionne*uarear(i,j)
      tensionnw = tensionnw*uarear(i,j)
      tensionsw = tensionsw*uarear(i,j)
      tensionse = tensionse*uarear(i,j)

      shearne = shearne*uarear(i,j)
      shearnw = shearnw*uarear(i,j)
      shearsw = shearsw*uarear(i,j)
      shearse = shearse*uarear(i,j)

      ! Used in the collisional rheology ! 1/s^2
      !!det=e_11*e_22-e_12^2, so these should be minus signs
      !![TW 20140317]
      !!detne = 0.25*(divune**2+tensionne**2+shearne**2)
      !!detnw = 0.25*(divunw**2+tensionnw**2+shearnw**2)
      !!detne = 0.25*(divusw**2+tensionsw**2+shearsw**2)!NB should also be sw on LHS
      !!detne = 0.25*(divuse**2+tensionse**2+shearse**2)
      detne = 0.25*(divune**2-tensionne**2-shearne**2)
      detnw = 0.25*(divunw**2-tensionnw**2-shearnw**2)
      detsw = 0.25*(divusw**2-tensionsw**2-shearsw**2)
      detne = 0.25*(divuse**2-tensionse**2-shearse**2)

      ! K1, K2 and K3 are needed when you run the model without advecting
      ! the perturbation velocity. To set up perturbation velocity, 
      ! you need a source (e.g. a wave parameterization).

      ! K2                               ! 1/s
      fac  = 4*sqrt(2.d0)*(1+e)/(pi**2) - sqrt(2.d0)*(1-e**2)/pi
      K2ne = fac*divune
      K2nw = fac*divunw
      K2sw = fac*divusw
      K2se = fac*divuse

      ! K3                               ! 1/s^2
      fac  = -4*(1+e)/(3*pi)
      K3ne = K1*(0.375*divune**2 - 0.5*detne) + fac*(divune**2 - detne)
      K3nw = K1*(0.375*divunw**2 - 0.5*detnw) + fac*(divunw**2 - detnw)
      K3sw = K1*(0.375*divusw**2 - 0.5*detsw) + fac*(divusw**2 - detsw)
      K3se = K1*(0.375*divuse**2 - 0.5*detse) + fac*(divuse**2 - detse)

      ! Diagnostic v'/D based on assumption of stationary field ! 1/s
      vdne = 0.5*(sqrt(K2ne**2 - 4*K1*K3ne) - K2ne)/K1
      vdnw = 0.5*(sqrt(K2nw**2 - 4*K1*K3nw) - K2nw)/K1
      vdsw = 0.5*(sqrt(K2sw**2 - 4*K1*K3sw) - K2sw)/K1
      vdse = 0.5*(sqrt(K2se**2 - 4*K1*K3se) - K2se)/K1

      minvpD=min(minvpD,vdne,vdnw,vdsw,vdse)
      maxvpD=max(maxvpD,vdne,vdnw,vdsw,vdse)

      ! Background values of vd
      minvp=1.e-8
      vdne=max(vdne,minvp/diam)
      vdnw=max(vdnw,minvp/diam)
      vdse=max(vdse,minvp/diam)
      vdsw=max(vdsw,minvp/diam)

      ! Final viscosities ---
      ! TODO: aice, vice to triangles?
      ! Thickness of an ice floe
      cice  = max(aice(i,j),0.01)      ! always > 0.01 
      hfloe = max(vice(i,j)/cice,0.10) ! always > 0.1

      ! Simplifying factor
      sep=cice**3.5/max(0.01,sqrt(C0)-sqrt(cice))     ! unitless
      gamp=rhoi*(diam**2)*hfloe*sqrt(2.0)*0.25*(1+e)*sep/(pi**2)  ! kg

      ! Shear viscosity
      ! etane = pi*gamp*vdne/(3*sqrt(2.0))   ! kg/s
      ! etanw = pi*gamp*vdnw/(3*sqrt(2.0))
      ! etasw = pi*gamp*vdsw/(3*sqrt(2.0))
      ! etase = pi*gamp*vdse/(3*sqrt(2.0))

      ! Bulk viscosity
      ! zetane = pi*gamp*vdne/sqrt(2.0)      ! kg/s
      ! zetanw = pi*gamp*vdnw/sqrt(2.0)
      ! zetasw = pi*gamp*vdsw/sqrt(2.0)
      ! zetase = pi*gamp*vdse/sqrt(2.0)

      ! if (j0+j==100) print *,'prssne',prssne
      ! if (j0+j==100) print *,'etane' ,etane
      
      if (evp_damping) then
        ! DD -- enforce damping criterion
        c0ne = min(gamp*vdne,rcon_miz(i,j))
        c0nw = min(gamp*vdnw,rcon_miz(i,j))
        c0sw = min(gamp*vdsw,rcon_miz(i,j))
        c0se = min(gamp*vdse,rcon_miz(i,j))
      else
        ! DD -- regular damping            ! kg/s
        c0ne = gamp*vdne
        c0nw = gamp*vdnw
        c0sw = gamp*vdsw
        c0se = gamp*vdse
      endif

      c1ne = c0ne*dte2T
      c1nw = c0nw*dte2T
      c1sw = c0sw*dte2T
      c1se = c0se*dte2T

      ! denom1 = 1.0/(1.0 + dte2T) is the same as in the evp model
      denom2 = 1.0/(1.0 + 3*dte2T)

      !-----------------------------------------------------------------
      ! the stresses                            ! kg/s^2
      ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
      !-----------------------------------------------------------------

      stressp_1(i,j) = (stressp_1(i,j) + 
     &                  c1ne*(pi*sqrt(2.0)*divune-vdne))*denom
      stressp_2(i,j) = (stressp_2(i,j) +
     &                  c1nw*(pi*sqrt(2.0)*divunw-vdnw))*denom
      stressp_3(i,j) = (stressp_3(i,j) +
     &                  c1sw*(pi*sqrt(2.0)*divusw-vdsw))*denom
      stressp_4(i,j) = (stressp_4(i,j) +
     &                  c1se*(pi*sqrt(2.0)*divuse-vdse))*denom

      stressm_1(i,j) = (stressm_1(i,j) +
     &                  c1ne*pi*sqrt(2.0)*tensionne)*denom2
      stressm_2(i,j) = (stressm_2(i,j) +
     &                  c1nw*pi*sqrt(2.0)*tensionnw)*denom2
      stressm_3(i,j) = (stressm_3(i,j) +
     &                  c1sw*pi*sqrt(2.0)*tensionsw)*denom2
      stressm_4(i,j) = (stressm_4(i,j) +
     &                  c1se*pi*sqrt(2.0)*tensionse)*denom2

      stress12_1(i,j) = (stress12_1(i,j) +
     &                   c1ne*pi*sqrt(2.0)*shearne*.5)*denom2
      stress12_2(i,j) = (stress12_2(i,j) +
     &                   c1nw*pi*sqrt(2.0)*shearnw*.5)*denom2
      stress12_3(i,j) = (stress12_3(i,j) +
     &                   c1sw*pi*sqrt(2.0)*shearsw*.5)*denom2
      stress12_4(i,j) = (stress12_4(i,j) +
     &                   c1se*pi*sqrt(2.0)*shearse*.5)*denom2


      maxstressm=max(stressm_1(i,j),stressm_2(i,j),stressm_3(i,j),
     &               stressm_4(i,j),maxstressm)
      maxstressp=max(stressp_1(i,j),stressp_2(i,j),stressp_3(i,j),
     &               stressp_4(i,j),maxstressp)
      maxstress12=max(stress12_1(i,j),stress12_2(i,j),
     &                stress12_3(i,j),stress12_4(i,j),maxstress12)

      minstressm=min(stressm_1(i,j),stressm_2(i,j),stressm_3(i,j),
     &               stressm_4(i,j),minstressm)
      minstressp=min(stressp_1(i,j),stressp_2(i,j),stressp_3(i,j),
     &               stressp_4(i,j),minstressp)
      minstress12=min(stress12_1(i,j),stress12_2(i,j),
     &                stress12_3(i,j),stress12_4(i,j),minstress12)

      ! if (j0+j==100) print *,'stressp',stressp_1(i,j)
      ! if (j0+j==100) print *,'stressm',stressm_1(i,j)
      ! if (j0+j==100) print *,'stress12',stress12_1(i,j)
      !     if (j+j0==101.and.i+i0==253.and.ksub==2) print *,
      ! &      'stressm(i=253,j=101,:):',stressm_1(i,j)

      ! Keep viscosities for study
      mizzeta(i,j) = zetane
      mizeta(i,j)  = etane
      mizprss(i,j) = gamp*vdne**2
      mizvd(i,j)   = vdne
      
      endif
      !!end of EVP/MIZ check
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      endif!!check for ice presence
      enddo!!i
      enddo!!j
!$OMP END PARALLEL DO 

!      if (diag) then
!      call xcmaxr(maxstressm)
!      call xcminr(minstressm)
!      if (mnproc==1) print *,'min/max stressm:',minstressm,maxstressm
!      call xcmaxr(maxstressp)
!      call xcminr(minstressp)
!      if (mnproc==1) print *,'min/max stressp:',minstressp,maxstressp
!      call xcmaxr(maxstress12)
!      call xcminr(minstress12)
!      if (mnproc==1) print *,'min/max stress12:',minstress12,maxstress12
!      call xcmaxr(maxvpD)
!      call xcminr(minvpD)
!      if (mnproc==1) print *,'min/max vpD:',minvpD,maxvpD
!      end if
!      !call xcstop('(evp_stress_miz)')


      end subroutine evp_stress_miz
!!******


!=======================================================================
!BOP
!
! !IROUTINE: stepu - surface stresses and integrates mom eqn for u,v
!
! !INTERFACE:
!
      subroutine evp_stepu
!
! !DESCRIPTION:
!
! Calculation of the surface stresses \\
! Integration of the momentum equation to find velocity (u,v)
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke
!         Fluid Dynamics Group, Los Alamos National Laboratory
!
! !USES:
!
      !KAL use ice_flux
      use mod_xc, only : mnproc
      use mod_evp
#if defined(ICE_DYN_DIAG)
      use mod_common_ice, only : stressp, stressm, stress12
#endif
      implicit none
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer :: i, j

      real :: 
     &   vrel,cca,ccb,ab2,cc1,cc2, taux,tauy
     &,  str(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,8)
     &,  worka(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ! KAL from ice_work 
     &,  workb(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ! KAL from ice_work 
     &,  ssigpn, ssigps, ssigpe, ssigpw
     &,  ssigmn, ssigms, ssigme, ssigmw
     &,  ssig12n, ssig12s, ssig12e, ssig12w
     &,  ssigp1, ssigp2, ssigm1, ssigm2, ssig121, ssig122
     &,  csigpne, csigpnw, csigpse, csigpsw
     &,  csigmne, csigmnw, csigmse, csigmsw
     &,  csig12ne, csig12nw, csig12se, csig12sw
     &,  str12ew, str12we, str12ns, str12sn
     &,  strp_tmp, strm_tmp, str12_tmp

      integer ::
     &   ij      ! loop index, combination of i and j loops

      real, parameter :: p111 = 1.0/9.0
      real, parameter :: p166 = 1.0/6.0
      real, parameter :: p222 = 2.0/9.0
      real, parameter :: p333 = 1.0/3.0
      real, parameter :: p055 = p111*.5
      real, parameter :: p027 = p055*0.5


Ctst  real :: frac, fcortmp, masstmp 
Ctst  integer :: tmpproc

      !-----------------------------------------------------------------
      ! combinations of the stresses for the momentum equation ! kg/s^2
      !-----------------------------------------------------------------

      str(:,:,:) = 0.0

Cold  !KAL - TODO fix for margin
Cold  do ij=1,icellt
Cold    i = indxti(ij)
Cold    j = indxtj(ij)
Cold  if (i.ge.1-imargin .and. i.le.ii+imargin .and.
Cold &    j.ge.1-imargin .and. j.le.jj+imargin ) then

!$OMP PARALLEL DO PRIVATE(j,i,
!$OMP&     ssigpn ,    ssigps ,   ssigpe ,   ssigpw ,
!$OMP&     ssigp1 ,    ssigp2 , 
!$OMP&     ssigmn ,    ssigms ,   ssigme ,   ssigmw ,
!$OMP&     ssigm1 ,    ssigm2 , 
!$OMP&     ssig12n,    ssig12s,   ssig12e,   ssig12w,
!$OMP&     ssig121,    ssig122, 
!$OMP&     Csigpne , Csigpnw , Csigpsw , Csigpse ,
!$OMP&     Csigmne , Csigmnw , Csigmsw , Csigmse ,
!$OMP&     Csig12ne, Csig12nw, Csig12sw, Csig12se,
!$OMP&      str12ew,  str12we,  str12ns,  str12sn,
!$OMP&      strp_tmp,  strm_tmp)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-imargin,jj+imargin
      do i=1-imargin,ii+imargin
      if (icetmask(i,j)) then

#if defined(ICE_DYN_DIAG)
        stressp(i,j) = (stressp_1(i,j) + stressp_2(i,j)
     &               + stressp_3(i,j) + stressp_4(i,j))/4.0
        stressm(i,j) = (stressm_1(i,j) + stressm_2(i,j)
     &               + stressm_3(i,j) + stressm_4(i,j))/4.0
        stress12(i,j) = (stress12_1(i,j) + stress12_2(i,j)
     &                + stress12_3(i,j) + stress12_4(i,j))/4.0
#endif

        ssigpn  = stressp_1(i,j) + stressp_2(i,j)
        ssigps  = stressp_3(i,j) + stressp_4(i,j)
        ssigpe  = stressp_1(i,j) + stressp_4(i,j)
        ssigpw  = stressp_2(i,j) + stressp_3(i,j)
        ssigp1  =(stressp_1(i,j) + stressp_3(i,j))*p055
        ssigp2  =(stressp_2(i,j) + stressp_4(i,j))*p055

        ssigmn  = stressm_1(i,j) + stressm_2(i,j)
        ssigms  = stressm_3(i,j) + stressm_4(i,j)
        ssigme  = stressm_1(i,j) + stressm_4(i,j)
        ssigmw  = stressm_2(i,j) + stressm_3(i,j)
        ssigm1  =(stressm_1(i,j) + stressm_3(i,j))*p055
        ssigm2  =(stressm_2(i,j) + stressm_4(i,j))*p055

        ssig12n = stress12_1(i,j) + stress12_2(i,j)
        ssig12s = stress12_3(i,j) + stress12_4(i,j)
        ssig12e = stress12_1(i,j) + stress12_4(i,j)
        ssig12w = stress12_2(i,j) + stress12_3(i,j)
        ssig121 =(stress12_1(i,j) + stress12_3(i,j))*p111
        ssig122 =(stress12_2(i,j) + stress12_4(i,j))*p111

        csigpne = p111*stressp_1(i,j) + ssigp2 + p027*stressp_3(i,j)
        csigpnw = p111*stressp_2(i,j) + ssigp1 + p027*stressp_4(i,j)
        csigpsw = p111*stressp_3(i,j) + ssigp2 + p027*stressp_1(i,j)
        csigpse = p111*stressp_4(i,j) + ssigp1 + p027*stressp_2(i,j)

        csigmne = p111*stressm_1(i,j) + ssigm2 + p027*stressm_3(i,j)
        csigmnw = p111*stressm_2(i,j) + ssigm1 + p027*stressm_4(i,j)
        csigmsw = p111*stressm_3(i,j) + ssigm2 + p027*stressm_1(i,j)
        csigmse = p111*stressm_4(i,j) + ssigm1 + p027*stressm_2(i,j)

        csig12ne = p222*stress12_1(i,j) + ssig122 + p055*stress12_3(i,j)
        csig12nw = p222*stress12_2(i,j) + ssig121 + p055*stress12_4(i,j)
        csig12sw = p222*stress12_3(i,j) + ssig122 + p055*stress12_1(i,j)
        csig12se = p222*stress12_4(i,j) + ssig121 + p055*stress12_2(i,j)

        str12ew = dxt2(i,j)*(p333*ssig12e + p166*ssig12w)
        str12we = dxt2(i,j)*(p333*ssig12w + p166*ssig12e)
        str12ns = dyt2(i,j)*(p333*ssig12n + p166*ssig12s)
        str12sn = dyt2(i,j)*(p333*ssig12s + p166*ssig12n)

      !-----------------------------------------------------------------
      ! for dF/dx (u momentum)                  ! N = kg m/s^2
      !-----------------------------------------------------------------
        strp_tmp  = dyt4(i,j)*(p333*ssigpn  + p166*ssigps)
        strm_tmp  = dyt4(i,j)*(p333*ssigmn  + p166*ssigms)

        ! northeast (i,j)
        str(i,j,1) = -strp_tmp - strm_tmp - str12ew 
     &  + dxhy(i,j)*(-csigpne + csigmne) + dyhx(i,j)*csig12ne

        ! northwest (i+1,j)
        str(i,j,2) = strp_tmp + strm_tmp - str12we 
     &  + dxhy(i,j)*(-csigpnw + csigmnw) + dyhx(i,j)*csig12nw

        strp_tmp  = dyt4(i,j)*(p333*ssigps  + p166*ssigpn)
        strm_tmp  = dyt4(i,j)*(p333*ssigms  + p166*ssigmn)

        ! southeast (i,j+1)
        str(i,j,3) = -strp_tmp - strm_tmp + str12ew 
     &  + dxhy(i,j)*(-csigpse + csigmse) + dyhx(i,j)*csig12se

        ! southwest (i+1,j+1)
        str(i,j,4) = strp_tmp + strm_tmp + str12we 
     &  + dxhy(i,j)*(-csigpsw + csigmsw) + dyhx(i,j)*csig12sw


      !-----------------------------------------------------------------
      ! for dF/dy (v momentum)
      !-----------------------------------------------------------------
        strp_tmp  = dxt4(i,j)*(p333*ssigpe  + p166*ssigpw)
        strm_tmp  = dxt4(i,j)*(p333*ssigme  + p166*ssigmw)

        ! northeast (i,j)
        str(i,j,5) = -strp_tmp + strm_tmp - str12ns 
     &  - dyhx(i,j)*(csigpne + csigmne) + dxhy(i,j)*csig12ne

        ! southeast (i,j+1)
        str(i,j,6) = strp_tmp - strm_tmp - str12sn 
     &  - dyhx(i,j)*(csigpse + csigmse) + dxhy(i,j)*csig12se

        strp_tmp  = dxt4(i,j)*(p333*ssigpw  + p166*ssigpe)
        strm_tmp  = dxt4(i,j)*(p333*ssigmw  + p166*ssigme)

        ! northwest (i+1,j)
        str(i,j,7) = -strp_tmp + strm_tmp + str12ns
     &  - dyhx(i,j)*(csigpnw + csigmnw) + dxhy(i,j)*csig12nw

        ! southwest (i+1,j+1)
        str(i,j,8) = strp_tmp - strm_tmp + str12sn 
     &  - dyhx(i,j)*(csigpsw + csigmsw) + dxhy(i,j)*csig12sw

Cold    end if
Cold
Cold  enddo                     ! ij

      endif
      enddo
      enddo
!$OMP END PARALLEL DO 

      ! - KAL TODO ?
      !KAL call bound_narr_ne(8,str)
      !call xctilr(str ( 1-nbdy,1-nbdy,  1),1,   8, 6,6, halo_ps)

      !Correct imargin in stead of re-tiling
      !imargin=imargin -1

      ! KAL - everything is calculated. Also for "wrong" points


      !-----------------------------------------------------------------
      ! set velocity and stress to zero on land and (nearly) open water
      ! use working arrays to avoid u and v vector dependency
      !_________________________________________________________________

!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-imargin,jj+imargin
      do i=1-imargin,ii+imargin
         if (.not.iceumask(i,j)) then
            worka(i,j) = 0.0
            workb(i,j) = 0.0
            strocnx(i,j) = 0.0
            strocny(i,j) = 0.0
            strintx(i,j) = 0.0
            strinty(i,j) = 0.0
         endif
      enddo
      enddo
!$OMP END PARALLEL DO 

      !-----------------------------------------------------------------
      ! integrate the momentum equation
      !-----------------------------------------------------------------

      !KAL - TODO fix for margin
Cold  do ij=1,icellu
Cold    i=indxui(ij)
Cold    j=indxuj(ij)
Cold  if (i.ge.1-imargin .and. i.le.ii+imargin .and.
Cold &    j.ge.1-imargin .and. j.le.jj+imargin ) then
!$OMP PARALLEL DO PRIVATE(j,i,vrel,taux,tauy,cca,ccb,ab2,cc1,cc2)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
      if (iceumask(i,j)) then

        ! (magnitude of relative ocean current)*rhow*drag*aice
        vrel = aiu(i,j)*dragw*sqrt((uocn(i,j) - uvel(i,j))**2 + 
     &                             (vocn(i,j) - vvel(i,j))**2)  ! m/s
        ! ice/ocean stress
        taux = vrel*waterx(i,j) ! NOTE this is not the entire
        tauy = vrel*watery(i,j) ! ocn stress term

        ! alpha, beta are defined in Hunke and Dukowicz (1997), section 3.2
        cca = umassdtei(i,j) + vrel * cosw         ! alpha, kg/m^2 s
        ccb = fm(i,j)        + vrel * sinw         ! beta,  kg/m^2 s
        ab2 = cca**2 + ccb**2

        ! divergence of the internal stress tensor
        strintx(i,j) = uarear(i,j)*
     &      (str(i,j,1) + str(i+1,j,2) + str(i,j+1,3) + str(i+1,j+1,4))
        strinty(i,j) = uarear(i,j)*
     &      (str(i,j,5) + str(i,j+1,6) + str(i+1,j,7) + str(i+1,j+1,8))

        ! finally, the velocity components
        cc1 = strintx(i,j) + forcex(i,j) + taux 
     &       + umassdtei(i,j)*uvel(i,j)
        cc2 = strinty(i,j) + forcey(i,j) + tauy 
     &       + umassdtei(i,j)*vvel(i,j)

Ctst    !KAL test
Ctst    frac=aiu(i,j)
Ctst    fcortmp=fm(i,j)
Ctst    masstmp=umassdtei(i,j)
Ctst    tmpproc=mnproc

        worka(i,j) = (cca*cc1 + ccb*cc2)/ab2              ! m/s
        workb(i,j) = (cca*cc2 - ccb*cc1)/ab2


      !-----------------------------------------------------------------
      ! ocean-ice stress for coupling
      ! scale to full grid cell
      !-----------------------------------------------------------------
        strocnx(i,j) = taux/aiu(i,j)  
        strocny(i,j) = tauy/aiu(i,j)  

      endif
      enddo
      enddo
!$OMP END PARALLEL DO 
Cold  endif
Cold  enddo                     ! ij

      !-----------------------------------------------------------------
      ! recover u and v arrays
      !-----------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-imargin,jj+imargin
      do i=1-imargin,ii+imargin
         uvel(i,j) = worka(i,j)
         vvel(i,j) = workb(i,j)
      enddo
      enddo
!$OMP END PARALLEL DO 

      end subroutine evp_stepu

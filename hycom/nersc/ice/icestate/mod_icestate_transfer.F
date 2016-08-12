      module mod_icestate_transfer
      contains

      subroutine hycom2istate(m,n)
      use mod_icestate , only : icestate
      use mod_icestate_fluxes 
      use mod_common_ice
      use mod_forcing_nersc
      implicit none
      integer, intent(in) :: m,n
      integer :: i,j, patch,l
      real    :: tmlfac,smlfac
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: work2d
      include 'common_blocks.h'

c --- slp -- mbar
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         Islp  (i,j)=slp(i,j,l0)*w0+slp(i,j,l1)*w1
     &              +slp(i,j,l2)*w2+slp(i,j,l3)*w3
      end do
      end do
!$OMP END PARALLEL DO
c
c --- radfl= net       radiative thermal flux (w/m**2) positive into ocean
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         Iradflx(i,j)=radflx(i,j,l0)*w0+radflx(i,j,l1)*w1
     &               +radflx(i,j,l2)*w2+radflx(i,j,l3)*w3
      end do
      end do
!$OMP END PARALLEL DO
c
c --- swfl = shortwave radiative thermal flux (w/m**2) positive into ocean
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         Iswflx(i,j)=swflx (i,j,l0)*w0+swflx (i,j,l1)*w1
     &              +swflx (i,j,l2)*w2+swflx (i,j,l3)*w3
      end do
      end do
c
!$OMP END PARALLEL DO
c --- airt = air temperature (C)
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         Iairtmp(i,j)=airtmp(i,j,l0)*w0+airtmp(i,j,l1)*w1
     &               +airtmp(i,j,l2)*w2+airtmp(i,j,l3)*w3
         Iairtmp(i,j)=Iairtmp(i,j)+t0deg ! Must be Kelvin in icestate
      end do
      end do
!$OMP END PARALLEL DO
c
c ---   vpmx = water vapor mixing ratio (kg/kg)
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         Ivapmix(i,j)=vapmix(i,j,l0)*w0+vapmix(i,j,l1)*w1
     &               +vapmix(i,j,l2)*w2+vapmix(i,j,l3)*w3
      end do
      end do
!$OMP END PARALLEL DO
c
c --- prcp = precipitation (m/sec; positive into ocean)
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         Iprecip(i,j)=precip(i,j,l0)*w0+precip(i,j,l1)*w1
     &               +precip(i,j,l2)*w2+precip(i,j,l3)*w3
      end do
      end do
!$OMP END PARALLEL DO
c
c --- relhum = relhum ()
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         Irelhum(i,j)=relhum(i,j,l0)*w0+relhum(i,j,l1)*w1
     &               +relhum(i,j,l2)*w2+relhum(i,j,l3)*w3
      end do
      end do
!$OMP END PARALLEL DO
c
c --- wind = windspeed ()
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         Iwndspd(i,j)=wndspd(i,j,l0)*w0+wndspd(i,j,l1)*w1
     &               +wndspd(i,j,l2)*w2+wndspd(i,j,l3)*w3
      end do
      end do
!$OMP END PARALLEL DO
c
c --- clouds 
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         Iclouds(i,j)=clouds(i,j,l0)*w0+clouds(i,j,l1)*w1
     &               +clouds(i,j,l2)*w2+clouds(i,j,l3)*w3
      end do
      end do
!$OMP END PARALLEL DO
c --- taux 
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         Itaux(i,j)=taux(i,j,l0)*w0+taux(i,j,l1)*w1
     &             +taux(i,j,l2)*w2+taux(i,j,l3)*w3
      end do
      end do
!$OMP END PARALLEL DO
c --- tauy 
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         Itauy(i,j)=tauy(i,j,l0)*w0+tauy(i,j,l1)*w1
     &              +tauy(i,j,l2)*w2+tauy(i,j,l3)*w3
      end do
      end do
!$OMP END PARALLEL DO
c --- tauxice 
      Itauxice=tauxice
c --- tauyice 
      Itauyice=tauyice

C$OMP PARALLEL DO PRIVATE(i,j,l,patch) 
C$OMP&SCHEDULE(STATIC,1)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         !icestate(i,j)%hml = dp  (i,j,1,n)  / onem
         !icestate(i,j)%sml = saln(i,j,1,n)
         !icestate(i,j)%tml = temp(i,j,1,n) + t0deg
         ! Avoids decoupling issues
         icestate(i,j)%hml = .5*(dp  (i,j,1,n) + dp  (i,j,1,m))  / onem
         icestate(i,j)%sml = .5*(saln(i,j,1,n) + saln(i,j,1,m))
         icestate(i,j)%tml = .5*(temp(i,j,1,n) + temp(i,j,1,m)) + t0deg
      enddo
      enddo
COMP END PARALLEL DO
c
c --- This may be unnecessary, but just in case
      call xctilr(Islp     ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(Iradflx  ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(Iswflx   ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(Iairtmp  ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(Ivapmix  ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(Ivapmix  ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(Iprecip  ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(Irelhum  ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(Iwndspd  ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(Iclouds  ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(Itaux    ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(Itauy    ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(Itauxice ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(Itauyice ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      end subroutine hycom2istate



c --- Transfer variables from ICESTATE to HYCOM. Ocean Surface fluxes
      subroutine istate2hycom
      use mod_xc
      use mod_icestate        , ONLY :  icestate
      use mod_icestate_fluxes
      use mod_forcing_nersc
      implicit none
c     TODO: cekman not  needed anymore ?
      !real, parameter :: cekman=0.7
      !real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: tst
      real, dimension(itdm,jtdm) :: tst
      integer i,j,l
      include 'common_blocks.h'


c$OMP PARALLEL DO PRIVATE(j,i) 
c$OMP&SCHEDULE (STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         salflx(i,j) = 0.
         sswflx(i,j) = 0.
         surflx(i,j) = 0.
         ustar (i,j) = 0.
         hekman(i,j) = 0.
      end do
      end do
c$OMP END PARALLEL DO

c --- Set fluxes used by HYCOM
c$OMP PARALLEL DO PRIVATE(j,l,i) 
c$OMP&SCHEDULE (STATIC,jblk)
      do j=1-margin,jj+margin
      do l=1,isp(j)
      do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      
c ---       Salinity flux
            salflx(i,j) = Isalflx(i,j) 
c
c ---       Heat fluxes
            sswflx(i,j) = Isswflx(i,j)
            surflx(i,j) = Isurflx(i,j)
c
c ---       Corio from HYCOM
            ustar (i,j) = Iustar (i,j)
            hekman(i,j) = cekman*ustar(i,j)/ abs(corio(i,j))
c
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

      if (mnproc==1) print *,'KAL:Istate2hycom :', cekman
      if (mnproc==1) print *,1-nbdy,idm+nbdy,1-nbdy,jdm+nbdy,itdm

c --- Tile HYCOM fluxes
      call xctilr(salflx   ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(sswflx   ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(surflx   ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(ustar    ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xctilr(hekman   ( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)


c --- Diagnostics
      call xcaget(tst,sswflx,0)
      if (mnproc==1) print *,'KAL:Istate2hycom sswflx :',
     &   maxval(tst),minval(tst)
      call xcaget(tst,surflx,0)
      if (mnproc==1) print *,'KAL:Istate2hycom surflx :',
     &   maxval(tst),minval(tst)
      call xcaget(tst,salflx,0)
      if (mnproc==1) print *,'KAL:Istate2hycom sswflx :',
     &   maxval(tst),minval(tst)


cdiag  if (mnproc==2) print *,'point ssw',i0+52,j0+82,52,82,
cdiag &   sswflx(52,82),Isswflx(52,82),radfl_day(52,82)



      end subroutine istate2hycom


#if defined (EVP)
c --- Transfer variables from EVP to ICESTATE. Only ice divergence 
      subroutine evp2istate
      use mod_xc
      use mod_icestate_fluxes , ONLY : ISTATE_divu=>divu
      use mod_evp          ,    ONLY : EVP_divu=>divu
      implicit none
      integer :: i,j
c      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: tst
      real, dimension(itdm,jtdm) :: tst

c$OMP PARALLEL DO PRIVATE(j,i) 
c$OMP&SCHEDULE (STATIC,jblk)
      do j=1-margin,jj+margin
      do i=1-margin,ii+margin
         ISTATE_divu(i,j)=EVP_divu(i,j)
      end do
      end do
cOMP END PARALLEL DO
      call xctilr(ISTATE_divu( 1-nbdy,1-nbdy),1,   1, 6,6, halo_ps)
      call xcaget(tst,ISTATE_divu,0)
      if (mnproc==1) print *,'KAL:evp2istate DIVU :',
     &   maxval(tst),minval(tst)
      end subroutine evp2istate
#endif



      end module mod_icestate_transfer

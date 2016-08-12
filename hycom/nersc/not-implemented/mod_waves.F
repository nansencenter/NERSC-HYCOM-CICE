      module mod_waves
c      use mod_xc
      implicit none
      private :: lambda,wc
c      include 'common_blocks.h'

      contains

c      subroutine waves_init
      ! Load the lookup table of attenuation coefficients
c      end subroutine waves_init

      subroutine waves_advect
      ! Routine for estimating waves-in-ice height from an incident wave
      ! field and using the attenuation coefficient of Kohout and Meylan (2008).
      use mod_xc
      use mod_forcing_nersc
      use mod_hycom_nersc
      use mod_advem
      use mod_evp, only : aice,vice
      include 'common_blocks.h'
      
      integer             :: nw=10

c      real, intent(out), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
c     &        dfloe
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)    ::
     &        iceh
     &       ,dummy1,dummy2
     &       ,uwavflx,vwavflx
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,nw) ::
     &        tmp,wamp                     ! Individual wave amplitude [m]

      real, dimension(nw) :: period,wlng,wspd

      real                :: alpha,aa,bb
      real                :: offset,flxdiv

      integer             :: i,j,k,l,n,w
      integer, parameter  :: nwavesteps=200
      real,    parameter  :: twave=400
      real,    parameter  :: ds=3.5e3       ! grid cell spacing (temporary)
      real,    parameter  :: gravity=9.81
      character(len=80)   :: ncfil

      do w = 1,nw
         period(w) = i+5
         wlng(w)   = gravity*period(w)**2/(2*pi)
         wspd(w)   = 0.5*sqrt(gravity*wlng(w)/(2*pi))
      end do

      wamp = 0.0
      tmp  = 0.0

      print*, 'Wave periods'
      print*, period
      print*, 'Wavelenght'
      print*, wlng
      print*, 'Wave speed'
      print*, wspd

      margin=nbdy
C$OMP PARALLEL DO PRIVATE(j,l,i) 
C$OMP&SCHEDULE(STATIC,jblk)
      do j=1-nbdy,jdm+nbdy
         do i=1-nbdy,idm+nbdy
            dummy1(i,j)  = tenm
            dummy2(i,j)  = dummy1(i,j)
            uwavflx(i,j) = 0.
            vwavflx(i,j) = 0.
         end do
      end do
C$OMP END PARALLEL DO
c --- Set up wave variables as "layered" - so that integrating variable 
c --- over a layer retrieves origial fields.
c --- Set up uwavflx and vwavflx
      margin=nbdy-1 ! due to i+1,i-1
C$OMP PARALLEL DO PRIVATE(j,l,i) 
C$OMP&SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
c ---    Set up u-fluxes
         do l=1,isu(j)
         do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            uwavflx(i,j) = uwave(i,j)*scuy(i,j)*dummy1(i-1,j)
         end do
         end do
c ---    Set up v-fluxes
         do l=1,isv(j)
         do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            vwavflx(i,j) = vwave(i,j)*scvx(i,j)*dummy1(i,j-1)
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
         flxdiv=((uwavflx(i+1,j)-uwavflx(i,j))
     &          +(vwavflx(i,j+1)-vwavflx(i,j)))*twave*scp2i(i,j)
         dummy2(i,j)=dummy1(i,j)-flxdiv
         offset=min(0.,dummy1(i,j),dummy2(i,j))
         dummy2(i,j)=dummy2(i,j)-offset
         dummy1(i,j)=dummy1(i,j)-offset
      end do
      end do
      end do
C$OMP END PARALLEL DO

      ! Construction of the wave spectrum
      ! Pierson-Moskowitz spectrum
      aa = 8.1e-3*gravity
      bb = 1.25
     
      do w = 1,nw 
         tmp(:,:,w) = aa*sqrt(period(w)/(2.0*pi))*
     &                   exp(-1.*bb*(period(w)/mwp(:,:))**4)
         wamp(:,:,w) = sqrt(4.0*tmp(:,:,w)*pi/period(w))
      end do

      do n = 1,nwavesteps
         do w = 1,nw
            call advem(1,wamp(:,:,w),wamp(:,:,w),wspd(w),wspd(w),
     &                 dummy1,dummy2,0.,scp2,scp2i,twave)
         end do
      end do

      !call xcaget(h,swh,0)
      !call xcaget(t,mwp,0)
      !call xcaget(lon,plon,0)
      !call xcaget(lat,plat,0)
      !call xcaget(icec,aice,0)
      !call xcaget(icev,vice,0)
 
      ! Advect waves (TEST)

      alpha = 0.002

      ! Dump waves-in-ice variables to verify the propagation scheme.
      !ncfil='wavesinice.nc'
      !call xcaget(icec,aice,0)
      !call ncdraft(trim(ncfil),icec,'icec','clobber')
      !call xcaget(t,tii,0)
      !call ncdraft(trim(ncfil),t,'T','')
      !call xcaget(waves,wamp,0)
      !call ncdraft(trim(ncfil),waves,'WAMP','')
      !call xcaget(d,dfloe,0)
      !call ncdraft(trim(ncfil),d,'dfloe','')

      ! Tile dfloe
      !call xcaput(floesize,dfloe(1-nbdy,1-nbdy),1)
 
      end subroutine waves_advect


      real function lambda(period)
      ! Calculates the waves-in-ice wavelength as a function of the period
      ! and possibly the ice thickness.
      use mod_evp
      implicit none
      real, intent(in) :: period
      lambda = gravit*period**2/(2.d0*3.141592) ! Deep water dispersion relation [m]
      end function lambda


      real function wc(h,l)
      ! Calculates the critical wave height over which ice fails in flexion and cracks.
      use mod_evp
      implicit none
      real, intent(in) :: h,l
      real             :: wc1,wc2
      real, parameter  :: rhoice=930.     ! Sea ice density       [kg/m^3]
      real, parameter  :: sigmac=0.9e9    ! Flexural ice strength [Pa]
      real, parameter  :: epsilonc=5.0e-5 ! Yield flexural strain [-]
      real, parameter  :: mu=0.6          ! Fatigue parameter     [-]
      wc1 = 2.d0*3.141592*h**2*mu*sigmac/(3.d0*gravit*rhoice*l**2)  ! Stress yield amplitude [m]
      wc2 = l**2*mu*epsilonc/(2.d0*3.141592**2*h)  ! Stress yield amplitude [m]
      wc  = min(wc1,wc2)
      end function wc

      end module mod_waves

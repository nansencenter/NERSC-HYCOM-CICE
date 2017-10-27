#if defined(ROW_LAND)
#define SEA_P .true.
#define SEA_U .true.
#define SEA_V .true.
#elif defined(ROW_ALLSEA)
#define SEA_P allip(j).or.ip(i,j).ne.0
#define SEA_U alliu(j).or.iu(i,j).ne.0
#define SEA_V alliv(j).or.iv(i,j).ne.0
#else
#define SEA_P ip(i,j).ne.0
#define SEA_U iu(i,j).ne.0
#define SEA_V iv(i,j).ne.0
#endif

module mod_hycom_fabm
#ifdef _FABM_
   use fabm
   use fabm_config

   use mod_xc         ! HYCOM communication interface
   use mod_cb_arrays  ! HYCOM saved arrays

   implicit none

   private

   public hycom_fabm_configure, hycom_fabm_initialize, hycom_fabm_update, hycom_fabm_read_relax
   public fabm_surface_state, fabm_bottom_state

   class (type_model), pointer, save, public :: fabm_model => null()
   real, allocatable :: swflx_fabm(:, :)
   real, allocatable :: bottom_stress(:, :)
   logical, allocatable :: mask(:, :, :)
   integer, allocatable :: kbottom(:, :)
   real, allocatable :: h(:, :, :)
   real, allocatable :: fabm_surface_state(:, :, :, :)
   real, allocatable :: fabm_bottom_state(:, :, :, :)

   logical :: do_interior_sources, do_bottom_sources, do_surface_sources, do_vertical_movement

contains

    subroutine hycom_fabm_configure()
      integer :: configuration_method
      logical :: file_exists
      integer, parameter :: namlst = 9000
      integer :: ios
      character(len=*), parameter :: path = '../hycom_fabm.nml'
      namelist /hycom_fabm/ do_interior_sources, do_bottom_sources, do_surface_sources, do_vertical_movement

      ! Read coupler configuration
      do_interior_sources = .true.
      do_bottom_sources = .true.
      do_surface_sources = .true.
      do_vertical_movement = .true.
      inquire(file='../hycom_fabm.nml', exist=file_exists)
      if (file_exists) then
        write (*,*) 'Reading HYCOM-FABM coupler configuration from '//path
        open(namlst, file=path, action='read', status='old', iostat=ios)
        if (ios/=0) stop 'error opening hycom_fabm.nml'
        read(namlst, nml=hycom_fabm, iostat=ios)
        if (ios/=0) stop 'error reading hycom_fabm.nml'
        close(namlst)
      end if

      configuration_method = 1
      inquire(file='../fabm.yaml', exist=file_exists)
      if (.not.file_exists) then
        inquire(file='../fabm.nml', exist=file_exists)
        if (file_exists) configuration_method = 0
      end if

      select case (configuration_method)
      case (0)
        ! From namelists in fabm.nml - DEPRECATED!!
        fabm_model => fabm_create_model_from_file(namlst, '../fabm.nml')
      case (1)
        ! From YAML file fabm.yaml
        allocate(fabm_model)
        call fabm_create_model_from_yaml_file(fabm_model, '../fabm.yaml')
      end select
    end subroutine hycom_fabm_configure

    subroutine hycom_fabm_initialize()

      integer :: j, k

        allocate(swflx_fabm(ii, jj))
        allocate(bottom_stress(ii, jj))
        allocate(mask(ii, jj, kk))
        allocate(kbottom(ii, jj))
        allocate(h(ii, jj, kk))
        allocate(fabm_surface_state(ii, jj, 2, size(fabm_model%surface_state_variables)))
        allocate(fabm_bottom_state(ii, jj, 2, size(fabm_model%bottom_state_variables)))

        ! Provide extents of the spatial domain (number of layers nz for a 1D column)
        call fabm_set_domain(fabm_model, ii, jj, kk, baclin)

        ! Send mask - see SEA_P preprocessor macro in trcupd.F
        call fabm_set_mask(fabm_model, mask, mask(:, :, 1))

        ! Specify vertical index of surface and bottom
        call fabm_model%set_surface_index(1)
        call fabm_model%set_bottom_index(kbottom)

        call fabm_model%link_interior_data(standard_variables%cell_thickness, h(1:ii, 1:jj, 1:kk))
        call fabm_model%link_horizontal_data(standard_variables%surface_downwelling_shortwave_flux, swflx_fabm(1:ii, 1:jj))
        call fabm_model%link_horizontal_data(standard_variables%bottom_stress, bottom_stress(1:ii, 1:jj))

        call update_fabm_data(1, initializing=.true.)  ! initialize the entire column of wet points, including thin layers

        ! Check whether FABM has all dependencies fulfilled
        ! (i.e., whether all required calls for fabm_link_*_data have been made)
        call fabm_check_ready(fabm_model)

        ! Initialize the tracers
        ! This sets the values of arrays sent to fabm_link_interior_state_data, in this case interior_state.
        tracer = 0
        do k=1,kk
          do j=1,jj
              call fabm_initialize_state(fabm_model, 1, ii, j, k)
          end do
        end do
        do j=1,jj
            call fabm_initialize_bottom_state(fabm_model, 1, ii, j)
            call fabm_initialize_surface_state(fabm_model, 1, ii, j)
        end do

        ! Copy state from time step = 1 to time step = 2
        tracer(:, :, :, 2, :) = tracer(:, :, :, 1, :)
        fabm_bottom_state(:, :, 2, :) = fabm_bottom_state(:, :, 1, :)
        fabm_surface_state(:, :, 2, :) = fabm_surface_state(:, :, 1, :)
    end subroutine hycom_fabm_initialize

    subroutine hycom_fabm_read_relax()
    end subroutine hycom_fabm_read_relax

    subroutine hycom_fabm_update(m, n, ibio)
      integer, intent(in) :: m, n, ibio
      integer :: i, k, j, ivar

      real :: extinction(ii)
      real :: sms(ii, size(fabm_model%state_variables))
      real :: flux(ii, size(fabm_model%state_variables))
      real :: sms_bt(ii, size(fabm_model%bottom_state_variables))
      real :: sms_sf(ii, size(fabm_model%surface_state_variables))

      write (*,*) 'hycom_fabm_update'
!
! --- leapfrog time step.
!
      ! TODO: send m or n state for computation of source terms? Leapfrog would need m, ECOSMO seems to do n
      ! Note: if we use n, then the bottom, surface and interior operations below each perform their own update
      ! before the next operation comes in, and that next one will use the updated value. This is in effect operator splitting...
      call update_fabm_data(n, initializing=.false.)  ! skipping thin layers

      call check_state('before vertical_movement', n)

      if (do_vertical_movement) then
        call vertical_movement(n, n, delt1)
        call check_state('after vertical_movement', n)
      end if

#ifdef FABM_CHECK_NAN
    do j=1,jj
        do i=1,ii
            if (SEA_P) then
                if (isnan(swflx_fabm(i,j))) then
                    write (*,*) 'NaN in swflx_fabm:', swflx_fabm(i,j), swflx (i,j,l0),w0,swflx (i,j,l1),w1,swflx (i,j,l2),w2,swflx (i,j,l3),w3
                    stop
                end if
            end if
        end do
    end do
#endif

      do k=1,kk
        do j=1,jj
          call fabm_get_light_extinction(fabm_model, 1, ii, j, k, extinction)
#ifdef FABM_CHECK_NAN
          if (any(isnan(extinction))) then
            write (*,*) 'NaN in extinction:', extinction
            stop
          end if
#endif
        end do
      end do

      ! Update light field
      do i=1,ii
        do j=1,jj
            call fabm_get_light(fabm_model, 1, kk, i, j)
        end do
      end do

      ! Compute bottom source terms
      if (do_bottom_sources) then
      do j=1,jj
        flux = 0
        sms_bt = 0
        call fabm_do_bottom(fabm_model, 1, ii, j, flux, sms_bt)
        do ivar=1,size(fabm_model%bottom_state_variables)
          fabm_bottom_state(1:ii, j, n, ivar) = fabm_bottom_state(1:ii, j, n, ivar) + delt1 * sms_bt(1:ii, ivar)
        end do
        do i=1,ii
          if (SEA_P) then
            tracer(i, j, kbottom(i, j), n, :) = tracer(i, j, kbottom(i, j), n, :) + delt1 * flux(i, :)/h(i, j, kbottom(i, j))
#ifdef FABM_CHECK_NAN
            if (any(isnan(tracer(i, j, kbottom(i, j), n, :)))) then
              write (*,*) 'NaN after do_bottom:', tracer(i, j, kbottom(i, j), n, :), flux(i, :), h(i, j, kbottom(i, j))
              stop
            end if
#endif
          end if
        end do
      end do
      call check_state('after bottom sources', n)
      end if

      ! Compute surface source terms
      if (do_surface_sources) then
      do j=1,jj
        flux = 0
        sms_sf = 0
        call fabm_do_surface(fabm_model, 1, ii, j, flux, sms_sf)
        do ivar=1,size(fabm_model%surface_state_variables)
          fabm_surface_state(1:ii, j, n, ivar) = fabm_surface_state(1:ii, j, n, ivar) + delt1 * sms_sf(1:ii, ivar)
        end do
        do i=1,ii
          if (SEA_P) then
            tracer(i, j, 1, n, :) = tracer(i, j, 1, n, :) + delt1 * flux(i, :)/h(i, j, 1)
#ifdef FABM_CHECK_NAN
            if (any(isnan(tracer(i, j, 1, n, :)))) then
              write (*,*) 'NaN after do_surface:', tracer(i, j, 1, n, :), flux(i, :), h(i, j, 1)
              stop
            end if
#endif
          end if
        end do
      end do
      call check_state('after surface sources', n)
      end if

      ! Compute source terms and update state
      if (do_interior_sources) then
      do k=1,kk
        do j=1,jj
            sms = 0
            call fabm_do(fabm_model, 1, ii, j, k, sms)
            do ivar=1,size(fabm_model%state_variables)
               tracer(1:ii, j, k, n, ivar) = tracer(1:ii, j, k, n, ivar) + delt1 * sms(1:ii, ivar)
            end do
#ifdef FABM_CHECK_NAN
            if (any(isnan(sms))) then
              do ivar=1,size(fabm_model%state_variables)
                if (any(isnan(sms(1:ii, ivar)))) write (*,*) 'NaN in sms:',ivar,sms(1:ii, ivar)
              end do
              write (*,*) 'NaN in sms'
              do ivar=1,size(fabm_model%state_variables)
                write (*,*) 'state:',ivar,tracer(1:ii, j, k, m, ivar)
              end do
              stop
            end if
#endif
        end do
      end do
      call check_state('after interior sources', n)
      end if

      ! Copy bottom value for pelagic tracers to all layers below bottom
      ! (currently masked, but could be revived later)
      do j=1,jj
        do i=1,ii
          if (SEA_P) then
            do k=kbottom(i, j)+1, kk
               tracer(i, j, k, n, :) = tracer(i, j, kbottom(i, j), n, :)
            end do
          end if
        end do
      end do

    end subroutine hycom_fabm_update

    subroutine check_state(location, n)
      use, intrinsic :: ieee_arithmetic
      character(len=*), intent(in) :: location
      integer, intent(in) :: n

      logical, parameter :: repair = .false.
      logical :: valid_int, valid_sf, valid_bt

      integer :: j, k, ivar

      do k=1,kk
        do j=1,jj
          call fabm_check_state(fabm_model, 1, ii, j, k, repair, valid_int)
          if (.not.(valid_int.or.repair)) then
            write (*,*) 'Invalid interior state '//location
            stop
          end if
        end do
      end do
      do j=1,jj
        call fabm_check_surface_state(fabm_model, 1, ii, j, repair, valid_sf)
        call fabm_check_bottom_state(fabm_model, 1, ii, j, repair, valid_bt)
        if (.not.(valid_sf.and.valid_bt).and..not.repair) then
          write (*,*) 'Invalid interface state '//location
          stop
        end if
      end do

      do ivar=1,size(fabm_model%state_variables)
        if (.not.all(ieee_is_finite(tracer(1:ii, 1:jj, 1:kk, n, ivar)))) then
          write (*,*) location, 'Interior state variable not finite:', ivar, 'range', minval(tracer(1:ii, 1:jj, 1:kk, n, ivar)), maxval(tracer(1:ii, 1:jj, 1:kk, n, ivar))
          stop
        end if
      end do
    end subroutine check_state

    subroutine vertical_movement(n, m, timestep)
      integer, intent(in) :: n, m
      real, intent(in) :: timestep

      real :: w(ii, kk, size(fabm_model%state_variables))
      real :: flux(ii, 0:kk)
      integer :: i, j, k, ivar, kabove
      real, parameter :: epsilon = 1e-8

      do j=1,jj
        ! Get vertical velocities per tracer (m/s, > 0 for floating, < 0  for sinking)
        do k=1,kk
          call fabm_get_vertical_movement(fabm_model, 1, ii, j, k, w(1:ii, k, :))
        end do

        do ivar=1,size(fabm_model%state_variables)
          ! Compute tracer flux over layer interfaces
          flux = 0
          do k=1,kk
            do i=1,ii
              if (w(i, k, ivar) > 0) then
                ! Floating: move tracer upward over top interface of the layer (flux > 0)
                flux(i, k-1) = flux(i, k-1) + min((1-epsilon)*h(i, j, k)/timestep, w(i, k, ivar))*tracer(i, j, k, m, ivar)
              else
                ! Sinking: move tracer downward over bottom interface of the layer (flux < 0)
                flux(i, k) = flux(i, k) + max(-(1-epsilon)*h(i, j, k)/timestep, w(i, k, ivar))*tracer(i, j, k, m, ivar)
              end if
            end do ! i
          end do ! k

          ! Update state
          do i=1,ii
            kabove = 0
            do k=1,kbottom(i, j)-1
              if (h(i, j, k) > 0) kabove = k
              if (flux(i, k) /= 0) then
                ! non-zero flux across interface
                if (h(i, j, k+1) == 0) then
                  ! layer below is collapsed (height = 0) - move flux to next interface
                  flux(i, k+1) = flux(i, k+1) + flux(i, k)
                else
                  ! layer below has non-zero height
                  tracer(i, j, kabove, n, ivar) = tracer(i, j, kabove, n, ivar) + flux(i, k)*timestep/h(i, j, kabove)
                  tracer(i, j, k+1, n, ivar) = tracer(i, j, k+1, n, ivar) - flux(i, k)*timestep/h(i, j, k+1)
                end if
              end if
            end do
          end do
        end do ! ivar
      end do ! j

    end subroutine vertical_movement

    subroutine update_fabm_data(index, initializing)
        integer, intent(in) :: index
        logical, intent(in) :: initializing

        integer :: i, j, k
        integer :: ivar

        real, parameter :: rho_0 = 1025.   ! [kg/m3]

        ! Update cell thicknesses (m)
        h(:, :, :) = dp(1:ii, 1:jj, 1:kk, index)/onem

        ! Update mask and kbottom
        kbottom = 0
        mask = .false.
        do j=1,jj
            do i=1,ii
              if (SEA_P) then
                 if (initializing) then
                   kbottom(i, j) = kk
                 else
                   do k=kk,1,-1
                     if (h(i, j, k)>0.1) exit
                   end do
                   kbottom(i, j) = max(k, 2)
                 end if
                 mask(i, j, 1:kbottom(i, j)) = .true.
              end if
            end do
        end do

        if (.not.initializing) then
          call fabm_update_time(fabm_model, real(nstep))

          ! Compute downwelling shortwave (from thermf.F)
          do j=1,jj
              do i=1,ii
                  if (SEA_P) then
                      if (natm.eq.2) then
                        swflx_fabm(i,j)=swflx (i,j,l0)*w0+swflx (i,j,l1)*w1
                      else
                        swflx_fabm(i,j)=swflx (i,j,l0)*w0+swflx (i,j,l1)*w1+swflx (i,j,l2)*w2+swflx (i,j,l3)*w3
                      endif !natm
                  end if
              end do
          end do

          ! Compute bottom stress (Pa)
          do j=1,jj
            do i=1,ii
              if (SEA_P) then
                bottom_stress(i, j) = ustarb(i, j)*ustarb(i, j)*rho_0
              end if
            end do
          end do
        end if

        ! Send pointers to state variable data to FABM
        do ivar=1,size(fabm_model%state_variables)
          call fabm_link_interior_state_data(fabm_model, ivar, tracer(1:ii, 1:jj, 1:kk, index, ivar))
        end do
        do ivar=1,size(fabm_model%bottom_state_variables)
          call fabm_link_bottom_state_data(fabm_model, ivar, fabm_bottom_state(1:ii, 1:jj, index, ivar))
        end do
        do ivar=1,size(fabm_model%surface_state_variables)
          call fabm_link_surface_state_data(fabm_model, ivar, fabm_surface_state(1:ii, 1:jj, index, ivar))
        end do

        ! Transfer pointer to environmental data
        ! Do this for all variables on FABM's standard variable list that the model can provide.
        ! For this list, visit http://fabm.net/standard_variables
        call fabm_model%link_interior_data(standard_variables%temperature, temp(1:ii, 1:jj, 1:kk, index))
        call fabm_model%link_interior_data(standard_variables%practical_salinity, saln(1:ii, 1:jj, 1:kk, index))
    end subroutine update_fabm_data
#endif
end module

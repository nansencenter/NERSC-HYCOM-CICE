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

   implicit none

   private

   public fabm_create_model_from_yaml_file, hycom_fabm_initialize, hycom_fabm_update, hycom_fabm_read_relax

   type (type_model), save, public :: fabm_model
   real, allocatable :: swflx_fabm(:, :)
   logical, allocatable :: mask(:, :, :)
   integer, allocatable :: kbottom(:, :)
   real, allocatable :: h(:, :, :)

   real,    parameter   :: onem=9806.0          ! g/thref

contains

    subroutine hycom_fabm_initialize()
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays

      integer :: j, k

        allocate(swflx_fabm(ii, jj))
        allocate(mask(ii, jj, kk))
        allocate(kbottom(ii, jj))
        allocate(h(ii, jj, kk))

        ! Provide extents of the spatial domain (number of layers nz for a 1D column)
        call fabm_set_domain(fabm_model, ii, jj, kk)

        ! Send mask - see SEA_P preprocessor macro in trcupd.F
        call fabm_set_mask(fabm_model, mask, mask(:, :, 1))

        ! Specify vertical index of surface and bottom
        call fabm_model%set_surface_index(1)
        call fabm_model%set_bottom_index(kbottom)

        call fabm_model%link_interior_data(standard_variables%cell_thickness, h(1:ii, 1:jj, 1:kk))
        call fabm_model%link_horizontal_data(standard_variables%surface_downwelling_shortwave_flux, swflx_fabm(1:ii, 1:jj))

        call update_fabm_data(1)

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

        tracer(:, :, :, 2, :) = tracer(:, :, :, 1, :)
        write (*,*) 'tracer min, max, sum:', minval(tracer), maxval(tracer), sum(tracer)
    end subroutine hycom_fabm_initialize

    subroutine hycom_fabm_read_relax()
    end subroutine

    subroutine hycom_fabm_update(m, n, ibio)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays

      integer, intent(in) :: m, n, ibio
      integer :: i, k, j, ivar

      real :: extinction(ii)
      real :: dy(ii, size(fabm_model%state_variables))
      write (*,*) 'hycom_fabm_update'
!
! --- leapfrog time step.
!
      ! TODO: send m or n state for computation of source terms? Leapfrog would need m, ECOSMO seems to do n
      call update_fabm_data(n)

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

      do k=1,kk
        do j=1,jj
          call fabm_get_light_extinction(fabm_model, 1, ii, j, k, extinction)
          if (any(isnan(extinction))) then
            write (*,*) 'NaN in extinction:', extinction
            stop
          end if
        end do
      end do

      ! Update light field
      do i=1,ii
        do j=1,jj
            call fabm_get_light(fabm_model, 1, kk, i, j)
        end do
      end do

      ! Compute source terms and update state
      do k=1,kk
        do j=1,jj
            dy = 0
            call fabm_do(fabm_model, 1, ii, j, k, dy)
            do ivar=1,size(fabm_model%state_variables)
               if (any(isnan(dy(1:ii, ivar)))) write (*,*) 'NaN in dy:',ivar,dy(1:ii, ivar)
               tracer(1:ii, j, k, n, ivar) = tracer(1:ii, j, k, n, ivar) + delt1 * dy(1:ii, ivar)
            end do
            if (any(isnan(dy))) then
              write (*,*) 'NaN in dy'
              stop
            end if
        end do
      end do
      stop
    end subroutine hycom_fabm_update

    subroutine update_fabm_data(index)
        use mod_cb_arrays  ! HYCOM saved arrays

        integer, intent(in) :: index

        integer :: i, j, k
        integer :: ivar

        ! TODO: update mask and kbottom
        kbottom = kk
        mask = .false.
        do j=1,jj
            do i=1,ii
               if (SEA_P) mask(i, j, :) = .true.
            end do
        end do

        ! Update cell thicknesses (m)
        h(:, :, :) = dp(1:ii, 1:jj, 1:kk, index)/onem

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

        ! Send pointers to state variable data to FABM
        do ivar=1,size(fabm_model%state_variables)
          call fabm_link_interior_state_data(fabm_model, ivar, tracer(1:ii, 1:jj, 1:kk, index, ivar))
        end do
        do ivar=1,size(fabm_model%bottom_state_variables)
          !call fabm_link_bottom_state_data(fabm_model, ivar, ???)
        end do
        do ivar=1,size(fabm_model%surface_state_variables)
          !call fabm_link_surface_state_data(fabm_model, ivar, ???)
        end do

        ! Transfer pointer to environmental data
        ! Do this for all variables on FABM's standard variable list that the model can provide.
        ! For this list, visit http://fabm.net/standard_variables
        call fabm_model%link_interior_data(standard_variables%temperature, temp(1:ii, 1:jj, 1:kk, index))
        call fabm_model%link_interior_data(standard_variables%practical_salinity, saln(1:ii, 1:jj, 1:kk, index))
    end subroutine
#endif
end module
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
   use fabm_types, only: attribute_length, output_none
   use fabm_standard_variables, only: type_global_standard_variable

   use mod_xc         ! HYCOM communication interface
   use mod_cb_arrays  ! HYCOM saved arrays
#ifdef CPL_OASIS_HYCOM
   use mod_cpl_oasis_init
#endif
   implicit none

   private

   public hycom_fabm_configure, hycom_fabm_allocate, hycom_fabm_initialize, hycom_fabm_initialize_state, hycom_fabm_update
   public hycom_fabm_relax_init, hycom_fabm_relax_rewind, hycom_fabm_relax_skmonth, hycom_fabm_relax_read, hycom_fabm_relax
   public hycom_fabm_input_init, hycom_fabm_input_update
   public hycom_fabm_nest_next, hycom_fabm_nest_read, hycom_fabm_nest_update
   public hycom_fabm_allocate_mean_output, hycom_fabm_zero_mean_output, hycom_fabm_increment_mean_output, hycom_fabm_end_mean_output, hycom_fabm_write_mean_output
   public fabm_surface_state, fabm_bottom_state

   class (type_model), pointer, save, public :: fabm_model => null()
   real, allocatable :: swflx_fabm(:, :)
   real, allocatable :: wspd_fabm(:,:)
   real, allocatable :: bottom_stress(:, :)
   real, allocatable :: atmco2(:), atmco2_fabm(:,:)
   real :: bottom_stress_keep 
   logical, allocatable :: mask(:, :, :, :)
   integer, allocatable :: kbottom(:, :, :)
   integer :: nbottom
   real :: hbottom
   real :: wndstr,strspd
   real, allocatable :: h(:, :,:),delZ(:),codepth(:,:,:),cotemp(:,:,:),cosal(:,:,:),codens(:,:,:)
   real, allocatable :: hriver(:, :)
   real, allocatable, target :: fabm_surface_state(:, :, :, :)
   real, allocatable, target :: fabm_bottom_state(:, :, :, :)
   real, allocatable :: fabm_surface_state_old(:, :, :)
   real, allocatable :: fabm_bottom_state_old(:, :, :)

   logical :: do_interior_sources, do_bottom_sources, do_surface_sources, do_vertical_movement, do_check_state
   integer, save :: current_time_index = -1

   type type_horizontal_output
      class (type_external_variable), pointer :: metadata => null()
      real, pointer :: data2d(:,:) => null()
      real, pointer :: data3d(:,:,:) => null()
      real, allocatable :: mean(:,:)
      type (type_horizontal_output), pointer :: next => null()
   end type
   type (type_horizontal_output), pointer, save :: first_horizontal_output => null()

   type type_interior_output
      class (type_external_variable), pointer :: metadata => null()
      real, pointer :: data3d(:,:,:) => null()
      real, pointer :: data4d(:,:,:,:) => null()
      real, allocatable :: mean(:,:,:)
      type (type_interior_output), pointer :: next => null()
   end type
   type (type_interior_output), pointer, save :: first_interior_output => null()

   integer, parameter :: role_prescribe = 0
   integer, parameter :: role_river = 1
   integer, save :: next_unit = 940  !ASJUN18 - increased since 218 is reserved in hycom

   integer, allocatable :: hycom_fabm_relax(:)

   integer :: m_clim0, m_clim1, m_clim2, m_clim3
   integer :: l_clim0, l_clim1, l_clim2, l_clim3
   type type_input
      integer :: file_unit = -1

      integer :: roleriver = role_prescribe
      integer :: ivariable = -1           ! state variable index (only used if role is role_river)
      real, allocatable :: data_ip(:,:,:), data_src(:,:,:,:)
      type (type_input), pointer :: next => null()
      integer :: mrec = -1
   end type type_input
   type (type_input), pointer, save :: first_input => null()

   integer, allocatable :: istate_dev(:)
   integer :: nested_number
   character(len=attribute_length), allocatable :: nested_name_dev(:)
   real, allocatable :: nested_data_dev(:,:,:,:,:) ! i,j,k,mn,state
   character(len=attribute_length) :: nested_variables(256)
   logical :: nested_bio

   integer :: pCO2unit, yCO2init, nyearCO2,modelyear,modelmonth
   real    :: co2_seasonality(12),modelday,modeltime,pair
   real    :: dew,atmco2_0,atmco2_1,atmco2_2,atmco2_3
contains

    subroutine hycom_fabm_configure()
      integer :: configuration_method
      logical :: file_exists
      integer, parameter :: namlst = 9000
      integer :: ios, ivar, istate, nestn
      character(len=*), parameter :: path = '../hycom_fabm.nml'
      namelist /hycom_fabm/ do_interior_sources, do_bottom_sources, do_surface_sources, do_vertical_movement, nested_variables

      ! Read coupler configuration
      do_interior_sources = .true.
      do_bottom_sources = .true.
      do_surface_sources = .true.
      do_vertical_movement = .true.
      do_check_state = .false.
      nested_variables = ''
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

      nested_number=0
      do ivar=1, size(nested_variables)
        if (nested_variables(ivar) .ne. '') nested_number=nested_number+1
        if (nested_variables(ivar) == '') cycle
      enddo

      nested_bio = .false. ! by default, biogeochemical nesting is turned off. If "hycom_fabm.nml" passes nested variables, sets .true.
      if ( nested_number > 0 ) nested_bio = .true.

      allocate(istate_dev(nested_number))
      do nestn = 1,nested_number
        do istate=1, size(fabm_model%state_variables)
        if (nested_variables(nestn) == fabm_model%state_variables(istate)%name) istate_dev(nestn) = istate
        enddo
      enddo

    ! read atmospheric CO2 time-series
    pCO2unit = 2640
    yCO2init = 1948
    nyearCO2 = 71 ! last data year = 2018
    co2_seasonality = [1.8552,2.7411,3.7847,4.6522,4.2706,1.1206,-4.3468,-8.6286,-7.9663,-3.7896,1.8954,3.0054]
      
    end subroutine hycom_fabm_configure

    subroutine hycom_fabm_allocate()

        allocate(swflx_fabm(ii, jj))
        allocate(wspd_fabm(ii,jj))
        allocate(bottom_stress(ii, jj))
        allocate(mask(ii, jj, kk, 2))
        allocate(kbottom(ii, jj, 2))
        allocate(h(ii, jj, kk))
        allocate(delZ(kk))
        allocate(codepth(ii,jj,kk))
        allocate(cotemp(ii,jj,kk))
        allocate(cosal(ii,jj,kk))
        allocate(codens(ii,jj,kk))
        allocate(hriver(ii, jj))
        allocate(fabm_surface_state(1-nbdy:idm+nbdy, 1-nbdy:jdm+nbdy, 2, size(fabm_model%surface_state_variables)))
        allocate(fabm_bottom_state(1-nbdy:idm+nbdy, 1-nbdy:jdm+nbdy, 2, size(fabm_model%bottom_state_variables)))
        allocate(fabm_surface_state_old(1-nbdy:idm+nbdy, 1-nbdy:jdm+nbdy, size(fabm_model%surface_state_variables)))
        allocate(fabm_bottom_state_old(1-nbdy:idm+nbdy, 1-nbdy:jdm+nbdy, size(fabm_model%bottom_state_variables)))
        allocate(nested_data_dev(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kknest,2,nested_number))
        allocate(atmco2(nyearCO2))
        allocate(atmco2_fabm(ii, jj))
    end subroutine hycom_fabm_allocate

    subroutine hycom_fabm_initialize()

      integer :: j, k, ivar
      type (type_interior_output),   pointer :: last_interior_output
      type (type_horizontal_output), pointer :: last_horizontal_output

      integer :: yCO2

        ! Provide extents of the spatial domain (number of layers nz for a 1D column)
        call fabm_set_domain(fabm_model, ii, jj, kk, baclin)

        ! Specify vertical index of surface (constant during simulation)
        call fabm_model%set_surface_index(1)

        call fabm_model%link_interior_data(standard_variables%cell_thickness, h(1:ii, 1:jj, 1:kk))
        call fabm_model%link_horizontal_data(standard_variables%surface_downwelling_shortwave_flux, swflx_fabm(1:ii, 1:jj))
        call fabm_model%link_horizontal_data(standard_variables%wind_speed,wspd_fabm(1:ii, 1:jj))
        call fabm_model%link_horizontal_data(standard_variables%mole_fraction_of_carbon_dioxide_in_air,atmco2_fabm(1:ii,1:jj))
        call fabm_model%link_horizontal_data(standard_variables%bottom_stress, bottom_stress(1:ii, 1:jj))
        call fabm_model%link_scalar(type_global_standard_variable(name='time_step', units='s'), delt1)

        call update_fabm_data(1, initializing=.true.)  ! initialize the entire column of wet points, including thin layers

        ! Check whether FABM has all dependencies fulfilled
        ! (i.e., whether all required calls for fabm_link_*_data have been made)
        call fabm_check_ready(fabm_model)

        last_interior_output => null()
        do ivar=1, size(fabm_model%state_variables)
          if (add_interior_output(fabm_model%state_variables(ivar))) &
            last_interior_output%data4d => tracer(1:ii, 1:jj, 1:kk, :, ivar)
        end do
        do ivar=1, size(fabm_model%diagnostic_variables)
          if (add_interior_output(fabm_model%diagnostic_variables(ivar))) &
            last_interior_output%data3d => fabm_get_interior_diagnostic_data(fabm_model, ivar)
        end do

        last_horizontal_output => null()
        do ivar=1, size(fabm_model%surface_state_variables)
          if (add_horizontal_output(fabm_model%surface_state_variables(ivar))) &
            last_horizontal_output%data3d => fabm_surface_state(1:ii, 1:jj, :, ivar)
        end do
        do ivar=1, size(fabm_model%bottom_state_variables)
          if (add_horizontal_output(fabm_model%bottom_state_variables(ivar))) &
            last_horizontal_output%data3d => fabm_bottom_state(1:ii, 1:jj, :, ivar)
        end do
        do ivar=1, size(fabm_model%horizontal_diagnostic_variables)
          if (add_horizontal_output(fabm_model%horizontal_diagnostic_variables(ivar))) &
            last_horizontal_output%data2d => fabm_get_horizontal_diagnostic_data(fabm_model, ivar)
        end do

        open(unit=pCO2unit, file='pCO2a_1948_2018',form='formatted')
        do yCO2=1,nyearCO2
          read(pCO2unit,*) atmco2(yCO2)
        end do
        close(pCO2unit) 
    contains

      function add_interior_output(variable) result(saved)
        class (type_external_variable), target, intent(in) :: variable
        logical :: saved

        type (type_interior_output), pointer :: interior_output

        saved = variable%output /= output_none
        if (.not.saved) return
        allocate(interior_output)
        interior_output%metadata => variable
        if (associated(last_interior_output)) then
          last_interior_output%next => interior_output
        else
          first_interior_output => interior_output
        end if
        last_interior_output => interior_output
      end function

      function add_horizontal_output(variable) result(saved)
        class (type_external_variable), target, intent(in) :: variable
        logical :: saved

        type (type_horizontal_output), pointer :: horizontal_output

        saved = variable%output /= output_none
        if (.not.saved) return
        allocate(horizontal_output)
        horizontal_output%metadata => variable
        if (associated(last_horizontal_output)) then
          last_horizontal_output%next => horizontal_output
        else
          first_horizontal_output => horizontal_output
        end if
        last_horizontal_output => horizontal_output
      end function

    end subroutine hycom_fabm_initialize

    subroutine hycom_fabm_initialize_state()
      integer :: k, j

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
    end subroutine hycom_fabm_initialize_state

    subroutine hycom_fabm_relax_init()
      use mod_za  ! HYCOM I/O interface

      integer :: ivar, k
      logical :: file_exists
      character preambl(5)*79

      ! Allocate array to holds units for relaxation files of every pelagic state variable
      allocate(hycom_fabm_relax(size(fabm_model%state_variables)))

      ! Default: no relaxation
      hycom_fabm_relax = -1

      if (mnproc.eq.1) write (lp,*) 'Looking for relaxation data for pelagic FABM state variables...'
      do ivar=1,size(fabm_model%state_variables)
        ! Check for existence of a file named "relax.<FABMNAME>.a". If present, this will contain the relaxation field (one variable; all k levels)
        inquire(file=trim(flnmforw)//'relax.'//trim(fabm_model%state_variables(ivar)%name)//'.a', exist=file_exists)
        if (file_exists) then
          if (mnproc.eq.1) write (lp,*) '  - '//trim(fabm_model%state_variables(ivar)%name)//': ON, ' &
            //trim(flnmforw)//'relax.'//trim(fabm_model%state_variables(ivar)%name)//'.a was found.'

          ! Relaxation file exist; assign next available unit.
          hycom_fabm_relax(ivar) = next_unit
          next_unit = next_unit + 1

          ! Open binary file (.a)
          call zaiopf(trim(flnmforw)//'relax.'//trim(fabm_model%state_variables(ivar)%name)//'.a', 'old', hycom_fabm_relax(ivar))

          ! Open metadata (.b)
          if (mnproc.eq.1) then  ! .b file from 1st tile only
            open (unit=uoff+hycom_fabm_relax(ivar),file=trim(flnmforw)//'relax.'//trim(fabm_model%state_variables(ivar)%name)//'.b', &
               status='old', action='read')
            read (uoff+hycom_fabm_relax(ivar),'(a79)') preambl
          end if !1st tile
          call preambl_print(preambl)

          ! ?? Not sure why we are reading here, copying from forfun.F
          do k=1,kk
            call hycom_fabm_rdmonthck(util1, hycom_fabm_relax(ivar), 0, is_2d=.false.)
          end do
        else
          ! Disable relaxation for this tracer
          if (mnproc.eq.1) write (lp,*) '  - '//trim(fabm_model%state_variables(ivar)%name)//': OFF, ' &
            //trim(flnmforw)//'relax.'//trim(fabm_model%state_variables(ivar)%name)//'.a not found.'
          rmutr(:,:,ivar) = 0.0
        end if
      end do
    end subroutine hycom_fabm_relax_init

    subroutine hycom_fabm_input_init(dtime, dyear0, dyear, dmonth)
      real, intent(in) :: dtime, dyear0, dyear, dmonth

      integer :: ivar
      logical :: file_exists
      type (type_input), pointer :: input

      modeltime=dtime
      call fabm_gettime()
      ! Months and slot indices for monthly climatological forcing
      ! (shared between all inputs that are defined on monthly climatological time scales)
      m_clim1=1.+mod(dtime+dyear0,dyear)/dmonth
      m_clim0=mod(m_clim1+10,12)+1
      m_clim2=mod(m_clim1,   12)+1
      m_clim3=mod(m_clim2,   12)+1
      l_clim0=1
      l_clim1=2
      l_clim2=3
      l_clim3=4

      ! Detect river forcing for pelagic state variables
      if (mnproc.eq.1) write (lp,*) 'Looking for river loadings for pelagic FABM state variables...'
      do ivar=1,size(fabm_model%state_variables)
        ! Check for existence of a file named "rivers.<FABMNAME>.a". If present, this will contain the river loadings for this variable across the entire model grid (one 2d variable; units <UNITS>*m/s)
        inquire(file=trim(flnmforw)//'rivers.'//trim(fabm_model%state_variables(ivar)%name)//'.a', exist=file_exists)
        if (file_exists) then
          if (mnproc.eq.1) write (lp,*) '  - '//trim(fabm_model%state_variables(ivar)%name)//': ON, ' &
            //trim(flnmforw)//'rivers.'//trim(fabm_model%state_variables(ivar)%name)//'.a was found.'
          input => add_input(trim(flnmforw)//'rivers.'//trim(fabm_model%state_variables(ivar)%name)//'.a', .true.)
          input%roleriver = role_river
          input%ivariable = ivar
        else
          ! Disable river input for this tracer
          if (mnproc.eq.1) write (lp,*) '  - '//trim(fabm_model%state_variables(ivar)%name)//': OFF, ' &
            //trim(flnmforw)//'rivers.'//trim(fabm_model%state_variables(ivar)%name)//'.a not found.'
        end if
      end do
    end subroutine hycom_fabm_input_init

    function add_input(path, is_2d) result(input)
      use mod_za  ! HYCOM I/O interface

      character(len=*), intent(in) :: path
      logical,          intent(in) :: is_2d
      type (type_input), pointer :: input

      integer :: k
      character preambl(5)*79

      allocate(input)
      if (is_2d) then
        allocate(input%data_ip(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1))
      else
        allocate(input%data_ip(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk))
      end if
      allocate(input%data_src(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,size(input%data_ip,3),4))
      input%file_unit = next_unit
      input%next => first_input
      first_input => input
      next_unit = next_unit + 1

      ! Open binary file (.a)
      call zaiopf(path, 'old', input%file_unit)

      ! Open metadata (.b)
      if (mnproc.eq.1) then  ! .b file from 1st tile only
        open (unit=uoff+input%file_unit, file=path(1:len(path)-2)//'.b', &
            status='old', action='read')
        read (uoff+input%file_unit,'(a79)') preambl
      end if !1st tile
      call preambl_print(preambl)
      input%mrec = 0  !ASJUN18 - initiate value of mrec

!ASJUN18 Removed since it is not doing anything.
!      ! ?? Not sure why we are reading here, copying from forfun.F
!      do k=1,size(input%data_ip, 3)
!        call hycom_fabm_rdmonthck(util1, input%file_unit, 0, is_2d)
!      end do

      call read_input(input,m_clim0,l_clim0)
      call read_input(input,m_clim1,l_clim1)
      call read_input(input,m_clim2,l_clim2)
      call read_input(input,m_clim3,l_clim3)
    end function add_input

    subroutine read_input(input, mrec, lslot)
      use mod_za  ! HYCOM I/O interface

      type (type_input), intent(inout) :: input
      integer,           intent(in)    :: mrec, lslot

      integer :: irec, k
      logical :: is_2d
      character preambl(5)*79

      if (mrec <= input%mrec) then
        ! Rewind
        if (mnproc.eq.1) then  ! .b file from 1st tile only
          rewind uoff+input%file_unit
          read (uoff+input%file_unit,'(a79)') preambl
!          read  (uoff+input%file_unit,*)
!          read  (uoff+input%file_unit,*)
!          read  (uoff+input%file_unit,*)
!          read  (uoff+input%file_unit,*)
!          read  (uoff+input%file_unit,*)
        end if
        call zaiorw(input%file_unit)
        input%mrec = 0
      end if

      do irec=input%mrec+1,mrec-1
        do k=1,size(input%data_src, 3)
          call skmonth(input%file_unit)
        end do
      end do

      is_2d = size(input%data_src,3) == 1
      do k= 1,size(input%data_src,3)
        call hycom_fabm_rdmonthck(input%data_src(1-nbdy,1-nbdy,k,lslot), input%file_unit, mrec, is_2d)
      end do

      input%mrec = mrec
    end subroutine read_input

    subroutine hycom_fabm_input_update(dtime, dyear0, dyear, dmonth)
      real, intent(in) :: dtime, dyear0, dyear, dmonth

      integer :: imonth
      type (type_input), pointer :: input
      real :: month, x, x1, w0, w1, w2, w3
      integer :: lt

      modeltime=dtime
      month=1.+mod(dtime+dyear0,dyear)/dmonth
      imonth=int(month)
      if (mnproc.eq.1) write(lp,*) 'update_inputs - month = ',month,imonth
      call xcsync(flush_lp)

      input => first_input
      do while (associated(input))
        if (imonth.ne.m_clim1) then
          m_clim1=imonth
          m_clim0=mod(m_clim1+10,12)+1
          m_clim2=mod(m_clim1,   12)+1
          m_clim3=mod(m_clim2,   12)+1
          lt = l_clim0
          l_clim0=l_clim1
          l_clim1=l_clim2
          l_clim2=l_clim3
          l_clim3=lt
          call read_input(input, m_clim3, l3)
        end if
        x=mod(month,1.)
        x1=1.-x
        w1=x1*(1.+x *(1.-1.5*x ))
        w2=x *(1.+x1*(1.-1.5*x1))
        w0=-.5*x *x1*x1
        w3=-.5*x1*x *x
        input%data_ip = input%data_src(:,:,:,l_clim0)*w0 + input%data_src(:,:,:,l_clim1)*w1 + input%data_src(:,:,:,l_clim2)*w2 + input%data_src(:,:,:,l_clim3)*w3
        input => input%next
      end do
    end subroutine hycom_fabm_input_update

    subroutine hycom_fabm_relax_rewind()
      use mod_za
      integer :: ivar

      do ivar=1,size(fabm_model%state_variables)
        if (hycom_fabm_relax(ivar) /= -1) then
          if (mnproc.eq.1) then
            rewind uoff+hycom_fabm_relax(ivar)
            read  (uoff+hycom_fabm_relax(ivar),*)
            read  (uoff+hycom_fabm_relax(ivar),*)
            read  (uoff+hycom_fabm_relax(ivar),*)
            read  (uoff+hycom_fabm_relax(ivar),*)
            read  (uoff+hycom_fabm_relax(ivar),*)
          end if
          call zaiorw(hycom_fabm_relax(ivar))
        end if
      end do
    end subroutine hycom_fabm_relax_rewind

    subroutine hycom_fabm_relax_skmonth()
      integer :: ivar, k

      do ivar=1,size(fabm_model%state_variables)
        if (hycom_fabm_relax(ivar) /= -1) then
          do k=1,kk
            call skmonth(hycom_fabm_relax(ivar))
          end do
        end if
      end do
    end subroutine hycom_fabm_relax_skmonth

    subroutine hycom_fabm_relax_read(lslot, mnth)
      use mod_za  ! HYCOM I/O interface

      integer, intent(in) :: lslot, mnth

      integer :: ivar, k

      do ivar=1,size(fabm_model%state_variables)
        if (hycom_fabm_relax(ivar) /= -1) then
          do k= 1,kk
            call hycom_fabm_rdmonthck(trwall(1-nbdy,1-nbdy,k,lslot,ivar), hycom_fabm_relax(ivar), mnth, is_2d=.false.)
          end do
        end if
      end do
    end subroutine hycom_fabm_relax_read

    subroutine hycom_fabm_rdmonthck(field, iunit, mnthck, is_2d)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      use mod_za         ! HYCOM I/O interface

      integer   iunit,mnthck
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: field
      logical, intent(in) :: is_2d

      integer   i,ios,layer,mnth
      real      denlay,hmina,hminb,hmaxa,hmaxb
      character cline*80

      call zagetc(cline,ios, uoff+iunit)
      if (ios.ne.0) then
        if (mnproc.eq.1) then
          write(lp,*)
          write(lp,*) 'error in hycom_fabm_rdmonthck - hit end of input'
          write(lp,*) 'iunit,ios = ',iunit,ios
          write(lp,*)
        end if !1st tile
        call xcstop('(hycom_fabm_rdmonthck)')
        stop '(hycom_fabm_rdmonthck)'
      end if
      if (mnproc.eq.1) write (lp,'(a)')  cline  !print input array info
      i = index(cline,'=')

      if (is_2d) then
        read (cline(i+1:),*) mnth,hminb,hmaxb
      else
        read (cline(i+1:),*) mnth,layer,denlay,hminb,hmaxb
      end if
      if (mnth.lt.1 .or. mnth.gt.12) then
        if (mnproc.eq.1) write(lp,'(/ a,i4,a /)') 'error on unit',iunit,' - not monthly relaxation data'
        call xcstop('(hycom_fabm_rdmonthck)')
        stop '(hycom_fabm_rdmonthck)'
      end if
      if (mnthck.gt.0 .and. mnth.ne.mnthck) then
        if (mnproc.eq.1) &
          write(lp,'(/ a,i4,a,2i4,a /)') 'error on unit',iunit,' - wrong relaxation month (expected,input =',mnthck,mnth,')'
        call xcstop('(hycom_fabm_rdmonthck)')
        stop '(hycom_fabm_rdmonthck)'
      end if

      if (hminb.eq.hmaxb) then  !constant field
        field(:,:) = hminb
        call zaiosk(iunit)
      else
        call zaiord(field,ip,.false., hmina,hmaxa,iunit)

        if (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or. abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
          if (mnproc.eq.1) then
          write(lp,'(/ a / a,i3 / a / a,1p3e14.6 / a,1p3e14.6 /)') &
            'error - .a and .b files not consistent:', &
            'iunit = ',iunit, &
            cline, &
            '.a,.b min = ',hmina,hminb,hmina-hminb, &
            '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
          end if !1st tile
          call xcstop('(hycom_fabm_rdmonthck)')
          stop '(hycom_fabm_rdmonthck)'
        end if
      end if

    end subroutine hycom_fabm_rdmonthck

    subroutine hycom_fabm_update(m, n, ibio)
      integer, intent(in) :: m, n, ibio
      integer :: i, k, j, ivar

      real :: extinction(ii)
      real :: sms(ii, size(fabm_model%state_variables))
      real :: flux(ii, size(fabm_model%state_variables))
      real :: sms_bt(ii, size(fabm_model%bottom_state_variables))
      real :: sms_sf(ii, size(fabm_model%surface_state_variables))
      type (type_input), pointer :: input

            if (mnproc.eq.1) write (lp,*) 'hycom_fabm_update', nstep, time
            call xcsync(flush_lp)

!
! --- leapfrog time step.
!
      ! Send state at midpoint time=t (time index m) to FABM.
      ! As per leapfrog spec, fluxes at this time are used to update the state at t-delta_t to t+delta_t (both stored at time index n)
      ! This also sets FABM's mask, which will exclude layers that are vanishingly thin at either time m or n (or both).
      call update_fabm_data(m, initializing=.false.)  ! skipping thin layers

      ! Get index of bottom layers at the next time step.
      ! For that we use time index n, which assumes dp(:,:,:,n) has already been updated!
      ! A (partial?) update seems to happen in cnuity, which is indeed called before trcupd.
      call get_mask(m, mask(:, :, :, m), kbottom(:, :, m))
      call get_mask(n, mask(:, :, :, n), kbottom(:, :, n))

      ! Store old surface/bottom state for later application of Robert-Asselin filter.
      fabm_surface_state_old = fabm_surface_state(:, :, n, :)
      fabm_bottom_state_old = fabm_bottom_state(:, :, n, :)
      ! Make sure the biogeochemical state is valid (uses clipping if necessary)
      call check_state('when entering fabm_hycom_update', current_time_index, .true.)

      ! Vertical movement (includes sinking and floating)
      if (do_vertical_movement) then
        if (do_check_state) call check_state('before vertical_movement', n, .false.)
        call vertical_movement(n, m, delt1)
        if (do_check_state) call check_state('after vertical_movement', n, .false.)
      end if

#ifdef FABM_CHECK_NAN
    do j=1,jj
        do i=1,ii
            if (SEA_P) then
                if (isnan(swflx_fabm(i,j))) then
                    write (*,*) 'NaN in swflx_fabm:', swflx_fabm(i,j), sswflx (i,j)
                    stop
                end if
            end if
        end do
    end do
#endif
!write (*,*) 'hycom_fabm_after_vertical', nstep, time
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
!write (*,*) 'hycom_fabm_after_light', nstep, time
      ! Compute bottom source terms
      if (do_bottom_sources) then
      do j=1,jj
        flux = 0
        sms_bt = 0
        call fabm_do_bottom(fabm_model, 1, ii, j, flux, sms_bt)
        do i=1,ii
          if (kbottom(i, j, n) > 0) then
            fabm_bottom_state(i, j, n, :) = fabm_bottom_state(i, j, n, :) + delt1 * sms_bt(i, :) ! update sediment layer
            if ( dp(i, j, kbottom(i,j,n), n)/onem >= 3.0 ) then ! check if the bottom layer is thicker than 10 meters, if so, apply the flux as usual (I will decrease the criteria in time)
              tracer(i, j, kbottom(i,j,n), n, :) = tracer(i, j, kbottom(i,j,n), n, :) + delt1 * flux(i, :)/dp(i, j, kbottom(i,j,n), n)*onem
            else ! in case less than 10 meters, to avoid accumulation at the bottom thin layers, distribute the flux into multiple layers that add up to > 10 meters thickness
              hbottom = 0
              nbottom = 0
              do k = kbottom(i,j,n),1,-1
                hbottom = hbottom + dp(i ,j , k, n)/onem
                if ( hbottom >= 3.0 ) exit
                nbottom = nbottom + 1
              end do
              do k = kbottom(i,j,n)-nbottom , kbottom(i,j,n) ! distribute the flux to total height, and to multiple layers
                tracer(i, j, k, n, :) = tracer(i, j, k, n, :) + delt1 * flux(i, :)/hbottom
              end do
            end if
#ifdef FABM_CHECK_NAN
            if (any(isnan(tracer(i, j, kbottom(i,j,n), n, :)))) then
              write (*,*) 'NaN after do_bottom:', tracer(i, j, kbottom(i,j,n), n, :), flux(i, :), dp(i, j, kbottom(i,j,n), n)/onem
              stop
            end if
#endif
          end if
        end do
      end do
      if (do_check_state) call check_state('after bottom sources', n, .false.)
      end if
!write (*,*) 'hycom_fabm_after_bottom', nstep, time
      ! Compute surface source terms
      if (do_surface_sources) then
      do j=1,jj
        flux = 0
        sms_sf = 0
        call fabm_do_surface(fabm_model, 1, ii, j, flux, sms_sf)
        do i=1,ii
          if (kbottom(i, j, n) > 0) then
            fabm_surface_state(i, j, n, :) = fabm_surface_state(i, j, n, :) + delt1 * sms_sf(i, :)
            tracer(i, j, 1, n, :) = tracer(i, j, 1, n, :) + delt1 * flux(i, :)/dp(i, j, 1, n)*onem
#ifdef FABM_CHECK_NAN
            if (any(isnan(tracer(i, j, 1, n, :)))) then
              write (*,*) 'NaN after do_surface:', tracer(i, j, 1, n, :), flux(i, :), dp(i, j, 1, n)/onem
              stop
            end if
#endif
          end if
        end do
      end do
      if (do_check_state) call check_state('after surface sources', n, .false.)
      end if
!write (*,*) 'hycom_fabm_after_surface', nstep, time
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
      if (do_check_state) call check_state('after interior sources', n, .false.)
      end if
!write (*,*) 'hycom_fabm_after_interior', nstep, time
      input => first_input
      do while (associated(input))
        if (input%roleriver == role_river) then
          ! River field
          do i=1,ii
            do j=1,jj
              hriver(i, j) = sum ( h(i, j, 1:5) )
              if (SEA_P) then
                if ( hriver(i,j) > 0.0 ) then
                  do k=1,5
                    tracer(i, j, k, n, input%ivariable) = tracer(i, j, k, n, input%ivariable) + delt1 * input%data_ip(i, j, 1)/hriver(i, j)
                  end do
                end if
              end if
            end do
          end do
        end if
        input => input%next
      end do
!write (*,*) 'hycom_fabm_after_correct', nstep, time
      call check_state('after hycom_fabm_update', n, .true.)
!write (*,*) 'hycom_fabm_after_check_state', nstep, time
      ! Apply the Robert-Asselin filter to the surface and bottom state.
      ! Note that RA will be applied to the pelagic tracers within mod_tsavc - no need to do it here!
!write (*,*) 'hycom_fabm_before_assalin', nstep, time
      fabm_surface_state(1:ii, 1:jj, m, :) = fabm_surface_state(1:ii, 1:jj, m, :) + 0.5*ra2fac*(fabm_surface_state_old(1:ii, 1:jj, :)+fabm_surface_state(1:ii, 1:jj, n, :)-2.0*fabm_surface_state(1:ii, 1:jj, m, :))
      fabm_bottom_state(1:ii, 1:jj, m, :) = fabm_bottom_state(1:ii, 1:jj, m, :) + 0.5*ra2fac*(fabm_bottom_state_old(1:ii, 1:jj, :)+fabm_bottom_state(1:ii, 1:jj, n, :)-2.0*fabm_bottom_state(1:ii, 1:jj, m, :))
!write (*,*) 'hycom_fabm_update_after', nstep, time
    end subroutine hycom_fabm_update

    subroutine check_state(location, index, repair)
      use, intrinsic :: ieee_arithmetic
      character(len=*), intent(in) :: location
      integer, intent(in) :: index
      logical, intent(in) :: repair

      logical :: valid_int, valid_sf, valid_bt

      integer :: i, j, k, ivar, old_index

      old_index = current_time_index
      call update_fabm_state(index)
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
        if (.not.all(ieee_is_finite(tracer(1:ii, 1:jj, 1:kk, index, ivar)))) then
          write (*,*) location, 'Interior state variable not finite:', ivar, 'range', minval(tracer(1:ii, 1:jj, 1:kk, index, ivar)), maxval(tracer(1:ii, 1:jj, 1:kk, index, ivar))
          stop
        end if
      end do

      if (repair) then
        ! FABM will have placed "missing value" for all state variables in all masked cells.
        ! However, as these can be revived later in the simulation, make sure their value is valid by
        ! copying bottom value for pelagic tracers to all layers below bottom.
        do j=1,jj
          do i=1,ii
            if (SEA_P) then
              do k=kbottom(i, j, index)+1, kk
                tracer(i, j, k, index, :) = tracer(i, j, kbottom(i, j, index), index, :)
              end do
            end if
          end do
        end do
      end if

      call update_fabm_state(old_index)
    end subroutine check_state

    subroutine vertical_movement(n, m, timestep)
      integer, intent(in) :: n, m
      real, intent(in) :: timestep

      real :: w(ii, kk, size(fabm_model%state_variables))
      real :: flux(ii, 0:kk)
      integer :: i, j, k, ivar, kabove, kb
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
                flux(i, k-1) = flux(i, k-1) + min((1-epsilon)*dp(i, j, k, n)/onem/timestep*tracer(i, j, k, n, ivar), w(i, k, ivar)*tracer(i, j, k, m, ivar))
              else
                ! Sinking: move tracer downward over bottom interface of the layer (flux < 0)
                flux(i, k) = flux(i, k) + max(-(1-epsilon)*dp(i, j, k, n)/onem/timestep*tracer(i, j, k, n, ivar), w(i, k, ivar)*tracer(i, j, k, m, ivar))
              end if
            end do ! i
          end do ! k

          ! Update state
          do i=1,ii
            kabove = 0
            do k=1,kbottom(i, j, n)-1
              if (dp(i, j, k, n) > 0) kabove = k
              if (flux(i, k) /= 0) then
                ! non-zero flux across interface
                if (dp(i, j, k+1, n)/onem < 0.1) then ! THIS IS NOT THE BEST SOLUTION BUT MAKES THINGS STABLE AT THE MOMENT
                  ! layer below is collapsed (height = 0) - move flux to next interface
                !  flux(i, k+1) = 0.!flux(i, k+1) + flux(i, k)
                flux(i, k) = 0
                else
                  ! Prevent accumulation of settling particles in thin layers 
                  if ( flux(i, k) < 0 .and. k == kbottom(i, j, n)-1 ) then ! if settling and if at the layer above the bottom
                    if ( dp(i, j, k+1, n)/onem >= 3.0 ) then ! check if the bottom layer is actually < 10 meters, if not, apply the regular flux additions
                      tracer(i, j, kabove, n, ivar) = tracer(i, j, kabove, n, ivar) + flux(i, k)*timestep/(dp(i, j, kabove, n)/onem)
                      tracer(i, j, k+1, n, ivar) = tracer(i, j, k+1, n, ivar) - flux(i, k)*timestep/(dp(i, j, k+1, n)/onem)
                      else ! if < 10 meters
                        hbottom = 0
                        nbottom = 0 
                        do kb = kbottom(i, j, n),1,-1 ! find number of layers that add up to > 10 meters, and store the total height
                          hbottom = hbottom + dp(i ,j , kb, n)/onem
                          if ( hbottom >= 3.0 ) exit
                          nbottom = nbottom + 1
                        end do
                        ! Settle the particles from kabove
                        tracer(i, j, kabove, n, ivar) = tracer(i, j, kabove, n, ivar) + flux(i, k)*timestep/(dp(i, j, kabove, n)/onem) 
                        ! and distribute that flux to multiple layers which the depths add up to > 10 meters
                        do kb = kbottom(i, j, n) - nbottom , kbottom(i, j, n)
                          tracer(i, j, kb, n, ivar) = tracer(i, j, kb, n, ivar) - flux(i, k)*timestep/hbottom
                        end do
                    end if  
                  else ! this applies to all other layers including floating particles
                    tracer(i, j, kabove, n, ivar) = tracer(i, j, kabove, n, ivar) + flux(i, k)*timestep/(dp(i, j, kabove, n)/onem)
                    tracer(i, j, k+1, n, ivar) = tracer(i, j, k+1, n, ivar) - flux(i, k)*timestep/(dp(i, j, k+1, n)/onem)
                  endif
                end if
              end if
            end do
          end do
        end do ! ivar
      end do ! j

    end subroutine vertical_movement

    subroutine get_mask(index, lmask, lkbottom)
        integer, intent(in)  :: index
        logical, intent(out) :: lmask(:, :, :)
        integer, intent(out) :: lkbottom(:, :)

        real, parameter :: h_min = 0.1
        integer :: i, j, k

        lkbottom = 0
       ! do j=1,jj
       !     do i=1,ii
       !       if (SEA_P) then
       !         do k = kk, 1, -1
       !           if (dp(i, j, k, index)/onem > h_min) exit
       !         end do
       !         kbottom(i, j) = max(k, 2)
       !       end if
       !     end do
       ! end do
        do j=1,jj       ! CAGLAR - I did it from top to bottom in order to avoid having < 0.1 m layer in the water column.
            do i=1,ii   !          Looking for other solutions for this
              if (SEA_P) then
                do k = 1,kk
                  if (dp(i, j, k, index)/onem <= h_min) exit
                  lkbottom(i, j) = k
                end do
                lkbottom(i, j) = max(lkbottom(i,j), 2)
              end if
            end do
        end do

        lmask = .false.
        do j=1,jj
            do i=1,ii
                lmask(i, j, 1:lkbottom(i, j)) = .true.
            end do
        end do
    end subroutine get_mask

    subroutine update_fabm_data(index, initializing)
        integer, intent(in) :: index
        logical, intent(in) :: initializing

        integer :: i, j, k
        integer :: ivar
        real, parameter :: rho_0 = 1025.   ! [kg/m3]
        ! Update cell thicknesses (m)
        h(:, :, :) = dp(1:ii, 1:jj, 1:kk, index)/onem

        if (initializing) then
          ! Make sure everything is unmasked, so that the state is initialized everywhere
          mask = .true.
          kbottom = kk
        end if

        if (.not.initializing) then
          call fabm_update_time(fabm_model, real(nstep))
          call fabm_gettime() 


          ! Update surface forcing and internal variables for biology 
          do j=1,jj
              do i=1,ii
   !               if (SEA_P) then
                       bottom_stress(i, j) = ustarb(i, j)*ustarb(i, j)*rho_0 !
                       swflx_fabm(i,j)= sswflx(i,j) ! ice-corrected downwelling shortwave flux

#ifdef CPL_OASIS_HYCOM
                       wspd_fabm(i,j) = cplts_recv(i,j,i2o_wspd)
#else
                       if     (flxflg.eq.6) then ! wind speed (taken from therm.F)
                          wspd_fabm(i,j) = wndocn(i,j)
                       else if (natm.eq.2) then
                          wspd_fabm(i,j) = wndspd(i,j,l0)*w0+wndspd(i,j,l1)*w1
                       else
                          wspd_fabm(i,j) =wndspd(i,j,l0)*w0+wndspd(i,j,l1)*w1+wndspd(i,j,l2)*w2+wndspd(i,j,l3)*w3    
                       end if
#endif

#ifdef CPL_OASIS_HYCOM
                       pair = cplts_recv(i,j,i2o_mslp)
#else
                       if     (mslprf .or. flxflg.eq.6) then
                          if     (natm.eq.2) then
                             pair=mslprs(i,j,l0)*w0 + mslprs(i,j,l1)*w1
                          else
                             pair=mslprs(i,j,l0)*w0 + mslprs(i,j,l1)*w1 + mslprs(i,j,l2) * w2+mslprs(i,j,l3)*w3
                          endif 
                       else
                          pair = 1013.0*100.0
                       endif
                       pair = pair * 0.01 ! convert Pa --> mBar (hPa) (result is of magnitude 1E+3) 
#endif

                       if     (natm.eq.2) then
                          dew = dewpt(i,j,l0)*w0 + dewpt(i,j,l1)*w1
                       else
                          dew = dewpt(i,j,l0)*w0 + dewpt(i,j,l1)*w1 + dewpt(i,j,l2)*w2 + dewpt(i,j,l3)*w3
                       end if
                       dew = dew - 273.15 ! Kelvin --> degC

                       atmco2_0 = atmco2( min(modelyear - yCO2init + 1 , nyearCO2 ) ) + co2_seasonality( modelmonth )
                       atmco2_1 = ( 10.**( (7.5*dew)/(237.3+dew) ) )
                       atmco2_2 = atmco2_1 / 9.81 * 10.**(-4.0)
                       atmco2_3 = pair / 9.81 * 10.**(-2.0)
                       atmco2_fabm(i,j) = atmco2_0 * (atmco2_3 - atmco2_2) * 0.997                     

!                       if (i==itest.and.j==jtest) write(*,*)'DEWPT',dewpt(i,j,l0),dewpt(i,j,l1),dewpt(i,j,l2),dewpt(i,j,l3)
                       if (i==itest.and.j==jtest) write(*,'(A4,4(1x,F8.3))')'PCO2',dew,pair,atmco2_0,atmco2_fabm(i,j)
!write(*,'(A4,4(1x,F8.3))')'PCO2',dew,pair,atmco2_0,atmco2_fabm(i,j)
                       do k=1,kk
                          delZ(k) = dp(i,j,k,index)/onem                    !
                          if(k.eq.1)then                                    !
                             codepth(i,j,k) = delZ(k)                       !
                          else                                              !
                             codepth(i,j,k) = codepth(i,j,k-1)+delZ(k)      ! water depth
                          end if                                            !
                          cotemp(i,j,k) = max(-3.999,temp(i, j, k, index))  ! water temparature
                          cosal(i,j,k)  = saln(i, j, k, index)              ! salinity
                          codens(i,j,k) = th3d(i, j, k, index)+thbase+1000. ! water density

                       end do
    !              end if
              end do
          end do
        if (mnproc.eq.1) write(lp,*) 'call update_fabm_state ', &
                                    minval(atmco2_fabm), maxval(atmco2_fabm)
        end if
        ! Transfer pointer to environmental data
        ! Do this for all variables on FABM's standard variable list that the model can provide.
        ! For this list, visit http://fabm.net/standard_variables
!        call fabm_model%link_interior_data(standard_variables%temperature, temp(1:ii, 1:jj, 1:kk, index))
!        call fabm_model%link_interior_data(standard_variables%practical_salinity, saln(1:ii, 1:jj, 1:kk, index))
        call fabm_model%link_interior_data(standard_variables%temperature,cotemp(1:ii,1:jj, 1:kk))
        call fabm_model%link_interior_data(standard_variables%practical_salinity,cosal(1:ii,1:jj,1:kk))
!        call fabm_model%link_interior_data(standard_variables%density, th3d(1:ii,1:jj, 1:kk, index)+thbase+1000.)
        call fabm_model%link_interior_data(standard_variables%density,codens(1:ii,1:jj, 1:kk))
        call fabm_model%link_interior_data(standard_variables%pressure,codepth(1:ii,1:jj, 1:kk))

        call update_fabm_state(index)
    end subroutine update_fabm_data

    subroutine update_fabm_state(index)
        integer, intent(in) :: index

        if (current_time_index == index) return

        ! Update mask and index of bottommost layer
        call fabm_set_mask(fabm_model, mask(:, :, :, index), mask(:, :, 1, index))
        call fabm_model%set_bottom_index(kbottom(:, :, index))

        ! Send pointers to state variable data to FABM
        call fabm_model%link_all_interior_state_data(tracer(1:ii, 1:jj, 1:kk, index, :))
        call fabm_model%link_all_bottom_state_data(fabm_bottom_state(1:ii, 1:jj, index, :))
        call fabm_model%link_all_surface_state_data(fabm_surface_state(1:ii, 1:jj, index, :))

        current_time_index = index
    end subroutine update_fabm_state

    function hycom_fabm_allocate_mean_output(idm, jdm, kdm) result(n)
      integer, intent(in) :: idm, jdm, kdm
      integer :: n

      type (type_interior_output),   pointer :: interior_output
      type (type_horizontal_output), pointer :: horizontal_output

      n = 0
      interior_output => first_interior_output
      do while (associated(interior_output))
        allocate(interior_output%mean(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm))
        n = n + (idm+2*nbdy)*(jdm+2*nbdy)*kdm
        interior_output => interior_output%next
      end do
      horizontal_output => first_horizontal_output
      do while (associated(horizontal_output))
        allocate(horizontal_output%mean(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
        n = n + (idm+2*nbdy)*(jdm+2*nbdy)
        horizontal_output => horizontal_output%next
      end do
    end function hycom_fabm_allocate_mean_output

    subroutine hycom_fabm_zero_mean_output()
      type (type_interior_output),   pointer :: interior_output
      type (type_horizontal_output), pointer :: horizontal_output

      interior_output => first_interior_output
      do while (associated(interior_output))
        interior_output%mean = 0
        interior_output => interior_output%next
      end do
      horizontal_output => first_horizontal_output
      do while (associated(horizontal_output))
        horizontal_output%mean = 0
        horizontal_output => horizontal_output%next
      end do
    end subroutine hycom_fabm_zero_mean_output

    subroutine hycom_fabm_increment_mean_output(n, s)
      integer, intent(in) :: n
      real, intent(in) :: s

      type (type_interior_output),   pointer :: interior_output
      type (type_horizontal_output), pointer :: horizontal_output
      real, pointer :: pdata2d(:,:), pdata3d(:,:,:)
      integer :: i, j, k

      interior_output => first_interior_output
      do while (associated(interior_output))
        if (associated(interior_output%data4d)) then
          pdata3d => interior_output%data4d(1:ii, 1:jj, 1:kk, n)
        else
          pdata3d => interior_output%data3d(1:ii, 1:jj, 1:kk)
        end if
        do k=1,kk
          do j=1,jj
            do i=1,ii
              if (associated(interior_output%data3d) .and. pdata3d(i, j, k) .lt. -1E18) pdata3d(i, j, k) = 0.0 ! the intention here is to remove the mask (-2E20) from fabm 
!              if (i==itest .and. j==jtest .and. associated(interior_output%data3d)) write(*,*)'HERE',interior_output%mean(i, j, k),k, pdata3d(i, j, k)
              if (SEA_P) interior_output%mean(i, j, k) = interior_output%mean(i, j, k) + s * dp(i, j, k, n) * pdata3d(i, j, k)
            end do
          end do
        end do
        interior_output => interior_output%next
      end do

      horizontal_output => first_horizontal_output
      do while (associated(horizontal_output))
        if (associated(horizontal_output%data3d)) then
          pdata2d => horizontal_output%data3d(1:ii, 1:jj, n)
        else
          pdata2d => horizontal_output%data2d(1:ii, 1:jj)
        end if
        do j=1,jj
          do i=1,ii
            if (SEA_P) horizontal_output%mean(i, j) = horizontal_output%mean(i, j) + s * pdata2d(i,j)
          end do
        end do
        horizontal_output => horizontal_output%next
      end do
    end subroutine hycom_fabm_increment_mean_output

    subroutine hycom_fabm_end_mean_output(q, dp_m, dpthin)
      real, intent(in) :: q
      real, intent(in) :: dp_m(:, :, :)
      real, intent(in) :: dpthin

      type (type_interior_output),   pointer :: interior_output
      type (type_horizontal_output), pointer :: horizontal_output
      integer :: i, j, k

      interior_output => first_interior_output
      do while (associated(interior_output))
        do k=1,kk
          do j=1,jj
            do i=1,ii
              if (SEA_P) then 
                if (dp_m(i, j, k) .ge. dpthin) then
                  interior_output%mean(i, j, k) = interior_output%mean(i, j, k)*q/dp_m(i, j, k)
                else
                  interior_output%mean(i, j, k) = interior_output%mean(i, j, k-1)
                end if
              end if
            end do
          end do
        end do
        interior_output => interior_output%next
      end do

      horizontal_output => first_horizontal_output
      do while (associated(horizontal_output))
        do j=1,jj
          do i=1,ii
            if (SEA_P) horizontal_output%mean(i, j) = horizontal_output%mean(i, j)*q
          end do
        end do
        horizontal_output => horizontal_output%next
      end do
    end subroutine hycom_fabm_end_mean_output

    subroutine hycom_fabm_write_mean_output(nop, nopa, nmean, time_ave)
      use mod_za  ! HYCOM I/O interface

      integer, intent(in) :: nop, nopa, nmean
      real(8), intent(in) :: time_ave

      real :: xmin, xmax, coord
      integer :: k

      type (type_interior_output),   pointer :: interior_output
      type (type_horizontal_output), pointer :: horizontal_output

      do k=1,kk
        coord = sigma(k)
        interior_output => first_interior_output
        do while (associated(interior_output))
          call zaiowr(interior_output%mean(1-nbdy,1-nbdy,k),ip,.true.,xmin,xmax, nopa, .false.)
          if     (mnproc.eq.1) then
            write (nop,117) interior_output%metadata%name(1:8),nmean,time_ave,k,coord,xmin,xmax
            call flush(nop)
          endif !1st tile
          interior_output => interior_output%next
        end do
      end do

      coord = 0.
      horizontal_output => first_horizontal_output
      do while (associated(horizontal_output))
        call zaiowr(horizontal_output%mean,ip,.true.,xmin,xmax, nopa, .false.)
        if     (mnproc.eq.1) then
          write (nop,117) horizontal_output%metadata%name(1:8),nmean,time_ave,0,coord,xmin,xmax
          call flush(nop)
        endif !1st tile
        horizontal_output => horizontal_output%next
      end do

      return
 117  format (a8,' =',i11,f11.3,i3,f7.3,1p2e16.7)
     end subroutine hycom_fabm_write_mean_output


    subroutine hycom_fabm_nest_next()
      integer :: i, j, k, nestn
    if (nested_bio) then
      do nestn = 1,nested_number
        do k= 1,kk
          do j= 1,jj
            do i= 1,ii
!              if (i==itest.and.j==jtest.and.k==1) then
!                write(*,*)'HEREnext',nestn,nested_number,nested_data_dev(i,j,k,1,nestn),nested_data_dev(i,j,k,2,nestn)
!              endif
              nested_data_dev(i,j,k,1,nestn) = nested_data_dev(i,j,k,2,nestn)
            enddo
          enddo
        enddo
      enddo
    end if ! nested_bio
    end subroutine hycom_fabm_nest_next

    subroutine hycom_fabm_nest_read(iyear,iday,ihour,lslot)
      use mod_za
      integer, intent(in) :: iyear,iday,ihour,lslot

      integer, parameter :: iunit = 920
      character(len=27) :: flnm
      integer :: iline
      character(len=80) :: cline
      character(len=6) :: cvarin
      integer :: ios,idmtst,jdmtst
      !integer :: nestn,i,j,k
    if (nested_bio) then
      write(flnm,'("nest/archv_fabm.",i4.4,"_",i3.3,"_",i2.2)') iyear, iday, ihour

      if (mnproc.eq.1) write (lp,*) 'hycom_fabm_nest_rdnest_in: ', flnm

      call xcsync(flush_lp)

      call zaiopf(flnm//'.a','old', iunit)
      if (mnproc.eq.1) then  ! .b file from 1st tile only
        open (unit=uoff+iunit, file=flnm//'.b', form='formatted', status='old', action='read')
        do iline=1,7
          read(uoff+iunit,'(a)') cline
        end do
      end if !1st tile

      call zagetc(cline, ios, uoff+iunit)
      read(cline,*) idmtst, cvarin
      if (cvarin.ne.'idm   ') then
        if (mnproc.eq.1) then
          write(lp,*)
          write(lp,*) 'error in hycom_fabm_nest_rdnest_in - input ',cvarin,' but should be idm   '
          write(lp,*)
        end if !1st tile
        call xcstop('(hycom_fabm_nest_rdnest_in)')
               stop '(hycom_fabm_nest_rdnest_in)'
      end if
      call zagetc(cline, ios, uoff+iunit)
      read(cline,*) jdmtst, cvarin
      if (cvarin.ne.'jdm   ') then
        if (mnproc.eq.1) then
          write(lp,*)
          write(lp,*) 'error in hycom_fabm_nest_rdnest_in - input ',cvarin,' but should be jdm   '
          write(lp,*)
        end if !1st tile
        call xcstop('(hycom_fabm_nest_rdnest_in)')
               stop '(hycom_fabm_nest_rdnest_in)'
      end if

      if (idmtst.ne.itdm .or. jdmtst.ne.jtdm) then
        if (mnproc.eq.1) then
          write(lp,*)
          write(lp,*) 'error in hycom_fabm_nest_rdnest_in - input idm,jdm not consistent with parameters'
          write(lp,*) 'idm,jdm = ',itdm,  jtdm,  '  (dimensions.h)'
          write(lp,*) 'idm,jdm = ',idmtst,jdmtst,'  (input)'
          write(lp,*)
        end if !1st tile
        call xcstop('(hycom_fabm_nest_rdnest_in)')
               stop '(hycom_fabm_nest_rdnest_in)'
      end if

      if (mnproc.eq.1) then  ! .b file from 1st tile only
        read (uoff+iunit,*)
      end if

      do while (read_next_field())
      end do

      if     (mnproc.eq.1) then  ! .b file from 1st tile only
      close( unit=uoff+iunit)
      endif
      call zaiocl(iunit)
    end if ! nested_bio
    contains

      function read_next_field() result(success)
        logical :: success

        integer :: ios,k,nnstep,i, nestn
        real :: hmina,hminb,hmaxa,hmaxb,timein,thet

        success = .true.
        call zagetc(cline, ios, uoff+iunit)
        if (ios < 0) then
          ! End of file reached
          success = .false.
          return
        elseif (ios > 0) then
          if (mnproc.eq.1) then
            write(lp,*)
            write(lp,*) 'error in mod_hycom_fabm::read_next_field - error reading next field'
            write(lp,*) 'iunit,ios = ',iunit,ios
            write(lp,*)
          end if !1st tile
          call xcstop('(rd_archive)')
                 stop '(rd_archive)'
        end if

        ! Look up FABM variable with the name found in the nesting input.
        do nestn=1,nested_number
          if (cline(1:8) == nested_variables(nestn)) exit
        enddo
        if (nestn > nested_number) then
          call zaiosk(iunit)
          return
        end if
!          write(*,'(A8,1x,I1,1x,I1,1x,A8)')'HEREread',nestn,nested_number,nested_variables(nestn)
          i = index(cline,'=')
          read(cline(i+1:),*) nnstep,timein,k,thet,hminb,hmaxb

          if (hminb.eq.hmaxb) then  !constant field
            nested_data_dev(:,:,k,lslot,nestn) = hminb
            call zaiosk(iunit)
          else
            call zaiord(nested_data_dev(:,:,k,lslot,nestn),ip,.false.,hmina,hmaxa,iunit)
            if (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4) then
              if (mnproc.eq.1) write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')'error - .a and .b files not consistent:', &
                  '.a,.b min = ',hmina,hminb,hmina-hminb,'.a,.b max =',hmaxa,hmaxb,hmaxa-hmaxb
              ! We could have stopped here, but that is commented out in
              ! forfun.F/rd_archive
            end if
          end if
      end function

    end subroutine hycom_fabm_nest_read

    subroutine hycom_fabm_nest_update(n)
      integer, intent(in) :: n
      integer :: i, j, k, nestn

    if (nested_bio) then
      do nestn=1,nested_number
        do k=1,kk
          do j=1,jj
            do i=1,ii
!              if (i==itest.and.j==jtest) write(*,*)'HEREupd',k,nestn,(delt1*rmunp(i,j)* &
!                (nested_data_dev(i,j,k,ln0,nestn)*wn0 +nested_data_dev(i,j,k,ln1,nestn)*wn1) )/(1.0 + delt1*rmunp(i,j))
              tracer(i,j,k,n,istate_dev(nestn))=(tracer(i,j,k,n,istate_dev(nestn)) + delt1*rmunp(i,j)* &
                (nested_data_dev(i,j,k,ln0,nestn)*wn0 +nested_data_dev(i,j,k,ln1,nestn)*wn1) )/(1.0 + delt1*rmunp(i,j))
            end do
          end do
        end do
      enddo
    end if ! nested_bio
    end subroutine hycom_fabm_nest_update

    subroutine fabm_gettime()!(modeltime,yrflag)

        integer :: i,dayy,iyr, nleap,days_in_year,days_in_month(12)
        real*8  :: dtim1, day

        iyr   = (modeltime-1.d0)/365.25d0
        nleap = iyr/4
        dtim1 = 365.d0*iyr + nleap + 1.d0
        day   = modeltime - dtim1 + 1.d0
        if     (dtim1.gt.modeltime) then
          iyr = iyr - 1
        elseif (day.ge.367.d0) then
          iyr = iyr + 1
        elseif (day.ge.366.d0 .and. mod(iyr,4).ne.3) then
          iyr = iyr + 1
        endif
        nleap = iyr/4
        dtim1 = 365.d0*iyr + nleap + 1.d0

        modelyear =  1901 + iyr
        modelday  =  modeltime - dtim1 + 1.001d0
        !ihour = (dtime - dtim1 + 1.001d0 - iday)*24.d0
        days_in_year =daysinyear  (modelyear,yrflag)
        days_in_month=monthsofyear(modelyear,yrflag)
        dayy = int(modelday-1)
        do i=1,12
           if (dayy >= days_in_month(i)) then
              dayy=dayy-days_in_month(i)
           else
              modelmonth=i
              exit
           endif
        enddo
!        write(*,*)'GETTIME',modelyear, modelday, modelmonth

    end subroutine fabm_gettime

    integer function daysinyear(year,yrflag)
        implicit none
        integer,intent(in) :: year,yrflag
        if (yrflag==0) then
           daysinyear=360
        elseif (yrflag==1) then
           daysinyear=365
        elseif (yrflag==2) then
           daysinyear=366
        elseif (yrflag==3) then
           if (mod(year,4)/=0 .or. (mod(year,100)==0.and.mod(year,400)/=0)) then
              daysinyear=365
           else
              daysinyear=366
           end if
        end if
    end function daysinyear

function monthsofyear(year,yrflag)
        implicit none
        integer, dimension(12),parameter :: months_standard = &
           (/31,28,31,30,31,30,31,31,30,31,30,31/)
        integer, dimension(12),parameter :: months_leapyear = &
           (/31,29,31,30,31,30,31,31,30,31,30,31/)
        integer, dimension(12),parameter :: months_360 = &
           (/30,30,30,30,30,30,30,30,30,30,30,30/)
        integer, dimension(12),parameter :: months_365 = &
           (/31,28,31,30,31,30,31,31,30,31,30,31/)
        integer, dimension(12),parameter :: months_366 = &
           (/31,29,31,30,31,30,31,31,30,31,30,31/)
        integer :: monthsofyear(12)
        integer, intent(in) :: year,yrflag
        if (yrflag==0) then
           monthsofyear=months_360
        elseif (yrflag==1) then
           monthsofyear=months_365
        elseif (yrflag==2) then
           monthsofyear=months_366
        elseif (yrflag==3) then
           if (daysinyear(year,3)==366) then
              monthsofyear=months_leapyear
           else
              monthsofyear=months_standard
           end if
        end if
end function
#endif
end module

module mod_NERSCnml
! allow for namelist inputs in order to 
  use mod_xc
  use mod_za  ! HYCOM I/O interface
  use mod_cb_arrays, only: priver
  implicit none
  private
 
  real, save, public   :: &
    sssrmx_scalar       ! Limit for salinity relax Only used if

  logical,save, public :: &
    write_arche, &      ! print arche files or not
    sss_underice, &     ! relaxtion under ice (false if no) 
    highfq_river        ! high frequency rivers (false if no) 

  public NERSC_init

  contains

  !----------------------------------------------------------
  subroutine NERSC_init
  !----------------------------------------------------------
    implicit none
    integer (4), parameter :: funi=503
    integer (4) :: nml_err
    namelist /hycom_nml/ write_arche,sss_underice,highfq_river,sssrmx_scalar
    !default values
    write_arche     = .false.
    sss_underice    = .false.
    highfq_river    = .false.
    sssrmx_scalar   = 99.

    ! read namelist
    open (funi, file='./hycom_opt', status='old',iostat=nml_err)
    if (nml_err .ne. 0) then
      if  (mnproc.eq.1) then
        write (lp,'(a)') &
          'NERSC HYCOM ERROR: hycom_nml. Did not read: ./hycom_opt'
        call flush(lp)
      endif
      call xcstop('(NERSC_nml)')
    endif
    do while (nml_err == 0)
      read(funi, nml=hycom_nml,iostat=nml_err)
      if (nml_err > 0) then
        if (mnproc.eq.1) then
          write (lp,'(a)') &
          'NERSC HYCOM ERROR: Can not read namelist: hycom_opt'
          call flush(lp)
        endif
        call xcstop('(NERSC_nml)')
      endif
    end do
    close(funi)
    if (mnproc.eq.1) then
      write (lp,*)'NERSC HYCOM: Reading hycom_nml from: ./hycom_opt'
      write (lp,*)'NERSC HYCOM: Write arche    = ',write_arche
      write (lp,*)'NERSC HYCOM: sss_underice   = ',sss_underice
      write (lp,*)'NERSC HYCOM: sssrmx_scalar  = ',sssrmx_scalar
      write (lp,*)'NERSC HYCOM: highfq_river   = ',highfq_river
    endif !1st tile
    call flush(lp)
    if (priver .and. highfq_river) then
        if (mnproc.eq.1) then
           write(lp,*) &
           'error - priver must be .false. for highfq_river=.true.'
           call flush(lp)
        endif !1st tile
        call xcstop('(NERSC_nml)')
        stop '(NERSC_nml)'
    endif
    call xcsync(flush_lp)
  end subroutine NERSC_init

  !---------------------------------------------------------------------+

end module mod_NERSCnml

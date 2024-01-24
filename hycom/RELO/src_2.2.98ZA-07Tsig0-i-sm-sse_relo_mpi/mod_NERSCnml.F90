module mod_NERSCnml
! allow for namelist inputs in order to 
  use mod_xc
  use mod_za  ! HYCOM I/O interface

  implicit none
  private
 

  logical,save, public :: &
    write_arche, &      ! print arche files or not
    sss_underice     ! relaxtion under ice (false if no) 

  public NERSC_init

  contains

  !----------------------------------------------------------
  subroutine NERSC_init
  !----------------------------------------------------------
    implicit none
    integer (4), parameter :: funi=503
    integer (4) :: nml_err
    namelist /hycom_nml/ write_arche,sss_underice
    !default values
    write_arche     = .false.
    sss_underice    = .false.
    ! read namelist
    open (funi, file='../hycom_opt', status='old',iostat=nml_err)
    if (nml_err .ne. 0) then
      if  (mnproc.eq.1) then
        write (lp,'(a)') &
          'NERSC HYCOM ERROR: WARNING: hycom_nml namelist not read from file: ../hycom_opt'
        call flush(lp)
      endif
    endif
    do while (nml_err == 0)
      read(funi, nml=hycom_nml,iostat=nml_err)
      if (nml_err > 0) then
        if (mnproc.eq.1) then
          write (lp,'(a)') &
          'NERSC HYCOM ERROR: Can not read namelist: hycom_opt'
          call flush(lp)
          call xcstop('(NERSC_nml)')
          stop '(NERSC_nml)'
        endif
      endif
    end do
    close(funi)
    if (mnproc.eq.1) then
      write (lp,*)'NERSC HYCOM: Reading hycom_nml from: ../hycom_opt'
      write (lp,*)'NERSC HYCOM: Write arche    = ',write_arche
      write (lp,*)'      HYCOM: sss_underice   = ',sss_underice
    endif !1st tile
    call xcsync(flush_lp)

  end subroutine NERSC_init

  !---------------------------------------------------------------------+

end module mod_NERSCnml

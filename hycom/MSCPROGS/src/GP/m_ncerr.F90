module m_ncerr
contains
   subroutine ncerr(errcode)
   use netcdf
   implicit none
   integer, intent(in) :: errcode
   if (errcode/=NF90_NOERR) then
      write(6,'(a)') NF90_STRERROR(errcode)
      stop '(ncerr)'
   end if
   end subroutine
end module m_ncerr

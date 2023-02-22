










module m_handle_err
contains

   subroutine handle_err(errcode)
      use netcdf
      implicit none
      integer, intent(in) :: errcode

      if (errcode/=NF90_NOERR) then
         write(6,'(a)') NF90_STRERROR(errcode)
         stop '(handle_err)'
      end if

   end subroutine



end module m_handle_err


module m_handle_err

contains


   subroutine handle_err(errcode,halt)
      use netcdf
      implicit none
      integer, intent(in) :: errcode
      logical, optional, intent(in) :: halt

      if (errcode/=NF90_NOERR) then
         write(6,'(a)') NF90_STRERROR(errcode)
      end if

      if (present(halt) .and. halt .and. errcode/=NF90_NOERR )  then
         print *, '(handle_err)'
         call exit(1)
      end if
   end subroutine handle_err

   subroutine stopon_err(errcode)
      implicit none
      integer,intent(in) :: errcode
      call handle_err(errcode,.true.)
   end subroutine stopon_err




end module m_handle_err














      
         

      


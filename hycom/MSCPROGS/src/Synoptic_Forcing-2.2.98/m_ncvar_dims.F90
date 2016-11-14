module m_ncvar_dims

      integer, parameter, private :: lprt=6

      private :: nf90_handle_err

contains

   subroutine ncvar_dims(filename,varname,dimsizes,ndims,recdim)
      use netcdf
      implicit none

      character(len=*), intent(in)  :: varname,filename
      integer         , intent(out) :: ndims
      integer         , intent(out) :: recdim
      integer         , intent(out) :: dimsizes(NF90_MAX_VAR_DIMS)

      integer :: i
      integer :: dimids(NF90_MAX_VAR_DIMS)
      integer :: rdim,varid,ncid


      ! Open NetCDF file
      call nf90_handle_err( &
           NF90_OPEN(filename,NF90_NOCLOBBER,ncid),&
           'Can not find '//filename)

      ! Inquire on variable -- get its id
      call nf90_handle_err( &
           nf90_inq_varid(ncid, varname, varid),   &
           filename//' does not contain variable '//varname)

      ! Get unlimited (record) dimension
      call nf90_handle_err(&
           nf90_Inquire(ncid, unlimitedDimId=rdim),&
           'Error On inquiring on record dimension ')

      !print *,'rdim=',rdim

      ! Inquire on variable -- dimensions and dimension ids
      call nf90_handle_err( &
           nf90_Inquire_Variable(ncid, varid,ndims=ndims,dimids=dimids),&
           'Error on inquiring variable '//varname)

      ! Inquire on dimensions -- sizes
      dimsizes=0
      recdim=0
      do i=1,ndims
         call nf90_handle_err( &
              nf90_Inquire_Dimension(ncid, dimids(i), len=dimsizes(i)), &
              ' Error on inquiring dimension for var '//varname)
         if (dimids(i)==rdim) recdim=i
      end do

      !print *,'dimids : ',recdim,dimids

      ! Close NetCDF file
      call nf90_handle_err(NF90_CLOSE(ncid),'Error closing file '//filename)
      !write(lprt,*)

   end subroutine ncvar_dims


   subroutine nf90_handle_err(errcode,info)
      use netcdf
      implicit none
      integer, intent(in) :: errcode
      character(len=*), intent(in) :: info

      if (errcode/=NF90_NOERR) then
         write(lprt,'(a)') NF90_STRERROR(errcode)
         write(lprt,'(a)') info
         stop '(ncvar_dims)'
      end if
   end subroutine nf90_handle_err

end module m_ncvar_dims














      
         

      


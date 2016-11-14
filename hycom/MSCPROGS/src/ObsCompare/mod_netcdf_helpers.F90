module mod_netcdf_helpers
contains

   subroutine ncerr(error)
      use netcdf
      implicit none
      integer, intent(in) :: error
      if (error/=nf90_noerr) then
         print *,'ncerr: '//nf90_strerror(error)
         stop
      end if
   end subroutine ncerr

   subroutine getdimms(ncid,dname,dmid,dmlen)
      use netcdf
      implicit none
      integer,          intent(in)  :: ncid
      character(len=*), intent(in)  :: dname
      integer,          intent(out) :: dmid,dmlen
      call ncerr( nf90_inq_dimid(ncid,dname,dmid) )
      call ncerr( nf90_inquire_dimension(ncid,dmid,len=dmlen) )
   end subroutine

   logical function probe_var(ncid,varname)
      use netcdf
      implicit none
      character(len=*), intent(in) :: varname
      integer         , intent(in) :: ncid
      integer :: test, varid
      test=nf90_inq_varid(ncid,varname,varid) 
      probe_var = test==nf90_noerr
   end function probe_var

   subroutine getvarstat(ncid,varname,varid,var_ndim,var_dimids)
      use netcdf
      implicit none
      integer,                          intent(in)  :: ncid
      character(len=*),                 intent(in)  :: varname
      integer,                          intent(out) :: varid,var_ndim
      integer, dimension(nf90_max_dims),intent(out) :: var_dimids
      call ncerr( nf90_inq_varid(ncid,varname,varid) )
      call ncerr( nf90_inquire_variable(ncid,varid,ndims=var_ndim,dimids=var_dimids) )
   end subroutine getvarstat

end module mod_netcdf_helpers

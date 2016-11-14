module m_dBM_MLD_read

contains

   subroutine dBM_MLD_read(mld,mldlon,mldlat,nblon,nblat,undef,jday)
   use netcdf
   use mod_netcdf_helpers
   implicit none

   integer, intent(out) :: nblon, nblat
   integer, intent(in ) :: jday
   real, pointer, dimension(:)   :: mldlon, mldlat
   real, pointer, dimension(:,:) :: mld
   real, intent(in) :: undef

   character(len=200) :: MLD_PATH
   character(len=*), parameter :: mldfile='mldT02_sk.nc'
   integer, dimension(nf90_max_dims) :: vdims ! Variable dimensions 
   integer                           :: vndim ! number of variable dims 
   integer :: ncid, dimid, lonvarid, latvarid, mldvarid
   integer :: stat


   call getenv('MLD_PATH',MLD_PATH)
   if (trim(MLD_PATH)=="") then
      print *,'Supply path to de Boyer Montegut climatology file ('//trim(mldfile)//')'
      print *,'Specify this in environment variable MLD_PATH'
      stop
   end if
   MLD_PATH=trim(MLD_PATH)//'/'

   print *,'reading MLD fields from '//trim(MLD_PATH)//trim(mldfile)
   stat = nf90_open(trim(MLD_PATH)//trim(mldfile),nf90_nowrite,ncid) 
   if (stat /= nf90_noerr) then
      call ncerr( stat)
   end if

   call getdimms(ncid,'lon',dimid,nblon)
   call getdimms(ncid,'lat',dimid,nblat)

   allocate(mldlon(nblon))
   allocate(mldlat(nblat))
   allocate(mld   (nblon,nblat)) 
   
   call getvarstat(ncid,'lon',lonvarid ,vndim  ,vdims)
   call getvarstat(ncid,'lat',latvarid ,vndim  ,vdims)
   call getvarstat(ncid,'mld',mldvarid ,vndim  ,vdims)

   
   call ncerr(nf90_get_var(ncid, lonvarid,mldlon))
   call ncerr(nf90_get_var(ncid, latvarid,mldlat))
   call ncerr(nf90_get_var(ncid, mldvarid,mld, &
                           start=(/1,1,1/),    &
                           count=(/nblon,nblat,1/)))   ! For testing - time interpolation later

   where(mld<0.)  mld=undef
   where(mld>1e8) mld=undef

   stat = nf90_close(ncid)
   if (stat /= nf90_noerr) then
      call ncerr( stat)
   end if

   end subroutine

end module


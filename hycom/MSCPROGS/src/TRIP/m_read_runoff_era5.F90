module m_read_runoff_era5

   character(len= 200), private :: ropath=''
   character(len=*), parameter, private :: lsmfile='6h.RO_1992.nc'
   real, dimension(:), allocatable :: lon, lat
   integer :: nlon=1440, nlat=721
   private :: setpath
contains


   subroutine init_runoff_era5()
   use netcdf
   use m_handle_err
   implicit none
   integer :: varid, dimsize, ndims, natts, ncid, ierr
   integer :: dimids(NF90_MAX_VAR_DIMS)
   integer :: dimsizes(NF90_MAX_VAR_DIMS),j
   real :: dlon,dlat,flon,flat

   ! Path of ERA-5 data
    call setpath_era5()

! Input data grid (for now - simple lat-lon (possibly gaussian) grid )
    ierr= nf90_open(trim(ropath)//lsmfile,NF90_NOCLOBBER,ncid)
    if (ierr/=NF90_NOERR) then
      print *,'Could not open '//trim(ropath)//lsmfile
      stop
   end if
     ! Open file
   call handle_err(nf90_open(trim(ropath)//lsmfile,NF90_NOCLOBBER,ncid))

   ! Get longitude
   call handle_err( nf90_inq_varid(ncid, 'longitude', varid))
   call handle_err( &
        nf90_Inquire_Variable(ncid,varid,ndims=ndims,dimids=dimids,natts=natts))
   call handle_err( &
        nf90_Inquire_Dimension(ncid, dimids(1), len=dimsize))
   nlon=dimsize 
   allocate(lon(nlon))
   call handle_err(nf90_Get_Var(ncid, varid,lon))
   print *,'min max lon',lon(1),lon(nlon),nlon

   ! Get latitude                                                                                           
   call handle_err( nf90_inq_varid(ncid, 'latitude', varid))
   call handle_err( &
        nf90_Inquire_Variable(ncid,varid,ndims=ndims,dimids=dimids,natts=natts))
   call handle_err( &
        nf90_Inquire_Dimension(ncid, dimids(1), len=dimsize))
   nlat=dimsize
   allocate(lat(nlat))
<<<<<<< HEAD
   flat=90
   dlat=-0.25
   do j= 1,nlat
      lat(j)=flat+(j-1)*dlat
   enddo
=======
   call handle_err(nf90_Get_Var(ncid, varid,lat))
>>>>>>> 2f21a4dd36298dc0c24d7bc87f73fc528e36c409
   print *,'min max lat',lat(1),lat(nlat),nlat
   call handle_err(nf90_Close(ncid))

   end subroutine 


  subroutine setpath_era5()
  implicit none
  character(len=200) :: cenv
  ! Path of ERA-5 data
  if (trim(ropath)=='') then
     call getenv('ERA5_PATH',cenv)
     if (trim(cenv)=='') then
        print *,'Could not get ERA5 path, make sure environment variable'
        print *,'ERA5_PATH points to location of ERA5 data'
        stop
     else
        ropath=trim(cenv)//'/'
     end if
  end if
  end subroutine

   subroutine read_runoff_era5(startyear,rtime,fld,nx,ny)
   use netcdf
   use mod_year_info
   use m_handle_err
   implicit none

   integer, intent(in)  :: startyear
   real   , intent(in)  :: rtime
   integer, intent(in ) :: nx,ny
   real   , intent(out) :: fld(nx,ny)

   ! For now
   character(len= *), parameter :: varname= 'RO'
   character(len=280):: rofile
   character(len= 4):: cyy


   integer :: ncid, varid, natts, ndims, irec
   integer :: dimx, dimy, dimtime
   integer :: thisyear,thismonth,thisday, thishour, jday

   integer :: dimids  (NF90_MAX_VAR_DIMS)
   integer :: dimsizes(NF90_MAX_VAR_DIMS)

   real :: offset, scale_factor

   ! Path of ERA-5 data
   call setpath_era5()

   ! Convert rtime to year record
   call juliantodate(floor(rtime),thisyear,thismonth,thisday,startyear,1,1)
   thishour=24*(rtime-floor(rtime))
   jday=datetojulian(thisyear,thismonth,thisday,thisyear,1,1)
   irec=(jday)*4 + thishour/6
   !print *,thisyear,thismonth,thisday, thishour,jday,irec

   !NB
   if (irec==0.and.thisday==1.and.thismonth==1.and.thishour==0) then
      thisyear=thisyear-1
      irec=1
   end if

   write(cyy,'(i4.4)') thisyear
   rofile=trim(ropath)//'6h.RO_'//cyy//'.nc'   
<<<<<<< HEAD
  
 ! Open file
=======
   print *,rofile

   ! Open file
>>>>>>> 2f21a4dd36298dc0c24d7bc87f73fc528e36c409
   call handle_err(nf90_open(trim(rofile),NF90_NOCLOBBER,ncid))

   ! Get variable id
   call handle_err( nf90_inq_varid(ncid, trim(varname), varid))

   ! Get dimensions of variable
   call handle_err( &
        nf90_Inquire_Variable(ncid,varid,ndims=ndims,dimids=dimids,natts=natts))

   ! Get the offset and scale factor
   call handle_err(&
        nf90_get_att(ncid, varid,'add_offset', offset))
   call handle_err(&
        nf90_get_att(ncid, varid,'scale_factor', scale_factor))
   print*, offset,scale_factor

   if (ndims/=3) then
      print *,'Expected 3D var '
      stop
   end if

   ! Inquire three var dims
   call handle_err( &
        nf90_Inquire_Dimension(ncid, dimids(1), len=dimx))
   call handle_err( &
        nf90_Inquire_Dimension(ncid, dimids(2), len=dimy))
   call handle_err( &
        nf90_Inquire_Dimension(ncid, dimids(3), len=dimtime))

   print *,dimx, dimy, dimtime
   if (dimx/=nx .or. dimy/=ny) then
      print *,'Error in dimensions'
      print *,dimx,nx
      print *,dimy,ny
      stop '(read_runoff)'
   end if


   ! Read variable
   call handle_err(NF90_GET_VAR(ncid,varid,fld,start=(/1,1,irec/)))
   call handle_err(NF90_close(ncid))

   ! Conversion from m to m/s
   fld = (fld*scale_factor + offset) / (6*3600.) 

   end subroutine read_runoff_era5






end module m_read_runoff_era5

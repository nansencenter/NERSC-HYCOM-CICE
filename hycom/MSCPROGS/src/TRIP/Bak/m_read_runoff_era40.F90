module m_read_runoff_era40

   character(len= 200), private :: ropath=''
   character(len=*), parameter, private :: lsmfile='ans.LSM.nc'
   real, dimension(:), allocatable :: lon, lat
   integer :: nlon, nlat
   private :: setpath
contains


   subroutine init_runoff_era40()
   use netcdf
   use m_handle_err
   implicit none
   integer :: varid, dimsize, ndims, natts, ncid, ierr
   integer :: dimids(NF90_MAX_VAR_DIMS)
   integer :: dimsizes(NF90_MAX_VAR_DIMS)

   ! Path of ERA-40 data
   call setpath_era40()

   ! Input data grid (for now - simple lat-lon (possibly gaussian) grid )
   ierr= nf90_open(trim(ropath)//lsmfile,NF90_NOCLOBBER,ncid)
   if (ierr/=NF90_NOERR) then
      print *,'Could not open '//trim(ropath)//lsmfile
      stop
   end if

   ! Get lon dim
   call handle_err( nf90_inq_varid(ncid, 'lon', varid))
   call handle_err( &
        nf90_Inquire_Variable(ncid,varid,ndims=ndims,dimids=dimids,natts=natts))
   call handle_err( &
        nf90_Inquire_Dimension(ncid, dimids(1), len=dimsize))
   nlon=dimsize 
   allocate(lon(nlon))
   call handle_err(nf90_Get_Var(ncid, varid,lon))

   ! Get lat dim
   call handle_err( nf90_inq_varid(ncid, 'lat', varid))
   call handle_err( &
        nf90_Inquire_Variable(ncid,varid,ndims=ndims,dimids=dimids,natts=natts))
   call handle_err( &
        nf90_Inquire_Dimension(ncid, dimids(1), len=dimsize))
   nlat=dimsize
   allocate(lat(nlat))
   call handle_err(nf90_Get_Var(ncid, varid,lat))

   call handle_err(nf90_Close(ncid))
   end subroutine 



  subroutine setpath_era40()
  implicit none
  character(len=200) :: cenv
  ! Path of ERA-40 data
  if (trim(ropath)=='') then
     call getenv('ERA40_PATH',cenv)
     if (trim(cenv)=='') then
        print *,'Could not get ERA40 path, make sure environment variable'
        print *,'ERA40_PATH points to location of ERA40 data'
        stop
     else
        ropath=trim(cenv)//'/'
     end if
  end if
  end subroutine

   subroutine read_runoff_era40(startyear,rtime,fld,nx,ny)
   use netcdf
   use mod_year_info
   use m_handle_err
   implicit none

   integer, intent(in)  :: startyear
   real   , intent(in)  :: rtime
   integer, intent(in ) :: nx,ny
   real   , intent(out) :: fld(nx,ny)

   ! For now
   character(len= *), parameter :: varname= 'RO_sfc'
   character(len=280):: rofile
   character(len= 4):: cyy


   integer :: ncid, varid, natts, ndims, irec
   integer :: dimx, dimy, dimtime
   integer :: thisyear,thismonth,thisday, thishour, jday

   integer :: dimids  (NF90_MAX_VAR_DIMS)
   integer :: dimsizes(NF90_MAX_VAR_DIMS)

   ! Path of ERA-40 data
   call setpath_era40()

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
   rofile = trim(ropath)//'fcs.6h.'//cyy//'.RO.nc'

   ! Open file
   call handle_err(nf90_open(trim(rofile),NF90_NOCLOBBER,ncid))

   ! Get variable id
   call handle_err( nf90_inq_varid(ncid, trim(varname), varid))

!   ! Get dimensions of variable
!   call handle_err( &
!        nf90_Inquire_Variable(ncid,varid,ndims=ndims,dimids=dimids,natts=natts))
!
!   if (ndims/=3) then
!      print *,'Expected 3D var '
!      stop
!   end if
!
!   ! Inquire three var dims
!   call handle_err( &
!        nf90_Inquire_Dimension(ncid, dimids(1), len=dimx))
!   call handle_err( &
!        nf90_Inquire_Dimension(ncid, dimids(2), len=dimy))
!   call handle_err( &
!        nf90_Inquire_Dimension(ncid, dimids(3), len=dimtime))
!
!   !print *,dimx, dimy, dimtime
!   if (dimx/=nx .or. dimy/=ny) then
!      print *,'Error in dimensions'
!      print *,dimx,nx
!      print *,dimy,ny
!      stop '(read_runoff)'
!   end if


   ! Read variable
   call handle_err(NF90_GET_VAR(ncid,varid,fld,start=(/1,1,irec/)))
   call handle_err(NF90_close(ncid))

   ! Conversion from m to m/s
   fld = fld / (6*3600.) 

   end subroutine read_runoff_era40






end module m_read_runoff_era40

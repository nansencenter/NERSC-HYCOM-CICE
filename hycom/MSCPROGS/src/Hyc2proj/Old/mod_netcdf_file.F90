module mod_netcdf_file
use netcdf

type file_state
   integer :: ncid=-1
   integer :: lon_dimension_id, lat_dimension_id, k_dimension_id, t_dimension_id
   integer :: vars_1d(1),vars_2d(2), vars_3d(3), vars_4d(4)
   character(len=20) ::  grid_mapping='', coords_mapping='', char_scale_fac='', &
                         cprojection=''
   character(len=200) :: ncfil=''
   logical :: stations=.false.
end type

!private :: handle_err
contains


   subroutine openNCFile(ncstate,ncfil)
   implicit none
   character(len=*), intent(in)  :: ncfil
   type(file_state), intent(out) :: ncstate
   ncstate%ncfil=ncfil
   if (ncstate%ncid/=-1) then
      print *,'Only one netcdf file open at a time'
      stop '(OpenNCFile)'
   end if
   if (NF90_CREATE(ncstate%ncfil,NF90_CLOBBER,ncstate%ncid) /= NF90_NOERR) then
      print *,'An error occured when opening the netcdf file'
      stop '(openNCFile)'
   end if
   end subroutine

   subroutine closeNCFile(ncstate)
   implicit none
   type(file_state), intent(inout) :: ncstate
   if (ncstate%ncid==-1) then
      print *,'Cant close a closed file..'
      stop '(CloseNCFile)'
   end if
   if (NF90_CLOSE(ncstate%ncid)/=NF90_NOERR) then
      print *,'An error occured when opening the netcdf file'
      stop '(closeNCFile)'
   end if
   ncstate%ncid=-1
   end subroutine

   subroutine DefNCDim(ncstate,cprojection,nlons,nlats,ndeep)
   implicit none
   type(file_state), intent(inout) :: ncstate
   character(len=*), intent(in) :: cprojection
   integer, intent(in) :: nlons, nlats, ndeep
   ncstate%cprojection=cprojection
   ! Define dimensions -- idm,jdm,kdm
   ! Following MERCATOR name definitions for netcdf var name
   if (trim(ncstate%cprojection)=='regular') then
      call handle_err(NF90_DEF_DIM(ncstate%ncid,'longitude',nlons,ncstate%lon_dimension_id))
      call handle_err(NF90_DEF_DIM(ncstate%ncid,'latitude' ,nlats,ncstate%lat_dimension_id))
   else if (trim(ncstate%cprojection)=='ps_amsr' .or. trim(ncstate%cprojection)=='polar_stereographic' .or. &
            trim(ncstate%cprojection)=='native') then
      call handle_err(NF90_DEF_DIM(ncstate%ncid,'x',nlons,ncstate%lon_dimension_id))
      call handle_err(NF90_DEF_DIM(ncstate%ncid,'y',nlats,ncstate%lat_dimension_id))
   else
      print *,'Unknown projection '//trim(ncstate%cprojection)
   end if
   call handle_err(NF90_DEF_DIM(ncstate%ncid,'depth'  ,ndeep,ncstate%k_dimension_id))
   call handle_err(NF90_DEF_DIM(ncstate%ncid,'time',NF90_Unlimited,ncstate%t_dimension_id))
   ! Define dimension "holders"
   ncstate%vars_4d=(/ncstate%lon_dimension_id,ncstate%lat_dimension_id,ncstate%k_dimension_id,ncstate%t_dimension_id/)
   ncstate%vars_3d=(/ncstate%lon_dimension_id,ncstate%lat_dimension_id,ncstate%t_dimension_id/)
   ncstate%vars_2d=(/ncstate%lon_dimension_id,ncstate%lat_dimension_id/)
   end subroutine




   subroutine DefNCDimStat(ncstate,nstat,ndeep)
   ! Define dimensions -- idm,jdm,kdm
   ! Following MERCATOR name definitions for netcdf var name
   type(file_state), intent(inout) :: ncstate
   integer, intent(in) :: nstat, ndeep
   ! "Cheat" - and use lon dimension id as station dim id
   call handle_err(NF90_DEF_DIM(ncstate%ncid,'station',nstat,ncstate%lon_dimension_id))
   call handle_err(NF90_DEF_DIM(ncstate%ncid,'depth'  ,ndeep,ncstate%k_dimension_id))
   call handle_err(NF90_DEF_DIM(ncstate%ncid,'time',1,ncstate%t_dimension_id))
   ncstate%vars_4d=(/ncstate%lon_dimension_id,ncstate%k_dimension_id,ncstate%t_dimension_id,-1/)
   ncstate%vars_3d=(/ncstate%lon_dimension_id,ncstate%t_dimension_id,-1/)
   ncstate%vars_2d=(/ncstate%lon_dimension_id,-1/)
   ncstate%stations=.true.
   end subroutine


   subroutine DefNCHeader(ncstate,hfile,cver)
   use mod_year_info
   use mod_hycomfile_io
   implicit none
   type(file_state) :: ncstate
   type(hycomfile)  :: hfile
   character(len=*), intent(in) :: cver

   type (year_info) :: rtd, rtb 
   character(len=8)   :: dateinfo

   call forecastDate(hfile,rtd)
   call startDate(hfile,rtb)

   call handle_err(NF90_PUT_ATT(ncstate%ncid,NF90_GLOBAL,'title', 'TOPAZ model results'))
   !call handle_err(NF90_PUT_ATT(ncstate%ncid,NF90_GLOBAL,'comment', trim(comment)))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,NF90_GLOBAL,'institution',  'NERSC, Thormoehlens gate 47, N-5006 Bergen, Norway'))
   call date_and_time(date=dateinfo)
   call handle_err(NF90_PUT_ATT(ncstate%ncid,NF90_GLOBAL,'history', dateinfo(1:8)//':Created by program hyc2proj, version '//cver))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,NF90_GLOBAL,'source', 'NERSC-HYCOM model fields'))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,NF90_GLOBAL,'references', 'http://topaz.nersc.no'))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,NF90_GLOBAL,'field_type',  'Files based on file type '//trim(hfile%ftype)))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,NF90_GLOBAL,'Conventions',  'CF-1.3'))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,NF90_GLOBAL,'field_date',    rtd%cyy//'-'//rtd%cmm//'-'//rtd%cdm))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,NF90_GLOBAL,'bulletin_date', rtb%cyy//'-'//rtb%cmm//'-'//rtb%cdm))
   end subroutine defNCHeader


   subroutine defNCProjectionVar(ncstate, proj_conv_fac)
   implicit none
      type(file_state) :: ncstate
   real            , intent(in) :: proj_conv_fac
   integer :: var_id

   if (trim(ncstate%cprojection)=='polar_stereographic') then
      if (abs(proj_conv_fac-1.0) < 1e-6) then
         ncstate%char_scale_fac = ''
      else
         write(ncstate%char_scale_fac,'(i4)') nint(1.0/proj_conv_fac)
      end if
   end if

   if (trim(ncstate%cprojection)=='polar_stereographic') then
      call handle_err(NF90_DEF_VAR(ncstate%ncid,'polar_stereographic',NF90_INT,var_id))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'grid_mapping_name', 'polar_stereographic'))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'latitude_of_projection_origin', 90.))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'longitude_of_projection_origin', 45.))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'scale_factor_at_projection_origin', 1.))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'false_easting', 0.))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'false_northing', 0.))
      ncstate%grid_mapping='polar_stereographic'
      ncstate%coords_mapping='longitude latitude'
   end if
   end subroutine DefNCProjectionVar
      


   subroutine NCTimeVar(ncstate,hfile)
   use mod_year_info
   use mod_hycomfile_io
   implicit none
   type  (file_state) :: ncstate
   type  (hycomfile ) :: hfile
   type  (year_info)  :: rtb, rtd
   integer :: hourref, hourval, hour_offset, var_id

   call forecastDate(hfile,rtd)
   call startDate(hfile,rtb)

   ! Calculate julian day since 1-1-1950
   hourref=datetojulian(rtb%iyy,rtb%imm,rtb%idm,1950,1,1)*24
   hourval=datetojulian(rtd%iyy,rtd%imm,rtd%idm,1950,1,1)*24
   hour_offset=hourval-hourref

   call handle_err(NF90_DEF_VAR(ncstate%ncid,'time',NF90_DOUBLE,ncstate%t_dimension_id,var_id))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'units', 'hour since 1950-1-1T00:00:00Z'))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'long_name', 'forecast time'))
   call handle_err(nf90_enddef(ncstate%ncid))
   call handle_err(NF90_PUT_VAR(ncstate%ncid,var_id,hourval))
   call handle_err(nf90_redef(ncstate%ncid))
   end subroutine NCTimeVar



   ! Define and put coordinate variables in a standard gridded netcdf file
   subroutine NCSpaceVar(ncstate,lons,lats,nlons,nlats,firstlon,firstlat,dlon,dlat,deep,ndeep) 
   !use mod_netcdf_tools
   implicit none
   type(file_state) :: ncstate
   integer :: nlons, nlats, ndeep
   real, dimension(nlons,nlats) :: lons,lats
   real, intent(in) :: firstlon, firstlat, dlon, dlat
   real, dimension(ndeep) :: deep

   character(len=20) :: vartitle
   integer :: var_id, ilon, jlat
   real :: tmpx(nlons), tmpy(nlats)
   integer :: i


   ! Put first Coordinate variable
   if (trim(ncstate%cprojection)=='regular') then
      vartitle='longitude'
      call handle_err(NF90_DEF_VAR(ncstate%ncid,'longitude',NF90_FLOAT,ncstate%lon_dimension_id,var_id))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'axis','X'))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'standard_name','longitude'))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'units','degrees_east'))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'valid_range',(/-180.0,180.0/)))
      call handle_err(nf90_enddef(ncstate%ncid))
      call handle_err(NF90_PUT_VAR(ncstate%ncid,var_id,lons(:,1)))
   else if (trim(ncstate%cprojection)=='ps_amsr' .or. trim(ncstate%cprojection)=='polar_stereographic' .or. &
            trim(ncstate%cprojection)=='native') then
      call handle_err(NF90_DEF_VAR(ncstate%ncid,'x',NF90_FLOAT,ncstate%lon_dimension_id,var_id))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'axis','X'))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'standard_name','projection_x_coordinate'))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'units',adjustl(trim(ncstate%char_scale_fac))//' km'))
      do ilon=1,nlons
         tmpx(ilon)=firstlon+(ilon-1)*dlon
      end do
      call handle_err(nf90_enddef(ncstate%ncid))
      call handle_err(NF90_PUT_VAR(ncstate%ncid,var_id,tmpx))

      call handle_err(nf90_redef(ncstate%ncid))
      call handle_err(NF90_DEF_VAR(ncstate%ncid,'longitude',NF90_FLOAT,ncstate%vars_2d,var_id))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'standard_name','longitude'))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'units','degrees_east'))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'valid_range',(/-180.0,180.0/)))
      call handle_err(nf90_enddef(ncstate%ncid))
      call handle_err(NF90_PUT_VAR(ncstate%ncid,var_id,lons))
   else 
      print *,'undef proj'
      stop '(NCSpaceVar)'
   end if

   call handle_err(nf90_redef(ncstate%ncid))
   if (trim(ncstate%cprojection)=='regular') then
      vartitle= 'latitude'
      call handle_err(NF90_DEF_VAR(ncstate%ncid,trim(vartitle),NF90_FLOAT,ncstate%lat_dimension_id,var_id))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'axis','Y'))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'units','degrees_north'))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'valid_range',(/-90.0, 90.0/)))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'standard_name','latitude'))
      call handle_err(nf90_enddef(ncstate%ncid))
      call handle_err(NF90_PUT_VAR(ncstate%ncid,var_id,real(lats(1,:),4)))
   else if (trim(ncstate%cprojection)=='ps_amsr' .or. trim(ncstate%cprojection)=='polar_stereographic' .or. &
            trim(ncstate%cprojection)=='native') then
      call handle_err(NF90_DEF_VAR(ncstate%ncid,'y',NF90_FLOAT,ncstate%lat_dimension_id,var_id))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'standard_name','projection_y_coordinate'))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'axis','Y'))
      do jlat=1,nlats
         tmpy(jlat)=firstlat+(jlat-1)*dlat
      end do
      call handle_err(nf90_enddef(ncstate%ncid))
      call handle_err(NF90_PUT_VAR(ncstate%ncid,var_id,tmpy))

      call handle_err(nf90_redef(ncstate%ncid))
      call handle_err(NF90_DEF_VAR(ncstate%ncid,'latitude',NF90_FLOAT,ncstate%vars_2d,var_id))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'standard_name','latitude'))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'units','degrees_north'))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'valid_range',(/-180.0,180.0/)))
      call handle_err(nf90_enddef(ncstate%ncid))
      call handle_err(NF90_PUT_VAR(ncstate%ncid,var_id,real(lats,kind=4)))
   else
      print *,'undef proj'
      stop '(NCSpaceVar)'
   end if

   call handle_err(nf90_redef(ncstate%ncid))
   vartitle= 'depth'
   call handle_err(NF90_DEF_VAR(ncstate%ncid,'depth',NF90_FLOAT,ncstate%k_dimension_id,var_id))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'long_name',trim(vartitle)))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'units','m'))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'valid_range',(/0., maxval(deep)+1/)))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'standard_name','depth'))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'positive','down'))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'axis','Z'))
   call handle_err(nf90_enddef(ncstate%ncid))
   call handle_err(NF90_PUT_VAR(ncstate%ncid,var_id ,deep(1:ndeep)))
   end subroutine NCSpaceVar


   subroutine NCSpaceVarStat(ncstate,nstat,deeps,ndeep)
   use mod_station
   implicit none
   type (file_state) :: ncstate
   integer, intent(in) :: nstat, ndeep
   real , intent(in) :: deeps(ndeep)

   integer :: var_id
   character(len=20) :: vartitle
   integer :: i

   ! Stations
   !call handle_err(nf90_redef(ncstate%ncid),critical=.false.)
   vartitle= 'station'
   call handle_err(NF90_DEF_VAR(ncstate%ncid,'station',NF90_INT,ncstate%lon_dimension_id,var_id))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'long_name',trim(vartitle)))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'units','[]'))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'axis','X'))
   call handle_err(nf90_enddef(ncstate%ncid))
   do i=1,nstat
      call handle_err(NF90_PUT_VAR(ncstate%ncid,var_id,i,start=(/i/)))
   end do
   !print *,'test1'
   !print *,'ok3'
   !print *,'ok2.4'

   ! Depth levels
   call handle_err(nf90_redef(ncstate%ncid))
   vartitle= 'depth'
   call handle_err(NF90_DEF_VAR(ncstate%ncid,'depth',NF90_Float,ncstate%k_dimension_id,var_id))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'long_name',trim(vartitle)))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'units','m'))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'valid_range',(/0., maxval(deeps)+1/)))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'standard_name','depth'))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'positive','down'))
   call handle_err(nf90_enddef(ncstate%ncid))
   call handle_err(NF90_PUT_VAR(ncstate%ncid,var_id  ,deeps))
   end subroutine NCSpaceVarStat



   subroutine putNCVar(ncstate,field,nx,ny,nz,vname,ndims,gridrotate)
   use mod_parameters
   use mod_hycomfile_io
   implicit none
   type (file_state) :: ncstate
   integer, intent(in) :: nx,ny,nz
   integer, intent(in) :: ndims     ! Dimension of this variable (space-time)
   real :: field(nx,ny,nz)
   character(len=*), intent(in) :: vname
   logical, intent(in) :: gridrotate

   integer :: var_dims(ndims)
   integer :: fieldint(nx,ny,nz)
   integer :: varid, countlo, counthi 
   character(len=80) :: stdname, units, vnamenc, cellmethod
   integer :: fillval
   real    :: scale_factor, add_offset, limits(2)


   ! NB - ndims are space-time dimensions
   if (ncstate%stations) then

      if (ndims==1) then
         var_dims=ncstate%vars_2d(1:1)
      else if (ndims==2) then
         var_dims=ncstate%vars_3d(1:2)
      else if (ndims==3) then
         var_dims=ncstate%vars_4d(1:3)
      else
         print *,'putNCVar: Wrong number of dims ',ndims
         stop
      end if

   else

      if (ndims==2) then
         var_dims=ncstate%vars_2d
      else if (ndims==3) then
         var_dims=ncstate%vars_3d
      else if (ndims==4) then
         var_dims=ncstate%vars_4d
      else
         print *,'putNCVar: Wrong number of dims ',ndims
         stop
      end if
   end if

   call handle_err(nf90_redef(ncstate%ncid))

   ! Retrieve variable info in CF formulation
   call  netcdfInfo(vname,gridrotate,stdname,units,vnamenc,cellmethod,limits)

   ! No variable packing in this case
   if (limits(2)-limits(1) ==0 ) then
   !if (.true.) then
      call handle_err(NF90_DEF_VAR(ncstate%ncid,trim(vnamenc),NF90_FLOAT,var_dims,varid))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,varid,'_FillValue'   ,real(undef,kind=4)))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,varid,'missing_value'   ,real(undef,kind=4)))

   ! Variable packing
   else
      add_offset  =(limits(2)+limits(1))/2.
      scale_factor=(limits(2)-limits(1))/2**16
      fillval     =-(2**16-2)/2
      call handle_err(NF90_DEF_VAR(ncstate%ncid,trim(vnamenc),NF90_SHORT,var_dims,varid))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,varid,'_FillValue'   ,int(fillval,kind=2)))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,varid,'missing_value',fillval))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,varid,'add_offset'   ,add_offset))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,varid,'scale_factor' ,scale_factor))
   end if
   if (trim(units)/='')  &
      call handle_err(NF90_PUT_ATT(ncstate%ncid,varid,'units' ,     trim(Units)))
   if (trim(StdName)/='')  &
      call handle_err(NF90_PUT_ATT(ncstate%ncid,varid,'standard_name', trim(StdName)))
   if (trim(ncstate%grid_mapping)/='') &
      call handle_err(NF90_PUT_ATT(ncstate%ncid,varid,'grid_mapping',trim(ncstate%grid_mapping)))
   if (trim(ncstate%coords_mapping)/='')  &
      call handle_err(NF90_PUT_ATT(ncstate%ncid,varid,'coordinates',trim(ncstate%coords_mapping)))
   if (trim(cellmethod)/='')  call handle_err(NF90_PUT_ATT(ncstate%ncid,varid,'cell_methods',trim(cellmethod)))
   call handle_err(nf90_enddef(ncstate%ncid))


   if (limits(2)-limits(1) ==0 ) then
      call handle_err(NF90_PUT_VAR(ncstate%ncid,varid,real(field,kind=4)))
   else
      ! Check for over/undershoot
      counthi=count(limits(2)<field.and.abs(field-undef)>1e-4)
      countlo=count(limits(1)>field.and.abs(field-undef)>1e-4)
      if (counthi>0 .and. countlo>0. ) then
         print *,'Over/Undershoot for variable '//vname
         print '(a,2e14.6,2i10)','min, max, countlo, counthi:', &
            minval(field,mask=abs(undef-field)>1e-4),maxval(field,mask=abs(undef-field)>1e-4), &
            countlo,counthi
      end if
      where (field/=undef)
         !fieldint = floor((max(min(field,limits(2)),limits(1)) - add_offset) / scale_factor)
         fieldint = min(max(fillval,nint((field - add_offset)/scale_factor)),-fillval)
      elsewhere
         fieldint=fillval
      endwhere
      call handle_err(NF90_PUT_VAR(ncstate%ncid,varid,fieldint))
   end if
   end subroutine putNCVar




   subroutine handle_err(errcode,critical)
      use netcdf
      implicit none
      integer, intent(in) :: errcode
      logical, optional, intent(in) :: critical
      logical :: crit

      crit=.true.
      if (present(critical)) crit=critical

      if (errcode/=NF90_NOERR) then
         write(6,'(a)') NF90_STRERROR(errcode)
         if (crit) stop '(handle_err)'
      end if
   end subroutine
end module

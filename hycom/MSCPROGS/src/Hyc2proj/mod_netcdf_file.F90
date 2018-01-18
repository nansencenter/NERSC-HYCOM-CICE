module mod_netcdf_file
use netcdf

type file_state
   integer :: ncid=-1
   integer :: lon_dimension_id, lat_dimension_id, k_dimension_id, t_dimension_id
   integer :: char_dimension_id
   integer :: vars_1d(1),vars_2d(2), vars_3d(3), vars_4d(4)
   character(len=20) ::  grid_mapping='', coords_mapping='', char_scale_fac=''
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

   subroutine DefNCDim(ncstate,ndeep)
   use mod_toproj, only : nxp, nyp, cprojection
   implicit none
   type(file_state), intent(inout) :: ncstate
   integer, intent(in) ::  ndeep
   if (trim(cprojection)=='regular' .or. trim(cprojection)=='mercator' ) then
      call handle_err(NF90_DEF_DIM(ncstate%ncid,'longitude',nxp,ncstate%lon_dimension_id))
      call handle_err(NF90_DEF_DIM(ncstate%ncid,'latitude' ,nyp,ncstate%lat_dimension_id))
   else if (trim(cprojection)=='ps_amsr' .or. trim(cprojection)=='polar_stereographic' .or. &
            trim(cprojection)=='native') then
      call handle_err(NF90_DEF_DIM(ncstate%ncid,'x',nxp,ncstate%lon_dimension_id))
      call handle_err(NF90_DEF_DIM(ncstate%ncid,'y',nyp,ncstate%lat_dimension_id))
   else
      print *,'Unknown projection '//trim(cprojection)
   end if
   call handle_err(NF90_DEF_DIM(ncstate%ncid,'depth'  ,ndeep,ncstate%k_dimension_id))
   call handle_err(NF90_DEF_DIM(ncstate%ncid,'time',NF90_Unlimited,ncstate%t_dimension_id))
   ! Define dimension "holders"
   ncstate%vars_4d=(/ncstate%lon_dimension_id,ncstate%lat_dimension_id,ncstate%k_dimension_id,ncstate%t_dimension_id/)
   ncstate%vars_3d=(/ncstate%lon_dimension_id,ncstate%lat_dimension_id,ncstate%t_dimension_id/)
   ncstate%vars_2d=(/ncstate%lon_dimension_id,ncstate%lat_dimension_id/)
   end subroutine




#if ! defined (MULT_GRP_PER_FILE)
   subroutine DefNCDimStat(ncstate,nstat,ndeep)
   implicit none
   type(file_state), intent(inout) :: ncstate
   integer, intent(in) :: nstat, ndeep
#else
   subroutine DefNCDimStat(ncstate,ngroup,nstat,ndeep)
   implicit none
   type(file_state), intent(inout) :: ncstate
   integer, intent(in) :: ngroup,nstat, ndeep
#endif
#if defined (MULT_GRP_PER_FILE)
   call handle_err(NF90_DEF_DIM(ncstate%ncid,'group',ngroup,ncstate%lon_dimension_id))
   call handle_err(NF90_DEF_DIM(ncstate%ncid,'station',nstat,ncstate%lat_dimension_id))
   call handle_err(NF90_DEF_DIM(ncstate%ncid,'depth'  ,ndeep,ncstate%k_dimension_id))
   call handle_err(NF90_DEF_DIM(ncstate%ncid,'time',NF90_Unlimited,ncstate%t_dimension_id))
   call handle_err(NF90_DEF_DIM(ncstate%ncid,'N_40',40,ncstate%char_dimension_id))
   ncstate%vars_4d=(/ncstate%lon_dimension_id,ncstate%lat_dimension_id,ncstate%k_dimension_id,ncstate%t_dimension_id/)
   ncstate%vars_3d=(/ncstate%lon_dimension_id,ncstate%lat_dimension_id,ncstate%t_dimension_id/)
   ncstate%vars_2d=(/ncstate%lon_dimension_id,ncstate%lat_dimension_id/)
#else
   call handle_err(NF90_DEF_DIM(ncstate%ncid,'station',nstat,ncstate%lat_dimension_id))  
   call handle_err(NF90_DEF_DIM(ncstate%ncid,'depth'  ,ndeep,ncstate%k_dimension_id))
   call handle_err(NF90_DEF_DIM(ncstate%ncid,'time',NF90_Unlimited,ncstate%t_dimension_id))
   ncstate%vars_4d=(/ncstate%lat_dimension_id,ncstate%k_dimension_id,ncstate%t_dimension_id,-1/)
   ncstate%vars_3d=(/ncstate%lat_dimension_id,ncstate%t_dimension_id,-1/)
   ncstate%vars_2d=(/ncstate%lat_dimension_id,-1/)
#endif
   ncstate%stations=.true.
   end subroutine


   subroutine DefNCHeader(ncstate,hfile,cver,igrp,grpname)
   use mod_year_info
   use mod_hycomfile_io
   implicit none
   type(file_state) :: ncstate
   type(hycomfile)  :: hfile
   character(len=*), intent(in) :: cver
   character(len=*), intent(in), optional :: grpname
   integer, intent(in), optional :: igrp 

   type (year_info) :: rtd, rtb 
   character(len=8)   :: dateinfo
   character(len=80)  :: c80

   call forecastDate(hfile,rtd)
   call startDate(hfile,rtb)

   call handle_err(NF90_PUT_ATT(ncstate%ncid,NF90_GLOBAL,'title', 'Pilot MyOcean reanalysis by TOPAZ4 (2003-2008)'))
   !call handle_err(NF90_PUT_ATT(ncstate%ncid,NF90_GLOBAL,'comment', trim(comment)))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,NF90_GLOBAL,'institution',  'NERSC, Thormoehlens gate 47, N-5006 Bergen, Norway'))
   call date_and_time(date=dateinfo)
   call handle_err(NF90_PUT_ATT(ncstate%ncid,NF90_GLOBAL,'history', dateinfo(1:8)//':Created by program hyc2proj, version '//cver))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,NF90_GLOBAL,'source', 'NERSC-HYCOM model fields'))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,NF90_GLOBAL,'references', 'http://topaz.nersc.no'))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,NF90_GLOBAL,'field_type',  'Files based on file type '//trim(hfile%ftype)))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,NF90_GLOBAL,'Conventions',  'CF-1.4'))

   if (hfile%ftype=='archv'.or.hfile%ftype=='archv_wav') then
      !put time in field date
      write(c80,'(i2.2,a,i2.2,a,i2.2,a)') hfile%ihour,':',hfile%imin,':',hfile%isec,'Z'
      c80   = rtd%cyy//'-'//rtd%cmm//'-'//rtd%cdm//'T'//trim(c80)
      call handle_err(NF90_PUT_ATT(ncstate%ncid,NF90_GLOBAL,'field_date',trim(c80)))
   else
      !no time in field date
      call handle_err(NF90_PUT_ATT(ncstate%ncid,NF90_GLOBAL,'field_date',    rtd%cyy//'-'//rtd%cmm//'-'//rtd%cdm))
   end if

   if (hfile%ftype=='restart') then
      !add field time to restarts
      write(c80,'(i2.2,a,i2.2,a,i2.2,a)') hfile%ihour,':',hfile%imin,':',hfile%isec,'Z'
      call handle_err(NF90_PUT_ATT(ncstate%ncid,NF90_GLOBAL,'field_time',trim(c80)))
   end if

   if (hfile%ftype=='archv'.or.hfile%ftype=='archv_wav') then
      !put time in bulletin date
      write(c80,'(i2.2,a,i2.2,a,i2.2,a)') hfile%start_ihour,':',hfile%start_imin,':',hfile%start_isec,'Z'
      c80   = rtb%cyy//'-'//rtb%cmm//'-'//rtb%cdm//'T'//trim(c80)
      call handle_err(NF90_PUT_ATT(ncstate%ncid,NF90_GLOBAL,'bulletin_date',trim(c80)))
   else
      !no time in bulletin date
      call handle_err(NF90_PUT_ATT(ncstate%ncid,NF90_GLOBAL,'bulletin_date', rtb%cyy//'-'//rtb%cmm//'-'//rtb%cdm))
   endif

   if (present(grpname) .and. present(igrp) )then
      call handle_err(NF90_PUT_ATT(ncstate%ncid,NF90_GLOBAL,'station_group_number',igrp))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,NF90_GLOBAL,'station_group_name', trim(grpname)))
   endif
   end subroutine defNCHeader


   subroutine defNCProjectionVar(ncstate)
   use mod_toproj
   implicit none
   type(file_state) :: ncstate
   integer :: var_id

   if (trim(cprojection)=='polar_stereographic') then
      if (abs(proj_conv_fac-1.0) < 1e-6) then
         ncstate%char_scale_fac = ''
      else
         write(ncstate%char_scale_fac,'(i4)') nint(1.0/proj_conv_fac)
      end if
   end if

   if (trim(cprojection)=='polar_stereographic') then
      call handle_err(NF90_DEF_VAR(ncstate%ncid,'stereographic',NF90_INT,var_id))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'grid_mapping_name',  &
           'polar_stereographic'))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'latitude_of_projection_origin',&
            polar_stereographic_central_lat))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'longitude_of_projection_origin', &
            polar_stereographic_central_lon))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'scale_factor_at_projection_origin', 1.))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'straight_vertical_longitude_from_pole',-45.))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'false_easting', 0.))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'false_northing', 0.))
      ncstate%grid_mapping='stereographic'
      ncstate%coords_mapping='longitude latitude'
   else if (trim(cprojection)=='mercator') then
      call handle_err(NF90_DEF_VAR(ncstate%ncid,'mercator',NF90_INT,var_id))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'grid_mapping_name',  &
           'mercator'))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'longitude_of_projection_origin', &
            mercator_central_lon))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'scale_factor_at_projection_origin', 1.))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'false_easting', 0.))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'false_northing', 0.))
      ncstate%grid_mapping='mercator'
   end if
   end subroutine DefNCProjectionVar
      


   subroutine NCTimeVar(ncstate,hfile)
   use mod_year_info
   use mod_hycomfile_io
   implicit none
   type  (file_state) :: ncstate
   type  (hycomfile ) :: hfile
   type  (year_info)  :: rtb, rtd
   integer :: hourref, hourval, hour_offset, var_id,timeval
   !!
   !character(len=*),parameter :: tunits = 'hours since 1950-1-1T00:00:00Z'
   !integer,parameter  :: yr0        = 1950   !reference year
   !integer,parameter  :: unit_fac   = 1 !hours to hours
   !character(len=*),parameter :: tunits = 'seconds since 1970-1-1T00:00:00Z'
   !integer,parameter  :: yr0        = 1970   !reference year
   !integer,parameter  :: unit_fac   = 60*60  !hours to seconds
   character(len=80) :: tunits   !string for time units
   integer           :: year0    !reference year for time value


   call forecastDate(hfile,rtd)
   call startDate(hfile,rtb)

   if (hfile%ftype=='archv') then
      tunits   = 'seconds since 1970-1-1T00:00:00Z'
      year0    = 1970   !reference year
   elseif (hfile%ftype=='archv_wav') then
      tunits   = 'seconds since 1970-1-1T00:00:00Z'
      year0    = 1970   !reference year
   else
      tunits   = 'hour since 1950-1-1T00:00:00Z'
      year0    = 1950   !reference year
   endif

   ! Calculate julian day since 1-1-1950
   !LB: TODO: Change to 'seconds since 1-1-1970' for compliance to standards
   hourref=datetojulian(rtb%iyy,rtb%imm,rtb%idm+1,year0,1,1)*24
   hourval=datetojulian(rtd%iyy,rtd%imm,rtd%idm+1,year0,1,1)*24
   hour_offset=hourval-hourref

   if (hfile%ftype=='archv') then
      !!add "time of day"
      timeval  = 3600*(hourval+hfile%ihour)+60*hfile%imin+hfile%isec
   elseif (hfile%ftype=='archv_wav') then
      !!add "time of day"
      timeval  = 3600*(hourval+hfile%ihour)+60*hfile%imin+hfile%isec
   elseif (hfile%ftype=='restart') then
      !!add "time of day"
      timeval  = hourval+hfile%ihour
   else
      !!otherwise start of day is the "time"
      timeval  = hourval
   endif

   call handle_err(NF90_DEF_VAR(ncstate%ncid,'time',NF90_DOUBLE,ncstate%t_dimension_id,var_id))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'units',tunits))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'long_name', 'forecast time'))
   call handle_err(nf90_enddef(ncstate%ncid))
   call handle_err(NF90_PUT_VAR(ncstate%ncid,var_id,timeval))
   call handle_err(nf90_redef(ncstate%ncid))
   end subroutine NCTimeVar



   ! Define and put coordinate variables in a standard gridded netcdf file
   subroutine NCSpaceVar(ncstate,deep,ndeep) 
   use mod_toproj, only : firstx, firsty, nxp, nyp, dx, dy, cprojection, lons, lats, &
                          isRectilinear
   implicit none
   type(file_state) :: ncstate
   integer :: ndeep
   real, dimension(ndeep) :: deep

   character(len=20) :: vartitle
   integer :: var_id, ix, jy
   real :: tmpx(nxp), tmpy(nyp)
   integer :: i


   ! Put first Coordinate variable for rectilinear grids - isRectilinear is
   ! function in mod_toproj
   if (isRectilinear()) then
      call handle_err(NF90_DEF_VAR(ncstate%ncid,'longitude',NF90_FLOAT,ncstate%lon_dimension_id,var_id))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'axis','X'))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'standard_name','longitude'))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'units','degrees_east'))
      call handle_err(nf90_enddef(ncstate%ncid))
      call handle_err(NF90_PUT_VAR(ncstate%ncid,var_id,lons(:,1)))
      call handle_err(nf90_redef(ncstate%ncid))
   end if

   if (trim(cprojection)=='ps_amsr'.or. trim(cprojection)=='polar_stereographic' .or. &
       trim(cprojection)=='native' .or. trim(cprojection)=='mercator' ) then
      call handle_err(NF90_DEF_VAR(ncstate%ncid,'x',NF90_FLOAT,ncstate%lon_dimension_id,var_id))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'axis','X'))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'standard_name','projection_x_coordinate'))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'units',adjustl(trim(ncstate%char_scale_fac))//' km'))
      do ix=1,nxp
         tmpx(ix)=firstx+(ix-1)*dx
      end do
      call handle_err(nf90_enddef(ncstate%ncid))
      call handle_err(NF90_PUT_VAR(ncstate%ncid,var_id,tmpx))


      ! Put longitude for non- rectilinear grids
      if (.not.isRectilinear()) then
         call handle_err(nf90_redef(ncstate%ncid))
         call handle_err(NF90_DEF_VAR(ncstate%ncid,'longitude',NF90_FLOAT,ncstate%vars_2d,var_id))
         call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'standard_name','longitude'))
         call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'units','degrees_east'))
         call handle_err(nf90_enddef(ncstate%ncid))
         call handle_err(NF90_PUT_VAR(ncstate%ncid,var_id,lons))
      end if
      call handle_err(nf90_redef(ncstate%ncid))
   end if

   ! Put 2nd Coordinate variable for rectilinear grids
   if (isRectilinear()) then
      call handle_err(NF90_DEF_VAR(ncstate%ncid,'latitude',NF90_FLOAT,ncstate%lat_dimension_id,var_id))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'axis','Y'))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'units','degrees_north'))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'standard_name','latitude'))
      call handle_err(nf90_enddef(ncstate%ncid))
      call handle_err(NF90_PUT_VAR(ncstate%ncid,var_id,real(lats(1,:),4)))
      call handle_err(nf90_redef(ncstate%ncid))
   end if

   if (trim(cprojection)=='ps_amsr' .or. trim(cprojection)=='polar_stereographic' .or. &
       trim(cprojection)=='native'  .or. trim(cprojection)=='mercator') then
      call handle_err(NF90_DEF_VAR(ncstate%ncid,'y',NF90_FLOAT,ncstate%lat_dimension_id,var_id))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'standard_name','projection_y_coordinate'))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'axis','Y'))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'units',adjustl(trim(ncstate%char_scale_fac))//' km'))
      do jy=1,nyp
         tmpy(jy)=firsty+(jy-1)*dy
      end do
      call handle_err(nf90_enddef(ncstate%ncid))
      call handle_err(NF90_PUT_VAR(ncstate%ncid,var_id,tmpy))

      ! Put latitude for non-rectilinear grids
      if (.not.isRectilinear()) then
         call handle_err(nf90_redef(ncstate%ncid))
         call handle_err(NF90_DEF_VAR(ncstate%ncid,'latitude',NF90_FLOAT,ncstate%vars_2d,var_id))
         call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'standard_name','latitude'))
         call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'units','degrees_north'))
         call handle_err(nf90_enddef(ncstate%ncid))
         call handle_err(NF90_PUT_VAR(ncstate%ncid,var_id,real(lats,kind=4)))
      end if
      call handle_err(nf90_redef(ncstate%ncid))
   end if

   call handle_err(NF90_DEF_VAR(ncstate%ncid,'depth',NF90_FLOAT,ncstate%k_dimension_id,var_id))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'long_name','depth'))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'units','m'))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'standard_name','depth'))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'positive','down'))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'axis','Z'))
   call handle_err(nf90_enddef(ncstate%ncid))
   call handle_err(NF90_PUT_VAR(ncstate%ncid,var_id ,deep(1:ndeep)))
   end subroutine NCSpaceVar


#if ! defined (MULT_GRP_PER_FILE)
   subroutine NCSpaceVarStat(ncstate,nstat,deeps,ndeep)
#else
   subroutine NCSpaceVarStat(ncstate,ngroup,nstat,deeps,ndeep)
#endif
   use mod_station , only: groupname
   implicit none
   type (file_state) :: ncstate
#if ! defined (MULT_GRP_PER_FILE)
   integer, intent(in) :: nstat, ndeep
#else
   integer, intent(in) :: nstat, ndeep, ngroup
#endif
   real , intent(in) :: deeps(ndeep)

   integer :: var_id
   character(len=20) :: vartitle
   integer :: i

   ! Stations
   call handle_err(NF90_DEF_VAR(ncstate%ncid,'station',NF90_INT,ncstate%lat_dimension_id,var_id))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'long_name','Station index'))
   call handle_err(nf90_enddef(ncstate%ncid))
   do i=1,nstat
      call handle_err(NF90_PUT_VAR(ncstate%ncid,var_id,i,start=(/i/)))
   end do


#if defined (MULT_GRP_PER_FILE)
   ! groups
   call handle_err(nf90_redef(ncstate%ncid),critical=.false.)
   call handle_err(NF90_DEF_VAR(ncstate%ncid,'groups',NF90_INT,ncstate%lon_dimension_id,var_id))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'long_name','Index of station group'))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'units','[]'))
   call handle_err(nf90_enddef(ncstate%ncid))
   do i=1,ngroup
      call handle_err(NF90_PUT_VAR(ncstate%ncid,var_id,i,start=(/i/)))
   end do

   call handle_err(nf90_redef(ncstate%ncid),critical=.false.)
   call handle_err(NF90_DEF_VAR(ncstate%ncid,'group_name',NF90_CHAR, &
      (/ncstate%char_dimension_id,ncstate%lon_dimension_id/),var_id))
   call handle_err(NF90_PUT_ATT(ncstate%ncid,var_id,'long_name','Name of station group'))
   call handle_err(nf90_enddef(ncstate%ncid))
   do i=1,ngroup
      call handle_err(NF90_PUT_VAR(ncstate%ncid,var_id,groupname(i)(1:40),start=(/1,i/)))
   end do
#endif


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
   character(len=85) :: stdname, longname, units, vnamenc, cellmethod
   integer :: fillval
   real    :: scale_factor, add_offset, limits(2)


   ! NB - ndims are space-time dimensions
   if (ncstate%stations) then

#if defined (MULT_GRP_PER_FILE)
   !NB - same as in hyc2proj case
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
#else
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
#endif
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
   call  netcdfInfo(vname,gridrotate,stdname,longname,units,vnamenc,cellmethod,limits)

   ! No variable packing in this case
   if (limits(2)-limits(1) ==0 ) then
      call handle_err(NF90_DEF_VAR(ncstate%ncid,trim(vnamenc),NF90_FLOAT,var_dims,varid))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,varid,'_FillValue'   ,real(undef,kind=4)))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,varid,'missing_value'   ,real(undef,kind=4)))

   ! Variable packing
   else
      add_offset  =(limits(2)+limits(1))/2.
      scale_factor=(limits(2)-limits(1))/(2**16-10)
      fillval     =-(2**16-2)/2
      call handle_err(NF90_DEF_VAR(ncstate%ncid,trim(vnamenc),NF90_SHORT,var_dims,varid))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,varid,'_FillValue'   ,int(fillval,kind=2)))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,varid,'missing_value',int(fillval,kind=2)))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,varid,'add_offset'   ,add_offset))
      call handle_err(NF90_PUT_ATT(ncstate%ncid,varid,'scale_factor' ,scale_factor))
   end if
   if (trim(units)/='')  &
      call handle_err(NF90_PUT_ATT(ncstate%ncid,varid,'units' ,     trim(Units)))
   if (trim(StdName)/='')  &
      call handle_err(NF90_PUT_ATT(ncstate%ncid,varid,'standard_name', trim(StdName)))
   if (trim(longName)/='')  &
      call handle_err(NF90_PUT_ATT(ncstate%ncid,varid,'long_name', trim(longName)))
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

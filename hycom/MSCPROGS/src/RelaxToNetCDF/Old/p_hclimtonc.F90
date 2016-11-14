program p_hclimtonc
use mod_xc
use mod_za
use mod_grid
use mod_toproj
use mod_confmap
use netcdf
use m_handle_err
implicit none

character(len=*), parameter :: cver='0.1'
real, parameter :: undef=-1e14

integer :: ios
integer :: k,k2
integer :: ilon,jlat
integer,allocatable :: kindx(:)
character(len=80) :: tmparg,infile,ncfile
character(len= 2) :: cmonth
character(len= 8) :: dateinfo
character(len=  4) :: char_fac
character(len=40)  :: grid_mapping,coords_mapping, vartitle

integer :: ncid
integer :: lon_dimension_id, lat_dimension_id, k_dimension_id
integer :: vars_2d(2),vars_3d(3), var_id

real, allocatable, dimension(:) :: tmpx, tmpy

real :: hmina,hmaxa

integer :: numdepths
integer :: month
real   :: hcdepth
real,allocatable   :: indepth(:), regu2d(:,:), hy2d(:,:), regubathy(:,:)
#if defined (IARGC)
integer*4, external :: iargc
#endif


call xcspmd()
call zaiost()
call initconfmap(idm,jdm)    ! Init conf map
call get_grid()       ! Get grid 
call proj()       ! Init projection
call proj_ini(undef)  ! Init projection


allocate(tmpx(nlons))
allocate(tmpy(nlats))
allocate(regu2d   (nlons,nlats))        ! 2D vars on regular gri
allocate(regubathy(nlons,nlats))        ! 2D vars on regular gri
allocate(hy2d     (idm,jdm))        ! Holds 3D vars on original grid




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! -- Get arguments
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (iargc()>=2) then
   call getarg(1,tmparg) ; read(tmparg,*)  month

   numdepths=iargc()-1
   allocate(kindx  (numdepths))
   allocate(indepth(numdepths))
   do k=1,numdepths
      call getarg(k+1,tmparg) ; read(tmparg,*)  indepth(k)
   end do

else
   print *,'Routine reads horizontally interpolated z-level climatology'
   print *,'produced by hycom "ALL" routines. The climatology is interpolated'
   print *,'to grid specified by proj.in - similar to "hyc2proj" routines'
   print *,'input is month and depthlevel(s)'
   print *
   print *,'Usage : hclimtonc month depthlevel(s)'
   stop 
end if


! open hclim_depthlevels to see if we have a match agains depth
open(10,file='hclim_depthlevels',status='old')

do k2=1,numdepths
   ios=0
   k=1
   kindx(k2)=0
   do while (ios==0 .and. kindx(k2)==0)
      read(10,*,iostat=ios) hcdepth
      if (abs(hcdepth-indepth(k2))<.01) then
         kindx(k2)=k
      end if
      k=k+1
   end do

   if (kindx(k2)==0) then
      print *,'Could not find depth index ',indepth(k2)
      print *,'Note that depth levels mus match exactly those of the climatology'
      print *,'See hclim_depthlevels for available depth levels'
      close(10)
      stop
   end if
   print *,'climatology k index is :',kindx(k2),' for depth ',indepth(k2)
   rewind(10)
end do
close(10)









! File names, input and netcdf file
write(cmonth,'(i2.2)') month
ncfile='clim_sig0_m'//cmonth//'.nc'

! Open netcdf file
if (NF90_CREATE(ncfile,NF90_CLOBBER,ncid) /= NF90_NOERR) then
   print *,'An error occured when opening the netcdf file'
   stop '(hclimtonc)'
end if

! Define dimensions -- idm,jdm,kdm
! Following MERCATOR name definitions for netcdf var name
if (trim(cprojection)=='regular') then
   call handle_err(NF90_DEF_DIM(ncid,'longitude',nlons,lon_dimension_id))
   call handle_err(NF90_DEF_DIM(ncid,'latitude',nlats,lat_dimension_id))
else if (trim(cprojection)=='ps_amsr' .or. trim(cprojection)=='polar_stereographic' .or. &
         trim(cprojection)=='native') then
   call handle_err(NF90_DEF_DIM(ncid,'x',nlons,lon_dimension_id))
   call handle_err(NF90_DEF_DIM(ncid,'y',nlats,lat_dimension_id))
else
   print *,'Unknown projection '//trim(cprojection)
   stop
end if
call handle_err(NF90_DEF_DIM(ncid,'depth'  ,numdepths,k_dimension_id))

! Define dimension "holders"
vars_3d=(/lon_dimension_id,lat_dimension_id,k_dimension_id/)
vars_2d=(/lon_dimension_id,lat_dimension_id/)

! Define variables -- put attributes
!Global attributes
grid_mapping=''
coords_mapping=''
call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'title',        &
      'GDEM climatology - on TOPAZ2 grid'))
call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'institution',        &
   'NERSC, Thormoehlens gate 47, N-5006 Bergen, Norway'))
call date_and_time(date=dateinfo)
call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'history', &
   dateinfo(1:8)//':Created by program hclimtonc, version '//cver))
call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'Conventions', &
                'CF-1.0'))

if (trim(cprojection)=='polar_stereographic') then
   if (abs(proj_conv_fac-1.0) < 1e-6) then
      char_fac = ''
   else
      write(char_fac,'(i4)') nint(1.0/proj_conv_fac)
   end if
end if

if (trim(cprojection)=='polar_stereographic') then
   call handle_err(NF90_DEF_VAR(ncid,'polar_stereographic',NF90_INT,var_id))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'grid_mapping_name', &
      'polar_stereographic'))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'latitude_of_projection_origin', &
      90.))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'longitude_of_projection_origin', &
      45.))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'scale_factor_at_projection_origin', &
      1.))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'false_easting', &
      0.))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'false_northing', &
      0.))
      grid_mapping='polar_stereographic'
      coords_mapping='longitude latitude'
end if

! --- Start defining vars - put them as we move along
! Longitude
where (lons<-180.) lons=lons+360.
where (lons> 180.) lons=lons-360.
if (trim(cprojection)=='regular') then
   vartitle='longitude'
   call handle_err(NF90_DEF_VAR(ncid,'longitude',NF90_Float,lon_dimension_id,var_id))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'axis','X'))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'standard_name','longitude'))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'units','degrees_east'))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'valid_range',(/-180.0,180.0/)))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'long_name',trim(vartitle)))
   call handle_err(nf90_enddef(ncid))
   call handle_err(NF90_PUT_VAR(ncid,var_id,lons(:,1)))
else if (trim(cprojection)=='ps_amsr' .or. trim(cprojection)=='polar_stereographic') then
   !print *,'ps_amsr x'
   vartitle='x'
   call handle_err(NF90_DEF_VAR(ncid,trim(vartitle),NF90_Float,lon_dimension_id,var_id))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'axis','X'))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'standard_name','projection_x_coordinate'))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'units',adjustl(trim(char_fac))//' km'))
   do ilon=1,nlons
      tmpx(ilon)=firstlon+(ilon-1)*dlon
   end do
   call handle_err(nf90_enddef(ncid))
   call handle_err(NF90_PUT_VAR(ncid,var_id,tmpx))

   call handle_err(nf90_redef(ncid))
   vartitle='longitude'
   call handle_err(NF90_DEF_VAR(ncid,trim(vartitle),NF90_Float,vars_2d,var_id))
   !call handle_err(NF90_PUT_ATT(ncid,var_id,'standard_name','grid_longitude'))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'units','degrees'))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'valid_range',(/-180.0,180.0/)))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'long_name',trim(vartitle)))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'missing_value',undef))
   call handle_err(nf90_enddef(ncid))
   call handle_err(NF90_PUT_VAR(ncid,var_id,lons(:,:)))
else if (trim(cprojection)=='native') then
   vartitle='x'
   call handle_err(NF90_DEF_VAR(ncid,trim(vartitle),NF90_Float,lon_dimension_id,var_id))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'axis','X'))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'standard_name','projection_x_coordinate'))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'units',adjustl(trim(char_fac))//' km'))
   do ilon=1,nlons
      tmpx(ilon)=firstlon+(ilon-1)*dlon
   end do
   call handle_err(nf90_enddef(ncid))
   call handle_err(NF90_PUT_VAR(ncid,var_id,tmpx))
else 
   print *,'Unknown projection '//trim(cprojection)
   stop
end if


! Latitude
call handle_err(nf90_redef(ncid))
if (trim(cprojection)=='regular') then
   vartitle= 'latitude'
   call handle_err(NF90_DEF_VAR(ncid,trim(vartitle),NF90_Float,lat_dimension_id,var_id))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'axis','Y'))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'units','degrees_north'))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'valid_range',(/-90.0, 90.0/)))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'standard_name','latitude'))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'long_name',trim(vartitle)))
   call handle_err(nf90_enddef(ncid))
   call handle_err(NF90_PUT_VAR(ncid,var_id,real(lats(1,:),4)))
else if (trim(cprojection)=='ps_amsr' .or. trim(cprojection)=='polar_stereographic') then
   !print *,'ps_amsr y'
   vartitle= 'y'
   call handle_err(NF90_DEF_VAR(ncid,trim(vartitle),NF90_Float,lat_dimension_id,var_id))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'standard_name','projection_y_coordinate'))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'axis','Y'))
   do jlat=1,nlats
      tmpy(jlat)=firstlat+(jlat-1)*dlat
   end do
   call handle_err(nf90_enddef(ncid))
   call handle_err(NF90_PUT_VAR(ncid,var_id,real(tmpy,4)))
   call handle_err(nf90_redef(ncid))
   vartitle= 'latitude'
   call handle_err(NF90_DEF_VAR(ncid,trim(vartitle),NF90_Float,vars_2d,var_id))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'units','degrees'))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'valid_range',(/-90.0, 90.0/)))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'standard_name',trim(vartitle)))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'long_name',trim(vartitle)))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'missing_value',undef))
   call handle_err(nf90_enddef(ncid))
   call handle_err(NF90_PUT_VAR(ncid,var_id,real(lats(:,:),4)))
else if (trim(cprojection)=='native') then
   !print *,'ps_amsr y'
   vartitle= 'y'
   call handle_err(NF90_DEF_VAR(ncid,trim(vartitle),NF90_Float,lat_dimension_id,var_id))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'standard_name','projection_y_coordinate'))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'axis','Y'))
   do jlat=1,nlats
      tmpy(jlat)=firstlat+(jlat-1)*dlat
   end do
   call handle_err(nf90_enddef(ncid))
   call handle_err(NF90_PUT_VAR(ncid,var_id,real(tmpy,4)))
else
   print *,'Unknown projection '//trim(cprojection)
   stop
end if

call handle_err(nf90_redef(ncid))
vartitle= 'model_depth'
call handle_err(NF90_DEF_VAR(ncid,trim(vartitle),NF90_Float,vars_2d,var_id))
call handle_err(NF90_PUT_ATT(ncid,var_id,'units','meters'))
call handle_err(NF90_PUT_ATT(ncid,var_id,'valid_range',(/0.0, 9000.0/)))
call handle_err(NF90_PUT_ATT(ncid,var_id,'long_name',trim(vartitle)))
call handle_err(NF90_PUT_ATT(ncid,var_id,'missing_value',undef))
call handle_err(NF90_PUT_ATT(ncid,var_id,'standard_name','sea_floor_depth_below_sea_level'))
call handle_err(NF90_PUT_ATT(ncid,var_id,'coordinates',trim(coords_mapping)))
call handle_err(nf90_enddef(ncid))
call to_proj(depths,regubathy)
call handle_err(NF90_PUT_VAR(ncid,var_id,regubathy))


! Define temperature variable
call handle_err(nf90_redef(ncid))
call handle_err(NF90_DEF_VAR(ncid,'temperature',NF90_Float,vars_3d,var_id))
call handle_err(NF90_PUT_ATT(ncid,var_id,'long_name','temperature'))
call handle_err(NF90_PUT_ATT(ncid,var_id,'units','degrees_celsius'))
call handle_err(NF90_PUT_ATT(ncid,var_id,'standard_name','sea_water_potential_temperature'))
call handle_err(NF90_PUT_ATT(ncid,var_id,'missing_value',undef))
call handle_err(nf90_enddef(ncid))

! Open Temperature infile
print *,'Processing temperature'
infile='temp_sig0_m'//cmonth//'.a'
do k=1,numdepths

   print *,'Processing depth level ',indepth(k)
   

   call zaiopf(trim(infile),'old',911)

   if (kindx(k)>1) then
      do k2=1,kindx(k)-1
         call zaiosk(911)
      end do
   end if

   ! Read infile record
   call zaiord(hy2d,ip,.false.,hmina,hmaxa,911)

   call to_proj(hy2d,regu2d)
   where (regubathy<indepth(k)) regu2d=undef

   call handle_err(NF90_PUT_VAR(ncid,var_id  ,regu2d,start=(/1,1,k/)))
   call zaiocl(911)
end do
   


! Define temperature variable
call handle_err(nf90_redef(ncid))
call handle_err(NF90_DEF_VAR(ncid,'salinity',NF90_Float,vars_3d,var_id))
call handle_err(NF90_PUT_ATT(ncid,var_id,'long_name','salinity'))
call handle_err(NF90_PUT_ATT(ncid,var_id,'units','psu'))
call handle_err(NF90_PUT_ATT(ncid,var_id,'standard_name','sea_water_salinity'))
call handle_err(NF90_PUT_ATT(ncid,var_id,'missing_value',undef))
call handle_err(nf90_enddef(ncid))

! Open Temperature infile
print *,'Processing salinity'
infile='saln_sig0_m'//cmonth//'.a'
do k=1,numdepths

   print *,'Processing depth level ',indepth(k)
   call zaiopf(trim(infile),'old',911)

   if (kindx(k)>1) then
      do k2=1,kindx(k)-1
         call zaiosk(911)
      end do
   end if

   ! Read infile record
   call zaiord(hy2d,ip,.false.,hmina,hmaxa,911)

   call to_proj(hy2d,regu2d)


   where (regubathy<indepth(k)) regu2d=undef

   call handle_err(NF90_PUT_VAR(ncid,var_id  ,regu2d,start=(/1,1,k/)))
   call zaiocl(911)
end do
call handle_err(NF90_CLOSE(ncid))
print *,'File '//trim(ncfile)//' created'


end program

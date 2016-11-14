!KAL -- takes norsex-type input files, and creates a climatology from these
!KAL -- Norsex files are ascii-type files from NERSC routines.
!KAL -- Location of NORSEX files specified in environment variable NORSEX_PATH
!KAL -- TODO: Make sure the same mod_toproj as in Hyc2proj is used
program norsexclim
use mod_parameters, only : undef
use mod_toproj
use mod_year_info
use netcdf
use m_handle_err
implicit none

integer, parameter :: easenx=304,easeny=448 ! Also se format 100 if you change this
character(len=*), parameter :: cver='0.6'

integer, dimension(easenx) :: line
real   , dimension(easenx,easeny) :: ssmilon, ssmilat
integer, dimension(easenx,easeny) :: ssmiip
integer :: k
real, dimension(easenx,easeny,12) :: cice
real, dimension(easenx,easeny)    :: cicetmp
integer                           :: cicecnt(12)
real, dimension(:,:), allocatable :: ofld
character(len=80) :: ncfile,vartitle,grid_mapping,coords_mapping,fname
character(len=  4) :: char_fac
character(len= 8) :: dateinfo
character(len=200) :: ns_path


integer :: xdimid,ydimid,ncid,tdimid,tdimid2
integer :: vars_2d(2), var_id, vars_3d(3),vars_2d_clim(2)
integer :: ilon,jlat,jday
integer :: baseyear,iyear,imonth
logical :: fileex

real, allocatable, dimension(:) :: tmpx, tmpy


! Get norsex grid
call norsex_getgrid(ssmilon,ssmilat,easenx,easeny)

print *,minval(ssmilon),maxval(ssmilon)
print *,minval(ssmilat),maxval(ssmilat)

! Get base path of NORSEX data
call getenv('NORSEX_PATH',ns_path)
if (trim(ns_path)=='') then
   print *,'Provide path to NOrSEX top directory in env variable NORSEX_PATH'
   stop
end if
ns_path=trim(ns_path)//'/'


! Test read of one ascii file
open(10,file=trim(ns_path)//'ICE_CONC/TOTAL/1978/7811T.ASC',form='formatted',status='old')
do k=1,easeny
   read(10,100) line
   cicetmp(:,k)=line
end do
close(10)
ssmiip=1
where (cicetmp(:,:)>100.) ssmiip=0 ! SSMI mask


! Loop through years and months
baseyear=1900
imonth=11 ! First month in 1978 to use
fileex=.true.
cice=0.
cicecnt=0
do iyear=78,104

   do while (imonth<=12 .and. fileex)

   
      ! Test read of one ascii file
      write(fname ,'("ICE_CONC/TOTAL/",i4.4,"/",i2.2,i2.2,"T.ASC")') iyear+baseyear,iyear-(iyear/100)*100,imonth

      inquire(exist=fileex,file=trim(ns_path)//trim(fname))
      if (.not.fileex) then
         write(fname,'("ICE_CONC/TOTAL/",i4.4,"/",i2.2,i2.2,"T.asc")') iyear+baseyear,iyear-(iyear/100)*100,imonth
         inquire(exist=fileex,file=trim(ns_path)//trim(fname))
      end if

      print *,trim(fname),fileex


      if (fileex) then
         !open(10,file='ICE_CONC/TOTAL/1978/7811T.ASC',form='formatted',status='old')
         open(10,file=trim(ns_path)//trim(fname),form='formatted',status='old')
         do k=1,easeny
            read(10,100) line
            cicetmp(:,k)=line
         end do
         close(10)
         where(cicetmp>100.) ssmiip=0
         cice(:,:,imonth)=cice(:,:,imonth)+cicetmp
         cicecnt(imonth)=cicecnt(imonth)+1
      end if

      imonth=imonth+1
   end do
   imonth=1
end do





! Initialize projection - NB: old version of the one in Hyc2proj
call get_proj()
call proj_ini(undef,ssmilon,ssmilat,ssmiip,easenx,easeny)

! For dumping files
allocate(tmpx(nlons))
allocate(tmpy(nlats))
allocate(ofld(nlons,nlats))
ncfile='norsex_cice_climatology_arctic.nc'
if (NF90_CREATE(ncfile,NF90_CLOBBER,ncid) /= NF90_NOERR) then
   print *,'An error occured when opening the netcdf file'
   stop '(norsexclim)'
end if

! Define dimensions -- idm,jdm,kdm
! Following MERCATOR name definitions for netcdf var name
if (trim(cprojection)=='regular') then
   call handle_err(NF90_DEF_DIM(ncid,'longitude',nlons,xdimid))
   call handle_err(NF90_DEF_DIM(ncid,'latitude',nlats,ydimid))
else if (trim(cprojection)=='polar_stereographic' .or. &
         trim(cprojection)=='native') then
   call handle_err(NF90_DEF_DIM(ncid,'x',nlons,xdimid))
   call handle_err(NF90_DEF_DIM(ncid,'y',nlats,ydimid))
else
   print *,'Unknown projection '//trim(cprojection)
end if
call handle_err(NF90_DEF_DIM(ncid,'time',12,tdimid))
call handle_err(NF90_DEF_DIM(ncid,'nv',2,tdimid2)) ! for clim bounds
vars_2d=(/xdimid,ydimid/)
vars_2d_clim=(/tdimid,tdimid2/)
vars_3d=(/xdimid,ydimid,tdimid/)

! Define variables -- put attributes
!Global attributes
grid_mapping=''
coords_mapping=''
call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'title',        &
      'GDEM climatology - on Arctic MERSEA grid'))
call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'institution',        &
   'NERSC, Thormoehlens gate 47, N-5006 Bergen, Norway'))
call date_and_time(date=dateinfo)
call handle_err(NF90_PUT_ATT(ncid,NF90_GLOBAL,'history', &
   dateinfo(1:8)//':Created by program norsexclim, version '//cver))
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
   call handle_err(NF90_DEF_VAR(ncid,'longitude',NF90_Float,xdimid,var_id))
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
   call handle_err(NF90_DEF_VAR(ncid,trim(vartitle),NF90_Float,xdimid,var_id))
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
   call handle_err(NF90_DEF_VAR(ncid,trim(vartitle),NF90_Float,xdimid,var_id))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'axis','X'))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'standard_name','projection_x_coordinate'))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'units',adjustl(trim(char_fac))//' km'))
   do ilon=1,nlons
      tmpx(ilon)=firstlon+(ilon-1)*dlon
   end do
   call handle_err(nf90_enddef(ncid))
   call handle_err(NF90_PUT_VAR(ncid,var_id,tmpx))
else 
   print *,'undef proj'
   stop
end if






! Latitude
call handle_err(nf90_redef(ncid))
if (trim(cprojection)=='regular') then
   vartitle= 'latitude'
   call handle_err(NF90_DEF_VAR(ncid,trim(vartitle),NF90_Float,ydimid,var_id))
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
   call handle_err(NF90_DEF_VAR(ncid,trim(vartitle),NF90_Float,ydimid,var_id))
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
   call handle_err(NF90_DEF_VAR(ncid,trim(vartitle),NF90_Float,ydimid,var_id))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'standard_name','projection_y_coordinate'))
   call handle_err(NF90_PUT_ATT(ncid,var_id,'axis','Y'))
   do jlat=1,nlats
      tmpy(jlat)=firstlat+(jlat-1)*dlat
   end do
   call handle_err(nf90_enddef(ncid))
   call handle_err(NF90_PUT_VAR(ncid,var_id,real(tmpy,4)))
else
   print *,'undef proj'
   stop
end if


! KAL -- time dimension
vartitle='time'
call handle_err(nf90_redef(ncid))
call handle_err(NF90_DEF_VAR(ncid,trim(vartitle),NF90_Float,tdimid,var_id))
call handle_err(NF90_PUT_ATT(ncid,var_id,'climatology','climatology_bounds'))
call handle_err(NF90_PUT_ATT(ncid,var_id,'units','days since 1979-01-01'))
call handle_err(nf90_enddef(ncid))
do k=1,12
   jday=datetojulian(1979,k,15,1979,1,1)
   call handle_err(NF90_PUT_VAR(ncid,var_id,jday,start=(/k/)))
end do

! CF - compliant spec of climatology
vartitle='climatology_bounds'
call handle_err(nf90_redef(ncid))
call handle_err(NF90_DEF_VAR(ncid,trim(vartitle),NF90_Float,vars_2d_clim,var_id))
call handle_err(NF90_PUT_ATT(ncid,var_id,'units','days since 1979-01-01'))
call handle_err(nf90_enddef(ncid))
do k=1,12
   jday=datetojulian(1979,k,1,1979,1,1)
   call handle_err(NF90_PUT_VAR(ncid,var_id,jday,start=(/k,1/)))
   if (k<12) then
      jday=datetojulian(2004,k+1,1,1979,1,1)
   else
      jday=datetojulian(2004,k,31,1979,1,1)
   end if
   call handle_err(NF90_PUT_VAR(ncid,var_id,jday,start=(/k,2/)))
end do





! Define temperature variable
call handle_err(nf90_redef(ncid))
!print *,'def cice'
call handle_err(NF90_DEF_VAR(ncid,'fice',NF90_Float,vars_3d,var_id))
call handle_err(NF90_PUT_ATT(ncid,var_id,'long_name','Sea Ice Concentration'))
call handle_err(NF90_PUT_ATT(ncid,var_id,'units',''))
call handle_err(NF90_PUT_ATT(ncid,var_id,'standard_name','sea_ice_area_fraction'))
call handle_err(NF90_PUT_ATT(ncid,var_id,'missing_value',undef))
call handle_err(NF90_PUT_ATT(ncid,var_id,'_FillValue',real(undef,kind=4)))
if (trim(grid_mapping)/='') call handle_err(NF90_PUT_ATT(ncid,var_id,'grid_mapping',trim(grid_mapping)))
if (trim(coords_mapping)/='') call handle_err(NF90_PUT_ATT(ncid,var_id,'coordinates',trim(coords_mapping)))
call handle_err(NF90_PUT_ATT(ncid,var_id,'cell_methods','time: mean within years time: mean over years'))
call handle_err(nf90_enddef(ncid))

do k=1,12
   cice(:,:,k)=cice(:,:,k)/(cicecnt(k)*100.)
   print *,minval(cice(:,:,k)),maxval(cice(:,:,k))
   call to_proj(cice(:,:,k),ofld,easenx,easeny)
   print *,minval(ofld),maxval(ofld)
   call handle_err(NF90_PUT_VAR(ncid,var_id  ,ofld,start=(/1,1,k/)))
end do


ilon=NF90_CLOSE(ncid)


100 format(304i4)
end program

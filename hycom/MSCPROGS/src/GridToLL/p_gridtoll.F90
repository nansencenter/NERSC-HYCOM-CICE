! Routine to interpolate from model grid (confmap) to longitude-latitude.
! uses netcdf files created by m2nc
!
! Final files are suitable for use with the gmt package 

program p_gridtoll
use mod_xc
use mod_confmap
use netcdf
use m_handle_err
implicit none


integer :: idm2, jdm2
integer :: nlon, nlat
real :: firstlon, lastlon, dxlon
real :: firstlat, lastlat, dxlat
real, allocatable :: gridlon(:), gridlat(:)
real, allocatable, dimension(:,:) :: fldin,fldout
character(len=10) :: cvar

character(len=80) :: ctmp,ctmp2
integer :: i,j, ipiv, jpiv
integer :: ncid,stat
integer :: idimid, jdimid, varid, var_ndim
integer :: var_dimids(nf90_max_dims)
integer :: londimid, latdimid, varidx, varidy, varidz
real :: lat_n, lon_n
integer :: ind

#if defined(IARGC)
integer*4, external :: iargc
#endif


if (iargc()==3) then
   call getarg(1,ctmp); 
   ind=index(ctmp,':'); read(ctmp(1:ind-1),*) firstlon ; ctmp=ctmp(ind+1:len_trim(ctmp))
   ind=index(ctmp,':'); read(ctmp(1:ind-1),*) dxlon    ; ctmp=ctmp(ind+1:len_trim(ctmp))
   read(ctmp,*) lastlon 
   print *,firstlon,dxlon,lastlon

   call getarg(2,ctmp); 
   ind=index(ctmp,':'); read(ctmp(1:ind-1),*) firstlat ; ctmp=ctmp(ind+1:len_trim(ctmp))
   ind=index(ctmp,':'); read(ctmp(1:ind-1),*) dxlat    ; ctmp=ctmp(ind+1:len_trim(ctmp))
   read(ctmp,*) lastlat 
   print *,firstlat,dxlat,lastlat
   
   !firstlon
   call getarg(3,cvar) 
else
   print *,'Routine reads from a netcdf file created with m2nc and interpolates'
   print *,'The data onto a rectilinear, uniformly spaced longitude/latitude grid'
   print *
   print *,'Routine should be called like this:'
   print *,'   gridtoll lon_spec lat_spec variable_name'
   print *,'lon_spec and lat_spec is a text string of the type '
   print *,'  first:increment:last '
   print *,'Denoting the first value of the regular grid, spacing, and the last value. ' 
   print *
   print *,'Example:'
   print *,'   gridtoll -98:.1:20 10:.1:60 temp01'
   print *,'Will read temp01 from "tmp1.nc" and interpolate it on the grid'
   stop
end if
call xcspmd()
call initconfmap(idm,jdm)

! Set up output grid
nlon=(lastlon-firstlon)/dxlon + 1
nlat=(lastlat-firstlat)/dxlat + 1
print *,nlon,nlat
allocate(gridlon(nlon))
allocate(gridlat(nlat))
allocate(fldout(nlon,nlat))

do i=1,nlon
   gridlon(i)=firstlon+(i-1)*dxlon
end do
do j=1,nlat
   gridlat(j)=firstlat+(j-1)*dxlat
end do

! Read input fields
stat = nf90_open('tmp1.nc',nf90_nowrite,ncid) 

! Get dimensions
stat = nf90_inq_dimid(ncid,'idim',idimid)
stat = nf90_inquire_dimension(ncid,idimid,len=idm2)
stat = nf90_inq_dimid(ncid,'jdim',jdimid) 
stat = nf90_inquire_dimension(ncid,jdimid,len=jdm2)
print *,'idm2, jdm2',idm2,jdm2
if (idm2/=idm .or. jdm2/=jdm) then
   print *,'Grid Size mismatch between netcdf file and grid'
   stop '(gritoll)'
end if

! Get variable
allocate(fldin(idm2,jdm2))
call handle_err(nf90_inq_varid(ncid,trim(cvar),varid) )
call handle_err(nf90_inquire_variable(ncid,varid,ndims=var_ndim,dimids=var_dimids))
call handle_err(nf90_get_var(ncid, varid, fldin, start=(/1,1,1/)))
call handle_err(nf90_close(ncid))

do j=1,nlat
do i=1,nlon
   ! corresponding pivot points
   call oldtonew(gridlat(j),gridlon(i),lat_n,lon_n)
   call pivotp(lon_n,lat_n,ipiv,jpiv)

   ! NB - no bilinear interpolation here
   if (ipiv<=idm .and. ipiv>1 .and.  jpiv<=jdm .and. jpiv>1 ) then 
       fldout(i,j)=fldin(ipiv,jpiv)
   else
      fldout(i,j)=-1e14
   end if
end do
end do


stat = nf90_create('tmp2.nc',nf90_clobber,ncid) 
stat = nf90_def_dim(ncid, 'x', nlon, londimid)
stat = nf90_def_dim(ncid, 'y', nlat, latdimid)
stat = nf90_def_var(ncid, 'x', nf90_float, (/londimid/), varidx)
stat = nf90_def_var(ncid, 'y', nf90_float, (/latdimid/), varidy)
stat = nf90_def_var(ncid, 'z', nf90_float, (/londimid, latdimid/), varidz)
stat = nf90_enddef(ncid)
stat = nf90_put_var(ncid, varidx, gridlon) ; print *,stat
stat = nf90_put_var(ncid, varidy, gridlat)
stat = nf90_put_var(ncid, varidz, fldout)
stat = nf90_close(ncid)


end program









! Routine to find nearest point in a model grid (confmap) to longitude-latitude.
program p_findnearest
use mod_xc
use mod_za
use mod_grid, only:get_grid,plon,plat
use mod_confmap
use m_spherdist
implicit none


integer :: idm2, jdm2
integer :: nlon, nlat
real :: tlon, lastlon, dxlon
real :: tlat, lastlat, dxlat
real, allocatable :: gridlon(:), gridlat(:)
real, allocatable, dimension(:,:) :: fldin,fldout
character(len=10) :: cvar

character(len=80) :: ctmp,ctmp2
integer :: i,j, ipiv, jpiv
integer :: ncid,stat
integer :: idimid, jdimid, varid, var_ndim
integer :: londimid, latdimid, varidx, varidy, varidz
real :: lat_n, lon_n
integer :: ind

#if defined(IARGC)
integer*4, external :: iargc
#endif
if (iargc()==2) then
   call getarg(1,ctmp); 
   read(ctmp,*) tlon 
   call getarg(2,ctmp); 
   read(ctmp,*) tlat
else
   print *,'Routine need a lon and lat and return the nearest model gridpoint'
   print *,'Example:'
   print *,'   findnearest -98.1 10.1' 
   stop
end if
call xcspmd()
call initconfmap(idm,jdm)
call zaiost()
call get_grid

call oldtonew(tlat,tlon,lat_n,lon_n)
call pivotp(lon_n,lat_n,ipiv,jpiv)
print *, 'Point ',ipiv,jpiv,'dist',spherdist(tlon,tlat,plon(ipiv,jpiv),plat(ipiv,jpiv))
print *, 'Point ',ipiv+1,jpiv,'dist',spherdist(tlon,tlat,plon(ipiv+1,jpiv),plat(ipiv+1,jpiv))
print *, 'Point ',ipiv+1,jpiv+1,'dist',spherdist(tlon,tlat,plon(ipiv+1,jpiv+1),plat(ipiv+1,jpiv+1))
print *, 'Point ',ipiv,jpiv+1,'dist',spherdist(tlon,tlat,plon(ipiv,jpiv+1),plat(ipiv,jpiv+1))
end program









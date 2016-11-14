module m_ll2gind
contains
subroutine ll2gind(lon_o,lat_o,x,y)
! This subroutine computes the pivot point of each of the observations
! in the temporary array tmpobs of type observation. The pivot point
! is the biggest i and the biggest j, (i,j) is the computation points/
! the grid, that is less than the position of the observation.

   use mod_confmap
   use m_oldtonew
   implicit none

   real, intent(in) ::  lon_o,lat_o
   real, intent(out) :: x,y

   real tmptan,tmp, lon, lat
   real lontmp

   call oldtonew(lat_o,lon_o,lat,lon)

! fix for wrap-around
! western limit in new coordinates can be east of eastern limit (sigh)....
! in that case di is < 0
   if (lon < wlim .and. di > 0. ) then
      lontmp=lon+360.0
   elseif (lon > wlim .and. di < 0. ) then
         lontmp=lon-360.0
!more fixes ...
   elseif ((lon-wlim)>360.) then
      lontmp=lon
      do while(lontmp-wlim>360.)
         lontmp=lontmp-360.
      end do
   else
      lontmp=lon
   endif

   !print '(a,2i5,4f10.2)','pivotp:',ipiv,jpiv,lontmp,lon,wlim,di
   x=(lontmp-wlim)/di+1

   if (mercator) then
      if (abs(lat) < 89.999) then
         tmptan=tan(0.5*rad*lat+0.25*pi_1)
         y=(log(tmptan)-slim*rad)/(rad*dj) +1
      else
         y=-999
      endif 
   else
      y=(lat-slim)/dj+1
   endif


! Inverse transformation to check pivot point jpiv
   tmp=slim+(y-1)*dj
   tmp=(2.*atan(exp(tmp*rad))-pi_1*.5)*deg

   !print *,'ll2gind',lon_o,lat_o,tmp

end subroutine ll2gind
end module m_ll2gind

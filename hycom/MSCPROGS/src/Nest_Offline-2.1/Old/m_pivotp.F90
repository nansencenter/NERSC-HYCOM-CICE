module m_pivotp
contains
subroutine pivotp(lon,lat,ipiv,jpiv)
! This subroutine computes the pivot point of each of the observations
! in the temporary array tmpobs of type observation. The pivot point
! is the biggest i and the biggest j, (i,j) is the computation points/
! the grid, that is less than the position of the observation.

   use mod_confmap
   implicit none

   real, intent(in) ::  lon,lat
   integer, intent(out) :: ipiv,jpiv

   real tmptan,tmp
   real lontmp

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
   ipiv=int((lontmp-wlim)/di)+1

   if (mercator) then
      if (abs(lat) < 89.999) then
         tmptan=tan(0.5*rad*lat+0.25*pi_1)
         jpiv=int( (log(tmptan)-slim*rad)/(rad*dj) ) +1
      else
         jpiv=-999
      endif 
   else
      jpiv=int((lat-slim)/dj)+1
   endif


! Inverse transformation to check pivot point jpiv
!   tmp=slim+(jpiv-1)*dj
!   tmp=(2.*atan(exp(tmp*rad))-pi_1*.5)*deg

end subroutine pivotp
end module m_pivotp
